#!/usr/bin/env python3
"""Build target/caller-specific consensus intervals with biological support."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merged_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{number}: expected at least three BED fields")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: invalid interval")
            by_chrom[fields[0]].append((start, end))
    result: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in by_chrom.items():
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        result[chrom] = [(start, end) for start, end in merged]
    return result


def consensus(files: list[tuple[str, Path]], minimum: int) -> list[tuple[str, int, int, int]]:
    events: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for sample, path in files:
        for chrom, intervals in merged_intervals(path).items():
            for start, end in intervals:
                events[chrom].append((start, 1, sample))
                events[chrom].append((end, -1, sample))
    result: list[tuple[str, int, int, int]] = []
    for chrom in sorted(events):
        by_position: dict[int, list[tuple[int, str]]] = defaultdict(list)
        for position, direction, sample in events[chrom]:
            by_position[position].append((direction, sample))
        active: set[str] = set()
        previous: int | None = None
        supported: list[list[int]] = []
        for position in sorted(by_position):
            if previous is not None and position > previous and len(active) >= minimum:
                if supported and supported[-1][1] == previous:
                    supported[-1][1] = position
                    supported[-1][2] = max(supported[-1][2], len(active))
                else:
                    supported.append([previous, position, len(active)])
            # End events precede start events at a boundary; zero-width overlap is not support.
            for direction, sample in sorted(by_position[position]):
                if direction < 0:
                    active.discard(sample)
            for direction, sample in sorted(by_position[position]):
                if direction > 0:
                    active.add(sample)
            previous = position
        result.extend((chrom, start, end, support) for start, end, support in supported)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--cohort-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--minimum-support", type=int, default=2)
    parser.add_argument("--allow-single", action="store_true")
    parser.add_argument("--require-all", action="store_true")
    args = parser.parse_args()
    samples = read_manifest(args.sample_manifest)
    cohorts = read_manifest(args.cohort_manifest)
    failures = 0
    summary: list[dict[str, object]] = []
    for cohort in cohorts:
        cohort_id = cohort["cohort_id"]
        members = [row for row in samples if row["cohort_id"] == cohort_id and row["is_control"] == "FALSE"]
        caller, peak_class = cohort["primary_peak_caller"], cohort["primary_peak_class"]
        outdir = args.output_root / cohort_id / caller / peak_class
        outdir.mkdir(parents=True, exist_ok=True)
        minimum = args.minimum_support
        if len(members) < minimum:
            if args.allow_single and len(members) == 1:
                minimum = 1
            else:
                reason = f"{len(members)} biological samples; minimum support is {minimum}"
                (outdir / "SKIPPED.json").write_text(json.dumps({"status": "SKIPPED", "reason": reason}, indent=2) + "\n")
                summary.append({"cohort_id": cohort_id, "status": "SKIPPED", "total_samples": len(members),
                                "successful_peak_samples": 0, "excluded_samples": 0,
                                "regions": 0, "reason": reason})
                failures += int(args.require_all)
                continue
        peak_files: list[tuple[str, Path]] = []
        excluded: list[dict[str, str]] = []
        inconsistent: list[str] = []
        for member in members:
            sample_key = member["sample_key"]
            sample_root = args.output_root.parent / "per_sample" / sample_key
            peak = sample_root / caller / f"{sample_key}.{caller}.{peak_class}.bed"
            metadata = sample_root / "peakcall_metadata.tsv"
            if not metadata.is_file():
                inconsistent.append(f"{sample_key}: missing peakcall metadata")
                continue
            rows = read_manifest(metadata)
            if len(rows) != 1:
                inconsistent.append(f"{sample_key}: invalid peakcall metadata")
                continue
            status = rows[0].get("status", "ERROR")
            reason = rows[0].get("reason", "unspecified")
            if status not in {"SUCCESS", "EMPTY", "ERROR"}:
                inconsistent.append(f"{sample_key}: invalid peakcall status {status!r}")
                continue
            if (rows[0].get("primary_caller") != caller or
                    rows[0].get("primary_class") != peak_class):
                inconsistent.append(f"{sample_key}: primary caller/class metadata disagrees with cohort")
                continue
            if status == "SUCCESS":
                if not peak.is_file() or peak.stat().st_size == 0:
                    inconsistent.append(f"{sample_key}: SUCCESS but primary peak file is missing/empty")
                else:
                    peak_files.append((sample_key, peak))
            else:
                excluded.append({"sample_key": sample_key, "status": status, "reason": reason})
        with (outdir / "excluded_peak_samples.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample_key", "status", "reason"],
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(excluded)
        if inconsistent:
            reason = "inconsistent primary peak outputs: " + ", ".join(inconsistent)
            (outdir / "FAILED.json").write_text(json.dumps({"status": "FAILED", "reason": reason}, indent=2) + "\n")
            summary.append({"cohort_id": cohort_id, "status": "FAILED", "total_samples": len(members),
                            "successful_peak_samples": len(peak_files), "excluded_samples": len(excluded),
                            "regions": 0, "reason": reason})
            failures += 1
            continue
        if len(peak_files) < minimum:
            if args.allow_single and len(peak_files) == 1:
                minimum = 1
            else:
                reason = (f"{len(peak_files)} successful primary peak samples after "
                          f"{len(excluded)} exclusions; minimum support is {minimum}")
                status = "FAILED" if args.require_all else "SKIPPED"
                (outdir / f"{status}.json").write_text(
                    json.dumps({"status": status, "reason": reason}, indent=2) + "\n"
                )
                summary.append({"cohort_id": cohort_id, "status": status, "total_samples": len(members),
                                "successful_peak_samples": len(peak_files), "excluded_samples": len(excluded),
                                "regions": 0, "reason": reason})
                failures += int(status == "FAILED")
                continue
        regions = consensus(peak_files, minimum)
        if not regions:
            reason = "no intervals meet biological support"
            status = "FAILED" if args.require_all else "SKIPPED"
            (outdir / f"{status}.json").write_text(json.dumps({"status": status, "reason": reason}, indent=2) + "\n")
            summary.append({"cohort_id": cohort_id, "status": status, "total_samples": len(members),
                            "successful_peak_samples": len(peak_files), "excluded_samples": len(excluded),
                            "regions": 0, "reason": reason})
            failures += int(status == "FAILED")
            continue
        bed = outdir / f"{cohort_id}.{caller}.{peak_class}.support-ge{minimum}.consensus.bed"
        with bed.open("w", encoding="utf-8", newline="\n") as handle:
            for index, (chrom, start, end, support) in enumerate(regions, 1):
                handle.write(f"{chrom}\t{start}\t{end}\t{cohort_id}.region{index:07d}\t{support}\n")
        with gzip.open(outdir / "consensus_support.tsv.gz", "wt", encoding="utf-8", newline="\n") as handle:
            handle.write("region_id\tchrom\tstart\tend\tmaximum_support\n")
            for index, (chrom, start, end, support) in enumerate(regions, 1):
                handle.write(f"{cohort_id}.region{index:07d}\t{chrom}\t{start}\t{end}\t{support}\n")
        summary.append({"cohort_id": cohort_id, "status": "SUCCESS", "total_samples": len(members),
                        "successful_peak_samples": len(peak_files), "excluded_samples": len(excluded),
                        "regions": len(regions), "reason": "."})
    summary_path = args.output_root / "consensus_status.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cohort_id", "status", "total_samples",
                                "successful_peak_samples", "excluded_samples", "regions", "reason"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print(f"CONSENSUS ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
