#!/usr/bin/env python3
"""Resolve CUT workflow bigWigs into the shared metagene track-manifest contract."""

from __future__ import annotations

import argparse
import csv
import shlex
import sys
from pathlib import Path


def read_config(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        key, value = raw.split("=", 1)
        parts = shlex.split(value)
        result[key] = parts[0] if parts else ""
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--resolved-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--track-family", choices=["auto", "cpm", "spikein"], default="auto")
    parser.add_argument("--allow-cpm-fallback", action="store_true")
    parser.add_argument("--include-controls", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        config = read_config(args.resolved_config)
        samples = read_tsv(args.sample_manifest)
        if not samples:
            raise ValueError("sample manifest has no rows")
        spike_mode = config.get("SPIKEIN_MODE", "none")
        family = args.track_family
        if family == "auto":
            family = "spikein" if spike_mode != "none" else "cpm"
        if family == "spikein" and spike_mode == "none":
            raise ValueError("spike-in track family requested but SPIKEIN_MODE=none")

        scaling_path = args.output_dir / "06_qc/spikein/spikein_scaling.tsv"
        scaling = {row["sample_key"]: row for row in read_tsv(scaling_path)} if scaling_path.is_file() else {}
        rows: list[dict[str, str]] = []
        for sample in samples:
            if sample["is_control"] == "TRUE" and not args.include_controls:
                continue
            sample_key = sample["sample_key"]
            genome = sample["genome"]
            suffix = genome.upper().replace(".", "_").replace("-", "_")
            chrom_sizes = config.get(f"CHROM_SIZES_{suffix}", "")
            if not chrom_sizes:
                raise ValueError(f"missing chromosome sizes for {genome}")
            normalization = "CPM"
            detail = "upstream_fragment_or_read_CPM"
            bigwig = args.output_dir / f"04_tracks/cpm/{sample_key}.CPM.bw"
            if family == "spikein":
                cohort = sample["cohort_id"]
                candidate = args.output_dir / f"04_tracks/spikein/{cohort}/{sample_key}.SpikeInScaled.{spike_mode}.bw"
                scale_row = scaling.get(sample_key)
                passed = scale_row is not None and scale_row.get("status") != "FAILED"
                if candidate.is_file() and passed:
                    bigwig = candidate
                    normalization = "spikein"
                    detail = f"reference={spike_mode};host_scale_factor={scale_row.get('host_scale_factor', '.')}"
                elif not args.allow_cpm_fallback:
                    raise ValueError(f"{sample_key}: spike-in track unavailable or failed and CPM fallback is disabled")
                else:
                    detail = f"CPM_fallback_from_failed_or_missing_{spike_mode}_spikein"
            if not bigwig.is_file():
                raise ValueError(f"{sample_key}: selected bigWig does not exist: {bigwig}")
            rows.append({
                "sample_id": sample_key, "assay": sample["assay_profile"], "genome": genome,
                "bigwig": str(bigwig.resolve()), "normalization": normalization,
                "normalization_detail": detail, "blacklist": str(Path(sample["blacklist"]).resolve()),
                "chrom_sizes": str(Path(chrom_sizes).resolve()),
            })
        if not rows:
            raise ValueError("no samples remain after control filtering")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fields = ["sample_id", "assay", "genome", "bigwig", "normalization", "normalization_detail", "blacklist", "chrom_sizes"]
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        print(f"Resolved {len(rows)} metagene tracks using {family} preference")
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"METAGENE TRACK ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
