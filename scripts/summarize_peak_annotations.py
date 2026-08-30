#!/usr/bin/env python3
"""Annotate per-sample and consensus peaks and summarize genomic feature composition."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


DEFAULT_CATEGORIES = (
    "promoter", "enhancer", "exon", "intron", "gene_end",
    "other_regulatory", "intergenic", "unclassified",
)
COLORS = {
    "promoter": "#E41A1C",
    "enhancer": "#FF7F00",
    "exon": "#377EB8",
    "intron": "#4DAF4A",
    "gene_end": "#984EA3",
    "other_regulatory": "#A65628",
    "intergenic": "#999999",
    "unclassified": "#000000",
}


@dataclass(frozen=True)
class Feature:
    chrom: str
    start: int
    end: int
    category: str
    feature_id: str
    gene_id: str = "."
    gene_name: str = "."
    strand: str = "."


class IntervalIndex:
    def __init__(self, features: Iterable[Feature]):
        self.features: dict[str, list[Feature]] = defaultdict(list)
        self.starts: dict[str, list[int]] = {}
        self.prefix_max_end: dict[str, list[int]] = {}
        for feature in features:
            if feature.end > feature.start:
                self.features[feature.chrom].append(feature)
        for chrom, values in self.features.items():
            values.sort(key=lambda item: (item.start, item.end, item.category, item.feature_id))
            self.starts[chrom] = [item.start for item in values]
            maximum = 0
            prefix: list[int] = []
            for item in values:
                maximum = max(maximum, item.end)
                prefix.append(maximum)
            self.prefix_max_end[chrom] = prefix

    def overlaps(self, chrom: str, start: int, end: int) -> list[tuple[Feature, int]]:
        values = self.features.get(chrom, [])
        if not values:
            return []
        index = bisect.bisect_left(self.starts[chrom], end) - 1
        result: list[tuple[Feature, int]] = []
        while index >= 0 and self.prefix_max_end[chrom][index] > start:
            feature = values[index]
            overlap = min(end, feature.end) - max(start, feature.start)
            if overlap > 0:
                result.append((feature, overlap))
            index -= 1
        return result


def attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+) "([^"]*)"', text)}


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[list[int]] = []
    for start, end in sorted(intervals):
        if not result or start > result[-1][1]:
            result.append([start, end])
        elif end > result[-1][1]:
            result[-1][1] = end
    return [(start, end) for start, end in result]


def parse_gtf(path: Path, promoter_upstream: int, promoter_downstream: int,
              gene_end_window: int) -> tuple[list[Feature], list[dict[str, object]]]:
    genes: dict[str, dict[str, object]] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "exon"}:
                continue
            chrom, kind, start, end, strand = fields[0], fields[2], int(fields[3]) - 1, int(fields[4]), fields[6]
            meta = attributes(fields[8])
            gene_id = meta.get("gene_id", ".")
            if gene_id == ".":
                continue
            if kind == "gene":
                genes[gene_id] = {
                    "chrom": chrom, "start": start, "end": end, "strand": strand,
                    "gene_id": gene_id, "gene_name": meta.get("gene_name", gene_id),
                }
            else:
                exons[gene_id].append((start, end))

    features: list[Feature] = []
    gene_rows: list[dict[str, object]] = []
    for gene_id, gene in sorted(genes.items()):
        chrom, start, end = str(gene["chrom"]), int(gene["start"]), int(gene["end"])
        strand, gene_name = str(gene["strand"]), str(gene["gene_name"])
        gene_rows.append(gene)
        if strand == "-":
            tss = end
            promoter_start, promoter_end = max(0, tss - promoter_downstream), tss + promoter_upstream
            tes_start, tes_end = max(0, start - gene_end_window), start
        else:
            tss = start
            promoter_start, promoter_end = max(0, tss - promoter_upstream), tss + promoter_downstream
            tes_start, tes_end = end, end + gene_end_window
        features.append(Feature(chrom, promoter_start, promoter_end, "promoter", f"{gene_id}:promoter",
                                gene_id, gene_name, strand))
        if tes_end > tes_start:
            features.append(Feature(chrom, tes_start, tes_end, "gene_end", f"{gene_id}:gene_end",
                                    gene_id, gene_name, strand))
        merged_exons = merge_intervals(exons.get(gene_id, []))
        for number, (exon_start, exon_end) in enumerate(merged_exons, 1):
            features.append(Feature(chrom, exon_start, exon_end, "exon", f"{gene_id}:exon:{number}",
                                    gene_id, gene_name, strand))
        cursor = start
        intron_number = 0
        for exon_start, exon_end in merged_exons:
            exon_start, exon_end = max(start, exon_start), min(end, exon_end)
            if exon_start > cursor:
                intron_number += 1
                features.append(Feature(chrom, cursor, exon_start, "intron",
                                        f"{gene_id}:intron:{intron_number}", gene_id, gene_name, strand))
            cursor = max(cursor, exon_end)
        if cursor < end:
            intron_number += 1
            features.append(Feature(chrom, cursor, end, "intron", f"{gene_id}:intron:{intron_number}",
                                    gene_id, gene_name, strand))
    return features, gene_rows


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_ccre(path: Path | None) -> list[Feature]:
    if path is None:
        return []
    result: list[Feature] = []
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            label = "|".join(fields[3:]).lower()
            if "pls" in label or "promoter" in label:
                category = "promoter"
            elif any(token in label for token in ("dels", "pels", "enhancer", " els")):
                category = "enhancer"
            else:
                category = "other_regulatory"
            identifier = fields[3] if len(fields) > 3 and fields[3] else f"ccre:{number}"
            result.append(Feature(fields[0], start, end, category, identifier))
    return result


def nearest_gene(chrom: str, start: int, end: int,
                 genes: dict[str, tuple[list[int], list[dict[str, object]]]]) -> tuple[str, str, str]:
    midpoint = (start + end) // 2
    positions, values = genes.get(chrom, ([], []))
    if not positions:
        return ".", ".", "."
    insertion = bisect.bisect_left(positions, midpoint)
    candidates = values[max(0, insertion - 1):min(len(values), insertion + 1)]
    best: tuple[int, str, str, int] | None = None
    for gene in candidates:
        gstart, gend = int(gene["start"]), int(gene["end"])
        strand = str(gene["strand"])
        tss = gend if strand == "-" else gstart
        signed = midpoint - tss
        if strand == "-":
            signed = -signed
        candidate = (abs(signed), str(gene["gene_id"]), str(gene["gene_name"]), signed)
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return best[1], best[2], str(best[3])


def peak_sets(output: Path, include_consensus: bool) -> list[dict[str, str]]:
    sets: list[dict[str, str]] = []
    manifest = output / "00_metadata/sample_manifest.tsv"
    sample_meta: dict[str, dict[str, str]] = {}
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_meta[row["sample_key"]] = row
    for sample, meta in sorted(sample_meta.items()):
        if meta.get("is_control") == "TRUE":
            continue
        status_file = output / "05_peaks/per_sample" / sample / "caller_status.tsv"
        if not status_file.is_file():
            continue
        with status_file.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = Path(row.get("peak_file", ""))
                if row.get("status") != "SUCCESS" or not path.is_file():
                    continue
                sets.append({
                    "entity_type": "sample", "entity_id": sample, "sample_key": sample,
                    "cohort_id": meta.get("cohort_id", "."), "genome": meta["genome"],
                    "factor": meta.get("factor", "."), "condition": meta.get("condition", "."),
                    "replicate": meta.get("replicate", "."), "caller": row["caller"],
                    "peak_class": row["peak_class"], "peak_file": str(path),
                })
    if include_consensus:
        cohort_manifest = output / "00_metadata/cohort_manifest.tsv"
        with cohort_manifest.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                cohort = row["cohort_id"]
                caller, peak_class = row["primary_peak_caller"], row["primary_peak_class"]
                directory = output / "05_peaks/consensus" / cohort / caller / peak_class
                matches = sorted(directory.glob("*.consensus.bed")) if directory.is_dir() else []
                if not matches:
                    continue
                sets.append({
                    "entity_type": "consensus", "entity_id": cohort, "sample_key": ".",
                    "cohort_id": cohort, "genome": row["genome"], "factor": row.get("factor", "."),
                    "condition": row.get("conditions", "."), "replicate": ".", "caller": caller,
                    "peak_class": peak_class, "peak_file": str(matches[0]),
                })
    return sets


def read_chrom_sizes(path: Path) -> set[str]:
    with path.open(encoding="utf-8") as handle:
        return {line.split("\t", 1)[0] for line in handle if line.strip()}


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_composition(outdir: Path, rows: list[dict[str, object]], categories: tuple[str, ...], formats: list[str]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for peak-feature plots") from exc
    sample_rows = [row for row in rows if row["entity_type"] == "sample"]
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in sample_rows:
        groups[(str(row["caller"]), str(row["peak_class"]))].append(row)
    for (caller, peak_class), values in sorted(groups.items()):
        by_entity: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
        for row in values:
            by_entity[str(row["entity_id"])][str(row["category"])] = row
        labels = sorted(by_entity)
        if not labels:
            continue
        height = max(4.0, 0.42 * len(labels) + 2.5)
        figure, axes = plt.subplots(1, 3, figsize=(18, height), sharey=True)
        panels = (("count", "Peak count"), ("fraction", "Fraction of peaks"), ("bp_fraction", "Fraction of peak-covered bp"))
        for axis, (field, title) in zip(axes, panels):
            left = [0.0] * len(labels)
            for category in categories:
                numbers = [float(by_entity[label].get(category, {}).get(field, 0) or 0) for label in labels]
                axis.barh(labels, numbers, left=left, color=COLORS[category], label=category.replace("_", " "))
                left = [old + value for old, value in zip(left, numbers)]
            axis.set_title(title)
            axis.grid(axis="x", alpha=0.2)
        axes[1].set_xlim(0, 1)
        axes[2].set_xlim(0, 1)
        handles, legend_labels = axes[-1].get_legend_handles_labels()
        figure.legend(handles, legend_labels, loc="lower center", ncol=4, frameon=False)
        figure.suptitle(f"Peak genomic-feature composition: {caller} {peak_class}")
        figure.tight_layout(rect=(0, 0.08, 1, 0.96))
        stem = outdir / f"peak_feature_composition.{caller}.{peak_class}"
        for extension in formats:
            figure.savefig(Path(f"{stem}.{extension}"), dpi=200, bbox_inches="tight")
        plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path, required=True)
    parser.add_argument("--gtf", action="append", default=[], metavar="GENOME=PATH")
    parser.add_argument("--chrom-sizes", action="append", default=[], metavar="GENOME=PATH")
    parser.add_argument("--ccre", action="append", default=[], metavar="GENOME=PATH")
    parser.add_argument("--promoter-upstream", type=int, default=2000)
    parser.add_argument("--promoter-downstream", type=int, default=500)
    parser.add_argument("--gene-end-window", type=int, default=2000)
    parser.add_argument("--precedence", default=",".join(DEFAULT_CATEGORIES))
    parser.add_argument("--plot-formats", default="png,pdf,svg")
    parser.add_argument("--include-consensus", action="store_true")
    parser.add_argument("--skip-plots", action="store_true",
                        help="write annotation tables without figures (test/recovery use)")
    args = parser.parse_args()
    try:
        def mapping(values: list[str]) -> dict[str, Path]:
            result: dict[str, Path] = {}
            for value in values:
                genome, separator, path = value.partition("=")
                if not separator or not genome or not path:
                    raise ValueError(f"expected GENOME=PATH, received {value!r}")
                result[genome] = Path(path)
            return result

        gtfs, chrom_files, ccres = mapping(args.gtf), mapping(args.chrom_sizes), mapping(args.ccre)
        categories = tuple(item.strip() for item in args.precedence.split(",") if item.strip())
        if set(categories) != set(DEFAULT_CATEGORIES) or len(categories) != len(DEFAULT_CATEGORIES):
            raise ValueError("precedence must contain each supported category exactly once")
        if min(args.promoter_upstream, args.promoter_downstream, args.gene_end_window) < 0:
            raise ValueError("annotation windows must be non-negative")
        output, outdir = args.output_dir.resolve(), args.annotation_dir.resolve()
        outdir.mkdir(parents=True, exist_ok=True)
        sets = peak_sets(output, args.include_consensus)
        required_genomes = {item["genome"] for item in sets}
        missing = sorted(genome for genome in required_genomes if genome not in gtfs or genome not in chrom_files)
        if missing:
            raise ValueError("missing GTF/chrom-sizes mapping for: " + ",".join(missing))

        indexes: dict[str, IntervalIndex] = {}
        gene_indexes: dict[str, dict[str, tuple[list[int], list[dict[str, object]]]]] = {}
        canonical: dict[str, set[str]] = {}
        enhancer_status: dict[str, str] = {}
        for genome in sorted(required_genomes):
            features, genes = parse_gtf(gtfs[genome], args.promoter_upstream,
                                        args.promoter_downstream, args.gene_end_window)
            ccre_path = ccres.get(genome)
            if ccre_path and ccre_path.is_file() and ccre_path.stat().st_size:
                features.extend(parse_ccre(ccre_path))
                enhancer_status[genome] = "evaluated"
            else:
                enhancer_status[genome] = "not_evaluated"
            indexes[genome] = IntervalIndex(features)
            per_chrom_unsorted: dict[str, list[dict[str, object]]] = defaultdict(list)
            for gene in genes:
                per_chrom_unsorted[str(gene["chrom"])].append(gene)
            per_chrom: dict[str, tuple[list[int], list[dict[str, object]]]] = {}
            for chrom, values in per_chrom_unsorted.items():
                ordered = sorted(values, key=lambda gene: (
                    int(gene["end"]) if str(gene["strand"]) == "-" else int(gene["start"]),
                    str(gene["gene_id"]),
                ))
                positions = [int(gene["end"]) if str(gene["strand"]) == "-" else int(gene["start"])
                             for gene in ordered]
                per_chrom[chrom] = (positions, ordered)
            gene_indexes[genome] = per_chrom
            canonical[genome] = read_chrom_sizes(chrom_files[genome])

        assignment_rows: list[dict[str, object]] = []
        overlap_rows: list[dict[str, object]] = []
        status_rows: list[dict[str, object]] = []
        aggregate: dict[tuple[str, ...], dict[str, float]] = defaultdict(lambda: defaultdict(float))
        precedence_rank = {category: index for index, category in enumerate(categories)}
        for item in sets:
            valid = invalid = unclassified = 0
            path = Path(item["peak_file"])
            with open_text(path) as handle:
                for line_number, raw in enumerate(handle, 1):
                    if not raw.strip() or raw.startswith(("#", "track", "browser")):
                        continue
                    fields = raw.rstrip("\n").split("\t")
                    try:
                        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                        if start < 0 or end <= start:
                            raise ValueError
                    except (IndexError, ValueError):
                        invalid += 1
                        continue
                    valid += 1
                    peak_id = fields[3] if len(fields) > 3 and fields[3] else f"{item['entity_id']}.peak{valid:07d}"
                    overlaps = indexes[item["genome"]].overlaps(chrom, start, end)
                    if chrom not in canonical[item["genome"]]:
                        selected_category = "unclassified"
                        selected: Feature | None = None
                        selected_overlap = 0
                    elif overlaps:
                        selected, selected_overlap = min(
                            overlaps,
                            key=lambda pair: (precedence_rank[pair[0].category], -pair[1], pair[0].feature_id),
                        )
                        selected_category = selected.category
                    else:
                        selected_category = "intergenic"
                        selected, selected_overlap = None, 0
                    if selected_category == "unclassified":
                        unclassified += 1
                    nearest_id, nearest_name, tss_distance = nearest_gene(
                        chrom, start, end, gene_indexes[item["genome"]]
                    )
                    overlapping_categories = sorted({feature.category for feature, _ in overlaps},
                                                    key=lambda value: precedence_rank[value])
                    assignment = dict(item)
                    assignment.update({
                        "peak_id": peak_id, "chrom": chrom, "start": start, "end": end,
                        "width": end - start, "primary_category": selected_category,
                        "overlapping_categories": ",".join(overlapping_categories) or ".",
                        "primary_feature_id": selected.feature_id if selected else ".",
                        "primary_gene_id": selected.gene_id if selected else ".",
                        "primary_gene_name": selected.gene_name if selected else ".",
                        "primary_overlap_bp": selected_overlap, "nearest_gene_id": nearest_id,
                        "nearest_gene_name": nearest_name, "nearest_tss_signed_distance": tss_distance,
                        "enhancer_annotation": enhancer_status[item["genome"]],
                    })
                    assignment_rows.append(assignment)
                    key = tuple(str(item[field]) for field in (
                        "entity_type", "entity_id", "sample_key", "cohort_id", "genome", "factor",
                        "condition", "replicate", "caller", "peak_class",
                    ))
                    aggregate[key][f"count:{selected_category}"] += 1
                    aggregate[key][f"bp:{selected_category}"] += end - start
                    if overlaps:
                        for feature, overlap_bp in overlaps:
                            overlap = dict(item)
                            overlap.update({
                                "peak_id": peak_id, "chrom": chrom, "start": start, "end": end,
                                "feature_category": feature.category, "feature_id": feature.feature_id,
                                "feature_start": feature.start, "feature_end": feature.end,
                                "overlap_bp": overlap_bp, "gene_id": feature.gene_id,
                                "gene_name": feature.gene_name, "strand": feature.strand,
                            })
                            overlap_rows.append(overlap)
            status_rows.append({
                **item, "input_valid_peaks": valid, "invalid_peaks": invalid,
                "unclassified_peaks": unclassified, "enhancer_annotation": enhancer_status[item["genome"]],
                "status": "SUCCESS" if valid else "EMPTY", "reason": "." if valid else "no_valid_peaks",
            })

        summary_rows: list[dict[str, object]] = []
        identity_fields = ["entity_type", "entity_id", "sample_key", "cohort_id", "genome", "factor",
                           "condition", "replicate", "caller", "peak_class"]
        for key, values in sorted(aggregate.items()):
            total = sum(values[f"count:{category}"] for category in categories)
            total_bp = sum(values[f"bp:{category}"] for category in categories)
            for category in categories:
                count, bp = int(values[f"count:{category}"]), int(values[f"bp:{category}"])
                row = dict(zip(identity_fields, key))
                row.update({
                    "category": category, "count": count,
                    "fraction": f"{count / total:.8f}" if total else "0.00000000",
                    "percentage": f"{100 * count / total:.4f}" if total else "0.0000",
                    "bp": bp, "bp_fraction": f"{bp / total_bp:.8f}" if total_bp else "0.00000000",
                    "total_peaks": int(total), "total_peak_bp": int(total_bp), "color": COLORS[category],
                })
                summary_rows.append(row)

        assignment_fields = ["entity_type", "entity_id", "sample_key", "cohort_id", "genome", "factor",
                             "condition", "replicate", "caller", "peak_class", "peak_file", "peak_id", "chrom",
                             "start", "end", "width", "primary_category", "overlapping_categories",
                             "primary_feature_id", "primary_gene_id", "primary_gene_name", "primary_overlap_bp",
                             "nearest_gene_id", "nearest_gene_name", "nearest_tss_signed_distance",
                             "enhancer_annotation"]
        overlap_fields = ["entity_type", "entity_id", "sample_key", "cohort_id", "genome", "factor",
                          "condition", "replicate", "caller", "peak_class", "peak_file", "peak_id", "chrom",
                          "start", "end", "feature_category", "feature_id", "feature_start", "feature_end",
                          "overlap_bp", "gene_id", "gene_name", "strand"]
        summary_fields = identity_fields + ["category", "count", "fraction", "percentage", "bp", "bp_fraction",
                                            "total_peaks", "total_peak_bp", "color"]
        status_fields = identity_fields + ["peak_file", "input_valid_peaks", "invalid_peaks", "unclassified_peaks",
                                           "enhancer_annotation", "status", "reason"]
        write_tsv_gz(outdir / "peak_feature_assignments.tsv.gz", assignment_fields, assignment_rows)
        write_tsv_gz(outdir / "peak_feature_all_overlaps.tsv.gz", overlap_fields, overlap_rows)
        write_tsv(outdir / "peak_feature_summary.tsv", summary_fields, summary_rows)
        count_fields = identity_fields + ["category", "count", "total_peaks"]
        fraction_fields = identity_fields + ["category", "fraction", "percentage", "total_peaks"]
        bp_fields = identity_fields + ["category", "bp", "bp_fraction", "total_peak_bp"]
        write_tsv(outdir / "peak_feature_counts.tsv", count_fields,
                  [{field: row[field] for field in count_fields} for row in summary_rows])
        write_tsv(outdir / "peak_feature_fractions.tsv", fraction_fields,
                  [{field: row[field] for field in fraction_fields} for row in summary_rows])
        write_tsv(outdir / "peak_feature_bp_coverage.tsv", bp_fields,
                  [{field: row[field] for field in bp_fields} for row in summary_rows])
        write_tsv(outdir / "peak_annotation_status.tsv", status_fields, status_rows)
        write_tsv(outdir / "peak_feature_colors.tsv", ["category", "color"],
                  [{"category": category, "color": COLORS[category]} for category in categories])
        plot_formats = [item.strip().lower() for item in args.plot_formats.split(",") if item.strip()]
        if not set(plot_formats) <= {"png", "pdf", "svg"}:
            raise ValueError("plot formats must be a subset of png,pdf,svg")
        if not args.skip_plots:
            plot_composition(outdir, summary_rows, categories, plot_formats)
        print(f"Annotated {len(assignment_rows)} peaks from {len(sets)} peak sets: {outdir}")
        return 0
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        print(f"PEAK ANNOTATION ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
