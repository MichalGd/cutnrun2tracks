#!/usr/bin/env python3
"""Join consensus gene/cCRE reference annotations to differential result tables."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def annotation_map(annotation_root: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = defaultdict(lambda: {
        "nearest_gene_name": ".", "nearest_gene_id": ".", "distance_to_gene": ".", "ccre_reference_overlaps": "."
    })
    nearest = next(annotation_root.glob("*.nearest_gene.tsv"), None)
    if nearest:
        with nearest.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 12:
                    result[fields[3]].update({"nearest_gene_name": fields[8], "nearest_gene_id": fields[9], "distance_to_gene": fields[11]})
    ccre = next(annotation_root.glob("*.ccre_reference_overlaps.tsv"), None)
    if ccre:
        labels: dict[str, set[str]] = defaultdict(set)
        with ccre.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 10 and fields[8] not in {".", "-1"}:
                    labels[fields[3]].add(fields[8])
        for region, values in labels.items():
            result[region]["ccre_reference_overlaps"] = ",".join(sorted(values))
    return result


def annotate(path: Path, annotations: dict[str, dict[str, str]]) -> None:
    destination = path.with_name(path.name.removesuffix(".tsv.gz") + ".annotated.tsv.gz")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or "region_id" not in reader.fieldnames:
            return
        extra = ["nearest_gene_name", "nearest_gene_id", "distance_to_gene", "ccre_reference_overlaps"]
        with gzip.open(destination, "wt", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames + extra, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                row.update(annotations.get(row["region_id"], {key: "." for key in extra}))
                writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("output_dir", type=Path); args = parser.parse_args()
    differential = args.output_dir / "08_differential"
    annotation = args.output_dir / "07_annotation"
    if not differential.is_dir():
        return 0
    for cohort_dir in differential.iterdir():
        if not cohort_dir.is_dir():
            continue
        mapping = annotation_map(annotation / cohort_dir.name / "consensus")
        for path in cohort_dir.rglob("*.tsv.gz"):
            if not path.name.endswith(".annotated.tsv.gz"):
                annotate(path, mapping)
    return 0


if __name__ == "__main__": raise SystemExit(main())
