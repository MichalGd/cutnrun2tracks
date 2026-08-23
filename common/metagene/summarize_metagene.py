#!/usr/bin/env python3
"""Aggregate per-task metagene metadata into reporting-friendly tables."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


FIELDS = [
    "task_id", "sample_id", "assay", "genome", "normalization", "normalization_detail",
    "gene_set", "mode", "n_genes_reference", "n_genes_matrix", "bigwig", "bed12",
    "bed12_sha256", "matrix", "profile_data", "profile_png", "profile_pdf",
    "heatmap_png", "heatmap_pdf", "sorted_regions", "status",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    try:
        rows: list[dict[str, str]] = []
        for path in sorted(args.output_dir.rglob("*.metadata.tsv")):
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                row = next(reader, None)
                if row is None:
                    raise ValueError(f"empty metadata file: {path}")
                missing = set(FIELDS) - row.keys()
                if missing:
                    raise ValueError(f"{path}: missing metadata columns: {','.join(sorted(missing))}")
                for artifact in ("matrix", "profile_data", "profile_png", "profile_pdf", "heatmap_png", "heatmap_pdf", "sorted_regions"):
                    if not Path(row[artifact]).is_file():
                        raise ValueError(f"{path}: missing declared artifact {artifact}: {row[artifact]}")
                rows.append({field: row[field] for field in FIELDS})
        if not rows:
            raise ValueError("no completed metagene metadata files found")
        task_ids = [row["task_id"] for row in rows]
        if len(task_ids) != len(set(task_ids)):
            raise ValueError("duplicate task IDs in metagene metadata")

        rows.sort(key=lambda row: (row["genome"], row["gene_set"], row["mode"], row["sample_id"]))
        with (args.output_dir / "artifacts.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)

        multiqc_dir = args.output_dir / "multiqc_data"
        multiqc_dir.mkdir(parents=True, exist_ok=True)
        with (multiqc_dir / "metagene_summary_mqc.tsv").open("w", encoding="utf-8", newline="") as handle:
            handle.write("# id: metagene_summary\n")
            handle.write("# section_name: Metagene aggregate-signal summary\n")
            handle.write("# description: Gene counts retained in each sample, gene-set, and plotting mode.\n")
            handle.write("# format: tsv\n")
            handle.write("# plot_type: table\n")
            handle.write("task\tn_genes_reference\tn_genes_matrix\tnormalization\n")
            for row in rows:
                handle.write(f"{row['task_id']}\t{row['n_genes_reference']}\t{row['n_genes_matrix']}\t{row['normalization']}\n")
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"METAGENE SUMMARY ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
