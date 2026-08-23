#!/usr/bin/env python3
"""Build human HPA broadly expressed and optional mouse ortholog BED12 sets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import sys
import zipfile
from pathlib import Path
from typing import TextIO


SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def stable(identifier: str) -> str:
    return re.sub(r"\.\d+(?=_PAR_Y$|$)", "", identifier.strip())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class InputTable:
    def __init__(self, path: Path):
        self.path = path
        self.archive: zipfile.ZipFile | None = None
        self.binary: io.BufferedReader | gzip.GzipFile | None = None
        self.text: TextIO | None = None

    def __enter__(self) -> TextIO:
        if self.path.suffix == ".zip":
            self.archive = zipfile.ZipFile(self.path)
            members = [name for name in self.archive.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise ValueError(f"ZIP must contain exactly one table: {self.path}")
            self.binary = self.archive.open(members[0])  # type: ignore[assignment]
            self.text = io.TextIOWrapper(self.binary, encoding="utf-8")
        elif self.path.suffix == ".gz":
            self.binary = gzip.open(self.path, "rb")
            self.text = io.TextIOWrapper(self.binary, encoding="utf-8")
        else:
            self.text = self.path.open(encoding="utf-8", newline="")
        return self.text

    def __exit__(self, *_: object) -> None:
        if self.text:
            self.text.close()
        if self.archive:
            self.archive.close()


def read_bed12(path: Path) -> tuple[list[str], dict[str, str]]:
    order: list[str] = []
    rows: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{number}: BED12 required")
            gene_id = stable(fields[3])
            if gene_id in rows:
                raise ValueError(f"{path}:{number}: duplicate gene ID {gene_id}")
            order.append(gene_id)
            rows[gene_id] = raw if raw.endswith("\n") else raw + "\n"
    if not rows:
        raise ValueError(f"BED12 contains no records: {path}")
    return order, rows


def find_column(fieldnames: list[str], requested: str) -> str:
    normalized = {name.strip().lower(): name for name in fieldnames}
    key = requested.strip().lower()
    if key not in normalized:
        raise ValueError(f"missing column {requested!r}; available columns: {', '.join(fieldnames)}")
    return normalized[key]


def hpa_ids(args: argparse.Namespace) -> set[str]:
    with InputTable(args.hpa_table) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("HPA table has no header")
        gene_column = find_column(reader.fieldnames, args.hpa_gene_column)
        specificity_column = find_column(reader.fieldnames, args.hpa_specificity_column)
        distribution_column = find_column(reader.fieldnames, args.hpa_distribution_column)
        result = {
            stable(row[gene_column])
            for row in reader
            if (row.get(specificity_column) or "").strip().casefold() == args.specificity.strip().casefold()
            and (row.get(distribution_column) or "").strip().casefold() == args.distribution.strip().casefold()
            and stable(row.get(gene_column) or "")
        }
    if not result:
        raise ValueError("HPA filters selected no genes")
    return result


def mouse_orthologs(args: argparse.Namespace, human_ids: set[str]) -> tuple[set[str], list[dict[str, str]]]:
    if not args.ortholog_table:
        return set(), []
    with InputTable(args.ortholog_table) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("ortholog table has no header")
        human_column = find_column(reader.fieldnames, args.ortholog_human_column)
        mouse_column = find_column(reader.fieldnames, args.ortholog_mouse_column)
        type_column = find_column(reader.fieldnames, args.ortholog_type_column)
        confidence_column = find_column(reader.fieldnames, args.ortholog_confidence_column)
        selected: list[dict[str, str]] = []
        mouse: set[str] = set()
        for row in reader:
            human = stable(row.get(human_column) or "")
            target = stable(row.get(mouse_column) or "")
            orthology_type = (row.get(type_column) or "").strip().lower()
            confidence = (row.get(confidence_column) or "").strip()
            if human not in human_ids or not target:
                continue
            if args.orthology_policy == "one2one" and orthology_type not in {"ortholog_one2one", "one2one", "1-to-1"}:
                continue
            if args.require_confidence and confidence not in {"1", "1.0", "high"}:
                continue
            mouse.add(target)
            selected.append({
                "human_gene_id": human, "mouse_gene_id": target,
                "orthology_type": orthology_type, "orthology_confidence": confidence or ".",
            })
    if not mouse:
        raise ValueError("orthology filters selected no mouse genes")
    return mouse, selected


def write_subset(path: Path, order: list[str], rows: dict[str, str], selected: set[str]) -> tuple[int, set[str]]:
    present = selected & rows.keys()
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for gene_id in order:
            if gene_id in present:
                handle.write(rows[gene_id])
    return len(present), selected - rows.keys()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hpa-table", type=Path, required=True)
    parser.add_argument("--human-gene-model", type=Path, required=True)
    parser.add_argument("--human-assembly", default="GRCh38")
    parser.add_argument("--hpa-release", required=True)
    parser.add_argument("--hpa-gene-column", default="Gene")
    parser.add_argument("--hpa-specificity-column", default="RNA tissue specificity")
    parser.add_argument("--hpa-distribution-column", default="RNA tissue distribution")
    parser.add_argument("--specificity", default="Low tissue specificity")
    parser.add_argument("--distribution", default="Detected in all")
    parser.add_argument("--gene-set-id", default="broadly_expressed")
    parser.add_argument("--label", default="HPA broadly expressed")
    parser.add_argument("--ortholog-table", type=Path)
    parser.add_argument("--mouse-gene-model", type=Path)
    parser.add_argument("--mouse-assembly", default="GRCm39")
    parser.add_argument("--ensembl-release", default="not_applicable")
    parser.add_argument("--ortholog-human-column", default="Gene stable ID")
    parser.add_argument("--ortholog-mouse-column", default="Mouse gene stable ID")
    parser.add_argument("--ortholog-type-column", default="Mouse homology type")
    parser.add_argument("--ortholog-confidence-column", default="Mouse orthology confidence [0 low, 1 high]")
    parser.add_argument("--orthology-policy", choices=["one2one", "all"], default="one2one")
    parser.add_argument("--require-confidence", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    try:
        if not SAFE_ID.fullmatch(args.gene_set_id):
            raise ValueError("gene-set ID must contain only letters, numbers, dot, underscore, and hyphen")
        for path in (args.hpa_table, args.human_gene_model):
            if not path.is_file():
                raise ValueError(f"input does not exist: {path}")
        if bool(args.ortholog_table) != bool(args.mouse_gene_model):
            raise ValueError("--ortholog-table and --mouse-gene-model must be supplied together")
        human_selected = hpa_ids(args)
        human_order, human_rows = read_bed12(args.human_gene_model)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        human_output = args.output_dir / f"{args.gene_set_id}.{args.human_assembly}.bed12"
        human_count, missing_human = write_subset(human_output, human_order, human_rows, human_selected)
        if human_count == 0:
            raise ValueError("none of the selected HPA IDs occur in the human gene model")

        manifest_rows = [{
            "gene_set_id": "protein_coding", "genome": args.human_assembly,
            "bed12": str(args.human_gene_model.resolve()), "label": "All filtered protein-coding genes",
            "source_release": "caller_provided_annotation", "sha256": sha256(args.human_gene_model),
            "n_genes": str(len(human_rows)),
        }, {
            "gene_set_id": args.gene_set_id, "genome": args.human_assembly,
            "bed12": str(human_output.resolve()), "label": args.label,
            "source_release": f"HPA_{args.hpa_release}", "sha256": sha256(human_output),
            "n_genes": str(human_count),
        }]
        with (args.output_dir / "hpa_unmatched_human_ids.tsv").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("gene_id\n")
            for gene_id in sorted(missing_human):
                handle.write(gene_id + "\n")

        if args.ortholog_table and args.mouse_gene_model:
            mouse_selected, mapping_rows = mouse_orthologs(args, human_selected)
            mouse_order, mouse_rows = read_bed12(args.mouse_gene_model)
            mouse_output = args.output_dir / f"{args.gene_set_id}.{args.mouse_assembly}.bed12"
            mouse_count, missing_mouse = write_subset(mouse_output, mouse_order, mouse_rows, mouse_selected)
            if mouse_count == 0:
                raise ValueError("none of the selected mouse orthologs occur in the mouse gene model")
            manifest_rows.append({
                "gene_set_id": "protein_coding", "genome": args.mouse_assembly,
                "bed12": str(args.mouse_gene_model.resolve()), "label": "All filtered protein-coding genes",
                "source_release": "caller_provided_annotation", "sha256": sha256(args.mouse_gene_model),
                "n_genes": str(len(mouse_rows)),
            })
            manifest_rows.append({
                "gene_set_id": args.gene_set_id, "genome": args.mouse_assembly,
                "bed12": str(mouse_output.resolve()),
                "label": f"Mouse one-to-one orthologs of {args.label}",
                "source_release": f"HPA_{args.hpa_release};Ensembl_{args.ensembl_release}",
                "sha256": sha256(mouse_output), "n_genes": str(mouse_count),
            })
            with (args.output_dir / "selected_human_mouse_orthologs.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["human_gene_id", "mouse_gene_id", "orthology_type", "orthology_confidence"], delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(mapping_rows)
            with (args.output_dir / "hpa_unmatched_mouse_ids.tsv").open("w", encoding="utf-8", newline="\n") as handle:
                handle.write("gene_id\n")
                for gene_id in sorted(missing_mouse):
                    handle.write(gene_id + "\n")

        with (args.output_dir / "gene_sets.tsv").open("w", encoding="utf-8", newline="") as handle:
            fields = ["gene_set_id", "genome", "bed12", "label", "source_release", "sha256", "n_genes"]
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(manifest_rows)
        with (args.output_dir / "hpa_source_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["source", "release", "path", "sha256"])
            writer.writerow(["HPA", args.hpa_release, str(args.hpa_table.resolve()), sha256(args.hpa_table)])
            if args.ortholog_table:
                writer.writerow(["Ensembl_BioMart", args.ensembl_release, str(args.ortholog_table.resolve()), sha256(args.ortholog_table)])
        print(f"Prepared {human_count} human HPA reference genes in {human_output}")
        return 0
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        print(f"HPA SUBSET ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
