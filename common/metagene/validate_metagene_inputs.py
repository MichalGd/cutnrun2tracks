#!/usr/bin/env python3
"""Validate metagene manifests and resolve sample-by-gene-set-by-mode tasks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path


SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
TRACK_FIELDS = {"sample_id", "assay", "genome", "bigwig", "normalization", "normalization_detail", "blacklist", "chrom_sizes"}
GENE_FIELDS = {"gene_set_id", "genome", "bed12", "label"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = required - fields
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def read_chrom_sizes(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{number}: invalid chrom sizes row")
            result[fields[0]] = int(fields[1])
    if not result:
        raise ValueError(f"{path}: empty chrom sizes")
    return result


def validate_bed12(path: Path, chrom_sizes: dict[str, int] | None = None) -> tuple[int, set[str]]:
    count = 0
    names: set[str] = set()
    contigs: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{number}: BED12 required")
            chrom, name, strand = fields[0], fields[3], fields[5]
            start, end, blocks = int(fields[1]), int(fields[2]), int(fields[9])
            sizes = [int(value) for value in fields[10].rstrip(",").split(",") if value]
            starts = [int(value) for value in fields[11].rstrip(",").split(",") if value]
            if start < 0 or end <= start or strand not in {"+", "-"}:
                raise ValueError(f"{path}:{number}: invalid BED coordinates or strand")
            if blocks <= 0 or len(sizes) != blocks or len(starts) != blocks:
                raise ValueError(f"{path}:{number}: blockCount does not match BED12 blocks")
            if any(size <= 0 for size in sizes) or any(offset < 0 for offset in starts):
                raise ValueError(f"{path}:{number}: invalid BED12 block")
            if any(start + offset + size > end for offset, size in zip(starts, sizes)):
                raise ValueError(f"{path}:{number}: exon block outside transcript span")
            if name in names:
                raise ValueError(f"{path}:{number}: duplicate region name {name}")
            if chrom_sizes and (chrom not in chrom_sizes or end > chrom_sizes[chrom]):
                raise ValueError(f"{path}:{number}: interval is incompatible with chrom sizes")
            names.add(name); contigs.add(chrom); count += 1
    if count == 0:
        raise ValueError(f"{path}: no BED12 records")
    return count, contigs


def validate_bigwig(path: Path, chrom_sizes: dict[str, int], gene_contigs: set[str]) -> list[str]:
    try:
        import pyBigWig  # type: ignore
    except ImportError as exc:
        raise ValueError("pyBigWig is required to validate bigWig inputs") from exc
    handle = pyBigWig.open(str(path))
    if handle is None or not handle.isBigWig():
        raise ValueError(f"not a readable bigWig: {path}")
    try:
        bigwig_chroms = handle.chroms()
    finally:
        handle.close()
    if not bigwig_chroms:
        raise ValueError(f"bigWig has no chromosomes: {path}")
    mismatched = [chrom for chrom, length in bigwig_chroms.items() if chrom in chrom_sizes and chrom_sizes[chrom] != length]
    unknown = sorted(set(bigwig_chroms) - chrom_sizes.keys())
    if mismatched:
        raise ValueError(f"{path}: chromosome length mismatch for {','.join(sorted(mismatched)[:10])}")
    if unknown:
        raise ValueError(f"{path}: chromosomes absent from chrom sizes: {','.join(unknown[:10])}")
    if not set(bigwig_chroms) & gene_contigs:
        raise ValueError(f"{path}: no chromosomes overlap the BED12 gene model")
    missing = sorted(gene_contigs - bigwig_chroms.keys())
    return [f"bigWig omits {len(missing)} gene-model contigs; missing values will follow configured policy"] if missing else []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--track-manifest", type=Path, required=True)
    parser.add_argument("--gene-set-manifest", type=Path, required=True)
    parser.add_argument("--gene-sets", required=True)
    parser.add_argument("--modes", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--warnings", type=Path)
    parser.add_argument("--skip-bigwig-validation", action="store_true")
    args = parser.parse_args()

    try:
        tracks = read_tsv(args.track_manifest, TRACK_FIELDS)
        gene_sets = read_tsv(args.gene_set_manifest, GENE_FIELDS)
        selected_sets = [item.strip() for item in args.gene_sets.split(",") if item.strip()]
        modes = [item.strip().lower() for item in args.modes.split(",") if item.strip()]
        if not selected_sets or any(not SAFE_ID.fullmatch(item) for item in selected_sets):
            raise ValueError("--gene-sets must be a nonempty comma-separated list of safe IDs")
        if not modes or not set(modes) <= {"tss", "tes", "gene_body"}:
            raise ValueError("--modes must be a subset of tss,tes,gene_body")
        selected_sets = list(dict.fromkeys(selected_sets)); modes = list(dict.fromkeys(modes))

        gene_index: dict[tuple[str, str], dict[str, str]] = {}
        bed_cache: dict[tuple[Path, Path], tuple[int, set[str]]] = {}
        for row in gene_sets:
            gene_set_id, genome = row["gene_set_id"], row["genome"]
            if not SAFE_ID.fullmatch(gene_set_id) or not SAFE_ID.fullmatch(genome):
                raise ValueError(f"unsafe gene-set or genome ID: {gene_set_id}/{genome}")
            key = (gene_set_id, genome)
            if key in gene_index:
                raise ValueError(f"duplicate gene-set/genome key: {gene_set_id}/{genome}")
            bed = Path(row["bed12"])
            if not bed.is_file():
                raise ValueError(f"gene-set BED12 missing: {bed}")
            expected = row.get("sha256", "")
            if expected and expected != "." and expected.lower() != sha256(bed):
                raise ValueError(f"gene-set checksum mismatch: {bed}")
            gene_index[key] = row

        task_rows: list[dict[str, str]] = []
        warning_rows: list[dict[str, str]] = []
        normalization_groups: dict[tuple[str, str], set[str]] = {}
        for track in tracks:
            key = (track["genome"], track["assay"])
            normalization_groups.setdefault(key, set()).add(track["normalization"])
        for (genome, assay), normalizations in sorted(normalization_groups.items()):
            if len(normalizations) > 1:
                warning_rows.append({
                    "sample_id": ".", "gene_set_id": ".",
                    "warning": f"{genome}/{assay} contains mixed normalization families: {','.join(sorted(normalizations))}",
                })
        seen_samples: set[str] = set()
        for track in tracks:
            sample = track["sample_id"]
            if not SAFE_ID.fullmatch(sample):
                raise ValueError(f"unsafe sample ID: {sample}")
            if sample in seen_samples:
                raise ValueError(f"duplicate sample ID in track manifest: {sample}")
            seen_samples.add(sample)
            if not track["normalization"]:
                raise ValueError(f"{sample}: normalization label is required")
            for field in ("bigwig", "blacklist", "chrom_sizes"):
                if not Path(track[field]).is_file():
                    raise ValueError(f"{sample}: missing {field}: {track[field]}")
            chrom_sizes_path = Path(track["chrom_sizes"])
            sizes = read_chrom_sizes(chrom_sizes_path)
            for gene_set_id in selected_sets:
                key = (gene_set_id, track["genome"])
                if key not in gene_index:
                    raise ValueError(f"no {gene_set_id} gene set for genome {track['genome']}")
                gene = gene_index[key]
                if any(character in gene["label"] for character in ("\t", "\r", "\n")):
                    raise ValueError(f"{gene_set_id}/{track['genome']}: label contains a tab or newline")
                bed_path = Path(gene["bed12"])
                cache_key = (bed_path.resolve(), chrom_sizes_path.resolve())
                if cache_key not in bed_cache:
                    bed_cache[cache_key] = validate_bed12(bed_path, sizes)
                n_genes, contigs = bed_cache[cache_key]
                if not args.skip_bigwig_validation:
                    for warning in validate_bigwig(Path(track["bigwig"]), sizes, contigs):
                        warning_rows.append({"sample_id": sample, "gene_set_id": gene_set_id, "warning": warning})
                for mode in modes:
                    task_id = ".".join([sample, track["assay"], track["normalization"], track["genome"], gene_set_id, mode])
                    if not all(SAFE_ID.fullmatch(item) for item in (track["assay"], track["normalization"], track["genome"])):
                        raise ValueError(f"{sample}: unsafe assay, normalization, or genome label")
                    task_rows.append({
                        "task_id": task_id, "sample_id": sample, "assay": track["assay"],
                        "genome": track["genome"], "bigwig": str(Path(track["bigwig"]).resolve()),
                        "normalization": track["normalization"],
                        "normalization_detail": track["normalization_detail"] or ".",
                        "blacklist": str(Path(track["blacklist"]).resolve()),
                        "chrom_sizes": str(chrom_sizes_path.resolve()), "gene_set_id": gene_set_id,
                        "gene_set_label": gene["label"], "bed12": str(bed_path.resolve()),
                        "bed12_sha256": sha256(bed_path), "n_genes_reference": str(n_genes), "mode": mode,
                    })

        args.output.parent.mkdir(parents=True, exist_ok=True)
        fields = ["task_id", "sample_id", "assay", "genome", "bigwig", "normalization", "normalization_detail", "blacklist", "chrom_sizes", "gene_set_id", "gene_set_label", "bed12", "bed12_sha256", "n_genes_reference", "mode"]
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(task_rows)
        if args.warnings:
            args.warnings.parent.mkdir(parents=True, exist_ok=True)
            with args.warnings.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sample_id", "gene_set_id", "warning"], delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(warning_rows)
        print(f"Validated {len(tracks)} tracks and resolved {len(task_rows)} metagene tasks")
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"METAGENE INPUT ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
