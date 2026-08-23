#!/usr/bin/env python3
"""Build deterministic, filtered BED12 metagene models from a GTF annotation."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO


ATTRIBUTE_RE = re.compile(r"(?:^|;\s*)([^\s;]+)\s+(?:\"([^\"]*)\"|([^;\s]+))")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def attributes(text: str) -> dict[str, list[str]]:
    result: dict[str, list[str]] = defaultdict(list)
    for match in ATTRIBUTE_RE.finditer(text):
        result[match.group(1)].append(match.group(2) if match.group(2) is not None else match.group(3))
    return dict(result)


def first(values: dict[str, list[str]], *keys: str) -> str:
    for key in keys:
        if values.get(key):
            return values[key][0]
    return ""


def stable(identifier: str) -> str:
    return re.sub(r"\.\d+(?=_PAR_Y$|$)", "", identifier)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))


def canonical_priority(tags: set[str]) -> int:
    normalized = {tag.lower().replace("-", "_") for tag in tags}
    if "mane_select" in normalized:
        return 5
    if "ensembl_canonical" in normalized:
        return 4
    if "appris_principal_1" in normalized:
        return 3
    if any(tag.startswith("appris_principal") for tag in normalized):
        return 2
    if "basic" in normalized:
        return 1
    return 0


def choose_transcript(candidates: list[dict[str, object]], policy: str) -> dict[str, object]:
    def measurements(item: dict[str, object]) -> tuple[int, int, int, int]:
        exons = item["exons"]
        cds = item["cds"]
        assert isinstance(exons, list) and isinstance(cds, list)
        spliced = interval_length(exons)
        cds_length = interval_length(cds)
        span = max(end for _, end in exons) - min(start for start, _ in exons)
        tags = item["tags"]
        assert isinstance(tags, set)
        return canonical_priority(tags), cds_length, spliced, span

    def key(item: dict[str, object]) -> tuple[object, ...]:
        priority, cds_length, spliced, span = measurements(item)
        transcript_id = str(item["transcript_id"])
        if policy == "canonical_then_longest":
            return (-priority, -cds_length, -spliced, -span, transcript_id)
        if policy == "longest_cds":
            return (-cds_length, -spliced, -span, -priority, transcript_id)
        return (-spliced, -cds_length, -span, -priority, transcript_id)

    return sorted(candidates, key=key)[0]


def read_contigs(path: Path) -> tuple[list[str], set[str]]:
    ordered: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#"):
                ordered.append(line.split()[0])
    if not ordered:
        raise ValueError(f"canonical contig file is empty: {path}")
    return ordered, set(ordered)


def read_chrom_sizes(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{number}: expected chromosome and length")
            result[fields[0]] = int(fields[1])
    return result


def read_blacklist(path: Path | None) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if path is None:
        return result
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{number}: expected BED3 or greater")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: invalid interval")
            result[fields[0]].append((start, end))
    for chrom in result:
        result[chrom].sort()
    return result


def overlaps(intervals: list[tuple[int, int]], blacklist: list[tuple[int, int]]) -> bool:
    if not blacklist:
        return False
    for start, end in intervals:
        for blocked_start, blocked_end in blacklist:
            if blocked_start >= end:
                break
            if blocked_end > start:
                return True
    return False


def parse_gtf(path: Path) -> dict[str, list[dict[str, object]]]:
    transcripts: dict[str, dict[str, object]] = {}
    protein_coding_genes: set[str] = set()
    with open_text(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected nine GTF columns")
            chrom, _, feature, raw_start, raw_end, _, strand, _, raw_attributes = fields
            if strand not in {"+", "-"}:
                continue
            attrs = attributes(raw_attributes)
            gene_id = first(attrs, "gene_id")
            gene_type = first(attrs, "gene_type", "gene_biotype")
            if not gene_id:
                continue
            if gene_type == "protein_coding":
                protein_coding_genes.add(gene_id)
            transcript_id = first(attrs, "transcript_id")
            if not transcript_id:
                continue
            transcript_type = first(attrs, "transcript_type", "transcript_biotype")
            record = transcripts.setdefault(transcript_id, {
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "gene_name": first(attrs, "gene_name") or stable(gene_id),
                "gene_type": gene_type,
                "transcript_type": transcript_type,
                "chrom": chrom,
                "strand": strand,
                "tags": set(),
                "exons": [],
                "cds": [],
            })
            if record["gene_id"] != gene_id or record["chrom"] != chrom or record["strand"] != strand:
                raise ValueError(f"{path}:{number}: inconsistent transcript {transcript_id}")
            tags = record["tags"]
            assert isinstance(tags, set)
            tags.update(attrs.get("tag", []))
            if feature in {"exon", "CDS"}:
                start, end = int(raw_start) - 1, int(raw_end)
                if start < 0 or end <= start:
                    raise ValueError(f"{path}:{number}: invalid GTF coordinates")
                target = record["exons"] if feature == "exon" else record["cds"]
                assert isinstance(target, list)
                target.append((start, end))

    by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for record in transcripts.values():
        exons = record["exons"]
        if not isinstance(exons, list) or not exons:
            continue
        gene_id = str(record["gene_id"])
        is_pc_gene = gene_id in protein_coding_genes or record["gene_type"] == "protein_coding"
        transcript_type = str(record["transcript_type"])
        if is_pc_gene and transcript_type in {"", "protein_coding"}:
            record["exons"] = merge_intervals(exons)
            cds = record["cds"]
            assert isinstance(cds, list)
            record["cds"] = merge_intervals(cds)
            by_gene[gene_id].append(record)
    if not by_gene:
        raise ValueError("no protein-coding genes with exon-bearing transcripts were found")
    return by_gene


def bed12(record: dict[str, object]) -> str:
    exons = record["exons"]
    assert isinstance(exons, list) and exons
    chrom_start = min(start for start, _ in exons)
    chrom_end = max(end for _, end in exons)
    block_sizes = ",".join(str(end - start) for start, end in exons) + ","
    block_starts = ",".join(str(start - chrom_start) for start, _ in exons) + ","
    return "\t".join([
        str(record["chrom"]), str(chrom_start), str(chrom_end), stable(str(record["gene_id"])),
        "0", str(record["strand"]), str(chrom_start), str(chrom_end), "0",
        str(len(exons)), block_sizes, block_starts,
    ])


def overlap_exclusions(models: list[dict[str, object]], adjacency: int, strand_policy: str) -> set[str]:
    excluded: set[str] = set()
    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    for model in models:
        by_chrom[str(model["chrom"])].append(model)
    for chrom_models in by_chrom.values():
        chrom_models.sort(key=lambda item: (int(item["span_start"]), int(item["span_end"]), str(item["gene_id"])))
        active: list[dict[str, object]] = []
        for current in chrom_models:
            current_start = int(current["span_start"])
            active = [item for item in active if int(item["span_end"]) + adjacency >= current_start]
            for other in active:
                if strand_policy == "same" and other["strand"] != current["strand"]:
                    continue
                excluded.add(str(other["gene_id"]))
                excluded.add(str(current["gene_id"]))
            active.append(current)
    return excluded


def write_bed(path: Path, models: list[dict[str, object]], contig_order: dict[str, int]) -> None:
    ordered = sorted(models, key=lambda item: (
        contig_order.get(str(item["chrom"]), len(contig_order)),
        int(item["span_start"]), int(item["span_end"]), stable(str(item["gene_id"])),
    ))
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for model in ordered:
            handle.write(bed12(model) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--annotation-source", choices=["gencode", "ensembl"], required=True)
    parser.add_argument("--annotation-release", required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--canonical-contigs", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path)
    parser.add_argument("--transcript-policy", choices=["canonical_then_longest", "longest_cds", "longest_spliced"], default="canonical_then_longest")
    parser.add_argument("--min-gene-span", type=int, default=1000)
    parser.add_argument("--min-spliced-length", type=int, default=500)
    parser.add_argument("--blacklist-policy", choices=["none", "gene_span", "exon"], default="gene_span")
    parser.add_argument("--overlap-policy", choices=["keep", "exclude"], default="keep")
    parser.add_argument("--overlap-strand-policy", choices=["any", "same"], default="any")
    parser.add_argument("--adjacency-bp", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    try:
        for path in (args.gtf, args.chrom_sizes, args.canonical_contigs):
            if not path.is_file():
                raise ValueError(f"input does not exist: {path}")
        if args.blacklist and not args.blacklist.is_file():
            raise ValueError(f"blacklist does not exist: {args.blacklist}")
        if min(args.min_gene_span, args.min_spliced_length, args.adjacency_bp) < 0:
            raise ValueError("length and adjacency thresholds must be non-negative")
        contigs, canonical = read_contigs(args.canonical_contigs)
        chrom_sizes = read_chrom_sizes(args.chrom_sizes)
        missing_sizes = canonical - chrom_sizes.keys()
        if missing_sizes:
            raise ValueError("canonical contigs absent from chrom sizes: " + ",".join(sorted(missing_sizes)))
        blocked = read_blacklist(args.blacklist)
        by_gene = parse_gtf(args.gtf)
    except (OSError, ValueError) as exc:
        print(f"REFERENCE ERROR: {exc}", file=sys.stderr)
        return 1

    models: list[dict[str, object]] = []
    seen_stable: dict[str, str] = {}
    for gene_id, candidates in sorted(by_gene.items()):
        selected = choose_transcript(candidates, args.transcript_policy)
        exons = selected["exons"]
        cds = selected["cds"]
        assert isinstance(exons, list) and isinstance(cds, list)
        model = dict(selected)
        model["span_start"] = min(start for start, _ in exons)
        model["span_end"] = max(end for _, end in exons)
        model["span_length"] = int(model["span_end"]) - int(model["span_start"])
        model["spliced_length"] = interval_length(exons)
        model["cds_length"] = interval_length(cds)
        model["selection_priority"] = canonical_priority(selected["tags"])  # type: ignore[arg-type]
        model["reasons"] = []
        stable_id = stable(gene_id)
        if stable_id in seen_stable and seen_stable[stable_id] != gene_id:
            print(f"REFERENCE ERROR: stable ID collision: {stable_id}", file=sys.stderr)
            return 1
        seen_stable[stable_id] = gene_id
        reasons = model["reasons"]
        assert isinstance(reasons, list)
        chrom = str(model["chrom"])
        if chrom not in canonical:
            reasons.append("noncanonical_contig")
        elif int(model["span_end"]) > chrom_sizes[chrom]:
            reasons.append("outside_chromosome_bounds")
        if int(model["span_length"]) < args.min_gene_span:
            reasons.append("below_min_gene_span")
        if int(model["spliced_length"]) < args.min_spliced_length:
            reasons.append("below_min_spliced_length")
        if args.blacklist_policy != "none" and chrom in blocked:
            query = [(int(model["span_start"]), int(model["span_end"]))] if args.blacklist_policy == "gene_span" else exons
            if overlaps(query, blocked[chrom]):
                reasons.append("blacklist_overlap")
        models.append(model)

    initially_eligible = [model for model in models if not model["reasons"]]
    if args.overlap_policy == "exclude":
        colliding = overlap_exclusions(initially_eligible, args.adjacency_bp, args.overlap_strand_policy)
        reason = "overlapping_or_adjacent_gene" if args.adjacency_bp else "overlapping_gene"
        for model in initially_eligible:
            if str(model["gene_id"]) in colliding:
                model["reasons"].append(reason)  # type: ignore[union-attr]
    eligible = [model for model in models if not model["reasons"]]
    if not eligible:
        print("REFERENCE ERROR: all protein-coding genes were excluded", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    contig_order = {contig: index for index, contig in enumerate(contigs)}
    all_bed = args.output_dir / "protein_coding.all.bed12"
    filtered_bed = args.output_dir / "protein_coding.filtered.bed12"
    write_bed(all_bed, models, contig_order)
    write_bed(filtered_bed, eligible, contig_order)

    metadata_fields = [
        "gene_id", "stable_gene_id", "gene_name", "transcript_id", "stable_transcript_id",
        "chrom", "strand", "span_start", "span_end", "span_length", "spliced_length",
        "cds_length", "selection_priority", "eligible", "exclusion_reasons",
    ]
    with (args.output_dir / "gene_model_metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for model in models:
            writer.writerow({
                "gene_id": model["gene_id"], "stable_gene_id": stable(str(model["gene_id"])),
                "gene_name": model["gene_name"], "transcript_id": model["transcript_id"],
                "stable_transcript_id": stable(str(model["transcript_id"])), "chrom": model["chrom"],
                "strand": model["strand"], "span_start": model["span_start"], "span_end": model["span_end"],
                "span_length": model["span_length"], "spliced_length": model["spliced_length"],
                "cds_length": model["cds_length"], "selection_priority": model["selection_priority"],
                "eligible": "TRUE" if not model["reasons"] else "FALSE",
                "exclusion_reasons": ",".join(model["reasons"]) if model["reasons"] else ".",
            })
    with (args.output_dir / "excluded_genes.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_name", "exclusion_reasons"])
        for model in models:
            if model["reasons"]:
                writer.writerow([stable(str(model["gene_id"])), model["gene_name"], ",".join(model["reasons"])])

    manifest_fields = ["assembly", "annotation_source", "annotation_release", "transcript_policy", "blacklist_policy", "overlap_policy", "adjacency_bp", "n_all", "n_filtered", "gtf", "gtf_sha256", "filtered_bed12", "filtered_bed12_sha256"]
    with (args.output_dir / "reference_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "assembly": args.assembly, "annotation_source": args.annotation_source,
            "annotation_release": args.annotation_release, "transcript_policy": args.transcript_policy,
            "blacklist_policy": args.blacklist_policy, "overlap_policy": args.overlap_policy,
            "adjacency_bp": args.adjacency_bp, "n_all": len(models), "n_filtered": len(eligible),
            "gtf": str(args.gtf.resolve()), "gtf_sha256": sha256(args.gtf),
            "filtered_bed12": str(filtered_bed.resolve()), "filtered_bed12_sha256": sha256(filtered_bed),
        })
    with (args.output_dir / "gene_sets.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["gene_set_id", "genome", "bed12", "label", "source_release", "sha256", "n_genes"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "gene_set_id": "protein_coding", "genome": args.assembly,
            "bed12": str(filtered_bed.resolve()), "label": "All filtered protein-coding genes",
            "source_release": f"{args.annotation_source}_{args.annotation_release}",
            "sha256": sha256(filtered_bed), "n_genes": len(eligible),
        })
    print(f"Prepared {len(eligible)} of {len(models)} protein-coding gene models: {filtered_bed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
