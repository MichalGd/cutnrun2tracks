#!/usr/bin/env python3
"""Record resolved immutable reference paths and checksums used by a run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shlex
import sys
from pathlib import Path


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("resolved_config", type=Path)
    parser.add_argument("sample_manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    config: dict[str, str] = {}
    for raw in args.resolved_config.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        key, value = raw.split("=", 1)
        parsed = shlex.split(value)
        config[key] = parsed[0] if parsed else ""
    with args.sample_manifest.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    genomes = sorted({row["genome"] for row in samples})
    paths: list[tuple[str, str, Path]] = []
    for blacklist in sorted({row["blacklist"] for row in samples}):
        paths.append(("samplesheet", "blacklist", Path(blacklist)))
    for genome in genomes:
        suffix = genome.upper().replace(".", "_").replace("-", "_")
        for kind in ("FASTA", "CHROM_SIZES", "CANONICAL_CONTIGS", "GTF", "TSS_BED", "CCRE_BED"):
            value = config.get(f"{kind}_{suffix}", "")
            if value:
                paths.append((genome, kind.lower(), Path(value)))
        index = config.get(f"INDEX_{suffix}", "")
        if index:
            for item in sorted(Path(index).parent.glob(Path(index).name + ".*.bt2*")):
                paths.append((genome, "bowtie2_index", item))
    if config.get("SPIKEIN_MODE") != "none":
        for kind in ("SPIKEIN_FASTA", "SPIKEIN_CHROM_SIZES", "SPIKEIN_ALLOWED_CONTIGS", "SPIKEIN_BLACKLIST"):
            value = config.get(kind, "")
            if value:
                paths.append((config.get("SPIKEIN_REFERENCE_ID", "."), kind.lower(), Path(value)))
        index = config.get("SPIKEIN_INDEX", "")
        if index:
            for item in sorted(Path(index).parent.glob(Path(index).name + ".*.bt2*")):
                paths.append((config.get("SPIKEIN_REFERENCE_ID", "."), "composite_bowtie2_index", item))
    if config.get("RUN_METAGENE") == "true":
        gene_manifest = Path(config.get("METAGENE_GENE_SET_MANIFEST", ""))
        if not gene_manifest.is_file():
            print(f"REFERENCE ERROR: missing metagene gene-set manifest: {gene_manifest}", file=sys.stderr)
            return 1
        paths.append(("metagene", "gene_set_manifest", gene_manifest))
        try:
            with gene_manifest.open(encoding="utf-8", newline="") as handle:
                gene_rows = list(csv.DictReader(handle, delimiter="\t"))
            for row in gene_rows:
                bed12 = Path(row.get("bed12", ""))
                paths.append((row.get("genome", "metagene"), f"metagene_{row.get('gene_set_id', 'gene_set')}", bed12))
        except (OSError, csv.Error) as exc:
            print(f"REFERENCE ERROR: invalid metagene gene-set manifest: {exc}", file=sys.stderr)
            return 1
    records = []
    for reference_id, kind, path in paths:
        if not path.is_file():
            print(f"REFERENCE ERROR: missing {kind}: {path}", file=sys.stderr)
            return 1
        records.append({"reference_id": reference_id, "kind": kind, "resolved_path": str(path.resolve()),
                        "size": path.stat().st_size, "sha256": checksum(path)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reference_id", "kind", "resolved_path", "size", "sha256"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
