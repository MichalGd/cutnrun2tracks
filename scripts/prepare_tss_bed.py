#!/usr/bin/env python3
"""Create a strand-aware, zero-based one-base TSS BED from GTF genes."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("gtf", type=Path); parser.add_argument("output", type=Path); args = parser.parse_args()
    records = set()
    with args.gtf.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            start, end, strand = int(fields[3]), int(fields[4]), fields[6]
            match = re.search(r'gene_id "([^"]+)"', fields[8]); gene = match.group(1) if match else "."
            position = start - 1 if strand == "+" else end - 1
            records.add((fields[0], position, position + 1, gene, strand))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for chrom, start, end, gene, strand in sorted(records, key=lambda item: (item[0], item[1], item[3])):
            handle.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t{strand}\n")
    if not records:
        raise ValueError("no gene records found in GTF")
    return 0


if __name__ == "__main__": raise SystemExit(main())
