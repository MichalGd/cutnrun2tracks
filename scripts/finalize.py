#!/usr/bin/env python3
"""Write checksums of retained result files after successful cleanup."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("output_dir", type=Path); args = parser.parse_args()
    destination = args.output_dir / "00_metadata/final_checksums.sha256"
    files = [path for path in args.output_dir.rglob("*") if path.is_file() and ".checkpoints" not in path.parts and path != destination]
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for path in sorted(files):
            handle.write(f"{digest(path)}  {path.relative_to(args.output_dir).as_posix()}\n")
    return 0


if __name__ == "__main__": raise SystemExit(main())
