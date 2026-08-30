#!/usr/bin/env python3
"""Write checksums of retained result files after successful cleanup."""

from __future__ import annotations

import argparse
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("output_dir", type=Path)
    parser.add_argument("--jobs", type=int, default=1); args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be a positive integer")
    destination = args.output_dir / "00_metadata/final_checksums.sha256"
    timing = args.output_dir / "00_metadata/stage_timing.tsv"
    events = args.output_dir / "00_metadata/workflow_events.tsv"
    console = args.output_dir / "logs/cutnrun2tracks.console.log"
    files = sorted(path for path in args.output_dir.rglob("*") if path.is_file()
                   and ".checkpoints" not in path.parts
                   and path not in {destination, timing, events, console})
    if args.jobs == 1:
        checksums = map(digest, files)
    else:
        executor = ThreadPoolExecutor(max_workers=args.jobs)
        checksums = executor.map(digest, files)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for path, checksum in zip(files, checksums, strict=True):
            handle.write(f"{checksum}  {path.relative_to(args.output_dir).as_posix()}\n")
    if args.jobs != 1:
        executor.shutdown()
    return 0


if __name__ == "__main__": raise SystemExit(main())
