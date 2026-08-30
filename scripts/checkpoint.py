#!/usr/bin/env python3
"""Content-aware workflow signatures and JSON checkpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expanded(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            files.extend(item for item in path.rglob("*") if item.is_file())
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(path)
    return sorted(set(item.resolve() for item in files), key=str)


def hashes(paths: list[Path], jobs: int) -> list[str]:
    if jobs == 1:
        return [sha256(path) for path in paths]
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        return list(executor.map(sha256, paths))


def signature(paths: list[Path], jobs: int = 1) -> str:
    digest = hashlib.sha256()
    files = expanded(paths)
    for path, checksum in zip(files, hashes(files, jobs), strict=True):
        digest.update(str(path).encode())
        digest.update(b"\0")
        digest.update(checksum.encode())
        digest.update(b"\0")
    return digest.hexdigest()


def snapshot(paths: list[Path], jobs: int = 1) -> list[dict[str, object]]:
    result = []
    files = expanded(paths)
    for path, checksum in zip(files, hashes(files, jobs), strict=True):
        stat = path.stat()
        result.append({"path": str(path), "size": stat.st_size, "sha256": checksum})
    if not result:
        raise ValueError("checkpoint has no output files")
    return result


def command_signature(args: argparse.Namespace) -> int:
    print(signature(args.paths, args.jobs))
    return 0


def command_write(args: argparse.Namespace) -> int:
    payload = {
        "schema": 1,
        "stage": args.stage,
        "signature": args.signature,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": snapshot(args.outputs, args.jobs),
    }
    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.checkpoint.with_suffix(args.checkpoint.suffix + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, args.checkpoint)
    return 0


def command_check(args: argparse.Namespace) -> int:
    try:
        payload = json.loads(args.checkpoint.read_text(encoding="utf-8"))
        if payload.get("signature") != args.signature or payload.get("stage") != args.stage:
            return 1
        outputs = payload.get("outputs", [])
        output_paths = [Path(item["path"]) for item in outputs]
        if any(not path.is_file() for path in output_paths):
            return 1
        for item, path, checksum in zip(outputs, output_paths, hashes(output_paths, args.jobs), strict=True):
            if path.stat().st_size != item["size"] or checksum != item["sha256"]:
                return 1
        return 0 if outputs else 1
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return 1


def command_adopt(args: argparse.Namespace) -> int:
    """Adopt validated prior-stage outputs into a new run signature."""
    try:
        payload = json.loads(args.checkpoint.read_text(encoding="utf-8"))
        if payload.get("stage") != args.stage:
            return 1
        outputs = payload.get("outputs", [])
        if not outputs:
            return 1
        output_paths = [Path(item["path"]) for item in outputs]
        if any(not path.is_file() for path in output_paths):
            return 1
        for item, path, checksum in zip(outputs, output_paths, hashes(output_paths, args.jobs), strict=True):
            if path.stat().st_size != item["size"] or checksum != item["sha256"]:
                return 1
        previous = payload.get("signature")
        adopted_utc = datetime.now(timezone.utc).isoformat()
        payload.setdefault("signature_adoptions", []).append(
            {
                "from": previous,
                "to": args.signature,
                "adopted_utc": adopted_utc,
                "reason": "explicit --from-stage prior-stage reuse",
            }
        )
        payload["signature"] = args.signature
        temporary = args.checkpoint.with_suffix(
            args.checkpoint.suffix + f".tmp.{os.getpid()}"
        )
        temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, args.checkpoint)
        return 0
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return 1


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    sig = sub.add_parser("signature")
    sig.add_argument("paths", nargs="+", type=Path)
    sig.add_argument("--jobs", type=int, default=1)
    sig.set_defaults(function=command_signature)
    write = sub.add_parser("write")
    write.add_argument("--checkpoint", type=Path, required=True)
    write.add_argument("--stage", required=True)
    write.add_argument("--signature", required=True)
    write.add_argument("--outputs", nargs="+", type=Path, required=True)
    write.add_argument("--jobs", type=int, default=1)
    write.set_defaults(function=command_write)
    check = sub.add_parser("check")
    check.add_argument("--checkpoint", type=Path, required=True)
    check.add_argument("--stage", required=True)
    check.add_argument("--signature", required=True)
    check.add_argument("--jobs", type=int, default=1)
    check.set_defaults(function=command_check)
    adopt = sub.add_parser("adopt")
    adopt.add_argument("--checkpoint", type=Path, required=True)
    adopt.add_argument("--stage", required=True)
    adopt.add_argument("--signature", required=True)
    adopt.add_argument("--jobs", type=int, default=1)
    adopt.set_defaults(function=command_adopt)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be a positive integer")
    try:
        return args.function(args)
    except (OSError, ValueError) as exc:
        print(f"CHECKPOINT ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
