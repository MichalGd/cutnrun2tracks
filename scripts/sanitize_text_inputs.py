#!/usr/bin/env python3
"""Back up and normalize BOM/CRLF workflow inputs to UTF-8/LF."""

from __future__ import annotations

import argparse
import codecs
import os
import shutil
import stat
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def normalize(payload: bytes) -> tuple[bytes, list[str]]:
    findings: list[str] = []
    if payload.startswith((codecs.BOM_UTF32_LE, codecs.BOM_UTF32_BE)):
        text = payload.decode("utf-32")
        findings.append("UTF-32 BOM")
    elif payload.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        text = payload.decode("utf-16")
        findings.append("UTF-16 BOM")
    elif payload.startswith(codecs.BOM_UTF8):
        text = payload.decode("utf-8-sig")
        findings.append("UTF-8 BOM")
    else:
        if b"\x00" in payload:
            raise ValueError("NUL byte without a recognized Unicode BOM")
        text = payload.decode("utf-8")
    if "\r\n" in text:
        findings.append("CRLF")
    text = text.replace("\r\n", "\n")
    if "\r" in text:
        findings.append("bare CR")
        text = text.replace("\r", "\n")
    if text.endswith("\x1a"):
        findings.append("DOS EOF")
        text = text[:-1]
    return text.encode("utf-8"), findings


def sanitize(path: Path) -> None:
    path = path.expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"not a regular file: {path}")
    original = path.read_bytes()
    corrected, findings = normalize(original)
    if not findings:
        print(f"[INPUT] OK: {path}")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    backup = path.with_name(f"{path.name}.windows-artifact-backup.{stamp}.{os.getpid()}")
    shutil.copy2(path, backup)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.sanitize.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(corrected)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, stat.S_IMODE(path.stat().st_mode))
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    print(f"[INPUT] Corrected: {path} ({', '.join(findings)}); backup={backup}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()
    failed = False
    for path in args.files:
        try:
            sanitize(path)
        except (OSError, UnicodeError, ValueError) as exc:
            print(f"[INPUT] ERROR: {path}: {exc}", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
