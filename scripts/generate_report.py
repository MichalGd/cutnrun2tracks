#!/usr/bin/env python3
"""Generate the dependency-light HTML and TSV workflow report."""

from __future__ import annotations

import argparse
import csv
import html
from datetime import datetime, timezone
from pathlib import Path


def table(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def html_table(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "<p>Not available.</p>"
    columns = list(rows[0])
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns) + "</tr>" for row in rows)
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    root = args.output_dir
    report_dir = root / "10_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    samples = table(root / "00_metadata/sample_manifest.tsv")
    cohorts = table(root / "00_metadata/cohort_manifest.tsv")
    counts = table(root / "06_qc/alignment_and_complexity/observation_counts.tsv")
    consensus = table(root / "05_peaks/consensus/consensus_status.tsv")
    peakcalls = table(root / "05_peaks/per_sample/peakcall_status.tsv")
    normalized = table(root / "04_tracks/normalized_track_family_status.tsv")
    differential = table(root / "08_differential/stage_status.tsv")
    metagene = table(root / "06_qc/metagene/artifacts.tsv")
    warnings = []
    for failure in root.rglob("FAILED.json"):
        warnings.append({"severity": "ERROR", "item": str(failure.relative_to(root))})
    for skipped in root.rglob("SKIPPED.json"):
        warnings.append({"severity": "INFO", "item": str(skipped.relative_to(root))})
    for row in peakcalls:
        if row.get("status") != "SUCCESS" or row.get("caller_warnings") not in {"", "."}:
            warnings.append({
                "severity": "WARNING",
                "item": (f"peakcalling:{row.get('sample_key', '?')}:"
                         f"{row.get('status', '?')}:{row.get('caller_warnings', '.')}")
            })
    for row in consensus:
        if row.get("status") != "SUCCESS":
            warnings.append({"severity": "WARNING", "item": (
                f"consensus:{row.get('cohort_id', '?')}:{row.get('status', '?')}:"
                f"{row.get('reason', '.')}")})
    for row in normalized:
        if row.get("status") != "SUCCESS":
            warnings.append({"severity": "WARNING", "item": (
                f"normalization:{row.get('cohort_id', '?')}:{row.get('policy', '?')}:"
                f"{row.get('status', '?')}:{row.get('reason', '.')}")})
    for row in differential:
        if row.get("status") != "SUCCESS":
            warnings.append({"severity": "WARNING", "item": (
                f"differential:{row.get('status', '?')}:failed_modules="
                f"{row.get('failed_modules', '?')}:skipped_cohorts={row.get('skipped_cohorts', '0')}")})
    preseq_dir = root / "06_qc/alignment_and_complexity"
    for log in sorted((root / "logs/qc").glob("*.preseq.log")):
        sample = log.name.removesuffix(".preseq.log")
        if not (preseq_dir / f"{sample}.preseq.txt").is_file():
            warnings.append({"severity": "WARNING", "item": f"preseq:{sample}:FAILED_OR_EMPTY"})

    warnings.sort(key=lambda row: (row["severity"], row["item"]))
    write_tsv(report_dir / "warning_summary.tsv", warnings, ["severity", "item"])
    summary = [
        {"metric": "biological_libraries", "value": str(len(samples))},
        {"metric": "target_cohorts", "value": str(len(cohorts))},
        {"metric": "consensus_success", "value": str(sum(row.get("status") == "SUCCESS" for row in consensus))},
        {"metric": "consensus_skipped_or_failed", "value": str(sum(row.get("status") != "SUCCESS" for row in consensus))},
        {"metric": "peakcall_samples_with_warnings", "value": str(sum(
            row.get("status") != "SUCCESS" or row.get("caller_warnings") not in {"", "."}
            for row in peakcalls
        ))},
        {"metric": "warnings_and_skips", "value": str(len(warnings))},
        {"metric": "normalized_track_families_skipped_or_failed", "value": str(sum(
            row.get("status") != "SUCCESS" for row in normalized
        ))},
        {"metric": "differential_stage_status", "value": differential[0].get("status", "NOT_AVAILABLE") if differential else "NOT_AVAILABLE"},
        {"metric": "metagene_plot_tasks", "value": str(len(metagene))},
    ]
    write_tsv(report_dir / "run_summary.tsv", summary, ["metric", "value"])
    document = f"""<!doctype html><html><head><meta charset='utf-8'><title>cutnrun2tracks report</title>
<style>body{{font-family:sans-serif;margin:2rem}}table{{border-collapse:collapse;font-size:.85rem}}th,td{{border:1px solid #ccc;padding:.3rem}}th{{background:#eee}}code{{background:#f4f4f4}}</style></head><body>
<h1>cutnrun2tracks report</h1><p>Generated {datetime.now(timezone.utc).isoformat()}.</p>
<p>QC thresholds are descriptive. Different antibody targets are never normalized together.</p>
<h2>Run summary</h2>{html_table(summary)}<h2>Cohorts</h2>{html_table(cohorts)}
<h2>Retained analysis observations</h2>{html_table(counts)}<h2>Consensus status</h2>{html_table(consensus)}
<h2>Per-sample peak-calling status</h2>{html_table(peakcalls)}
<h2>Normalized-track families</h2>{html_table(normalized)}
<h2>Differential analysis</h2>{html_table(differential)}
<h2>Metagene aggregate-signal outputs</h2>{html_table(metagene)}
<h2>Warnings and skips</h2>{html_table(warnings)}</body></html>"""
    temporary_html = report_dir / "pipeline_report.html.tmp"
    temporary_html.write_text(document, encoding="utf-8")
    temporary_html.replace(report_dir / "pipeline_report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
