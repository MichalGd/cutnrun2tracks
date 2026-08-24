#!/usr/bin/env python3
"""Generate a dependency-light HTML and TSV workflow report."""

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
    with (report_dir / "warning_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["severity", "item"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(warnings)
    summary = [
        {"metric": "biological_libraries", "value": str(len(samples))},
        {"metric": "target_cohorts", "value": str(len(cohorts))},
        {"metric": "consensus_success", "value": str(sum(row.get("status") == "SUCCESS" for row in consensus))},
        {"metric": "peakcall_samples_with_warnings", "value": str(sum(
            row.get("status") != "SUCCESS" or row.get("caller_warnings") not in {"", "."}
            for row in peakcalls
        ))},
        {"metric": "warnings_and_skips", "value": str(len(warnings))},
        {"metric": "metagene_plot_tasks", "value": str(len(metagene))},
    ]
    with (report_dir / "run_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary)
    document = f"""<!doctype html><html><head><meta charset='utf-8'><title>cutnrun2tracks report</title>
<style>body{{font-family:sans-serif;margin:2rem}}table{{border-collapse:collapse;font-size:.85rem}}th,td{{border:1px solid #ccc;padding:.3rem}}th{{background:#eee}}code{{background:#f4f4f4}}</style></head><body>
<h1>cutnrun2tracks report</h1><p>Generated {datetime.now(timezone.utc).isoformat()}.</p>
<p>QC thresholds are descriptive. Different antibody targets are never normalized together.</p>
<h2>Run summary</h2>{html_table(summary)}<h2>Cohorts</h2>{html_table(cohorts)}
<h2>Observation counts</h2>{html_table(counts)}<h2>Consensus status</h2>{html_table(consensus)}
<h2>Per-sample peak-calling status</h2>{html_table(peakcalls)}
<h2>Metagene aggregate-signal outputs</h2>{html_table(metagene)}
<h2>Warnings and skips</h2>{html_table(warnings)}</body></html>"""
    (report_dir / "pipeline_report.html").write_text(document, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
