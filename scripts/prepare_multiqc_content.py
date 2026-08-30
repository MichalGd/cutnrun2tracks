#!/usr/bin/env python3
"""Prepare deterministic MultiQC custom content from cutnrun2tracks outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def aggregate(pattern: str, root: Path, source_column: str | None = None) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(root.glob(pattern)):
        for row in read_tsv(path):
            if source_column:
                row[source_column] = path.relative_to(root).as_posix()
            rows.append(row)
    return rows


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "item"
    if len(cleaned) > 160:
        digest = hashlib.sha256(cleaned.encode()).hexdigest()[:12]
        cleaned = f"{cleaned[:147]}_{digest}"
    return cleaned


def write_custom_table(
    destination: Path,
    *,
    identifier: str,
    section_name: str,
    description: str,
    rows: list[dict[str, str]],
    columns: list[str] | None = None,
    key_columns: tuple[str, ...] = (),
) -> bool:
    if not rows:
        return False
    selected = columns or list(rows[0])
    if not selected:
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"# id: {identifier}\n")
        handle.write(f"# section_name: {section_name}\n")
        handle.write(f"# description: {description}\n")
        handle.write("# format: tsv\n")
        handle.write("# plot_type: table\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["row_id", *selected])
        for index, row in enumerate(rows, 1):
            key = " | ".join(row.get(column, "") for column in key_columns).strip(" |")
            writer.writerow([key or f"row_{index}", *(row.get(column, "") for column in selected)])
    return True


def stage_images(root: Path, destination: Path) -> list[tuple[str, str]]:
    patterns = (
        "06_qc/controls/*.target_control_fingerprint.png",
        "06_qc/tss_signal_profile/*.descriptive_TSS_profile.png",
        "06_qc/correlation_pca_fingerprint/**/spearman_heatmap.png",
        "06_qc/correlation_pca_fingerprint/**/pca.png",
        "08_differential/**/pca.png",
        "08_differential/**/dispersion.png",
        "07_annotation/feature_summary/peak_feature_composition.*.png",
    )
    staged: list[tuple[str, str]] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for source in sorted(root.glob(pattern)):
            if not source.is_file() or source.stat().st_size == 0 or source in seen:
                continue
            seen.add(source)
            relative = source.relative_to(root).as_posix()
            name = f"cutnrun2tracks_{safe_id(relative.removesuffix('.png'))}_mqc.png"
            shutil.copy2(source, destination / name)
            staged.append((relative, name))
    return staged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("custom_dir", type=Path)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    custom = args.custom_dir.resolve()
    try:
        if not (root / "00_metadata/sample_manifest.tsv").is_file():
            raise ValueError(f"sample manifest missing under output directory: {root}")
        custom.mkdir(parents=True, exist_ok=True)
        tables: list[tuple[str, str, str, list[dict[str, str]], list[str] | None, tuple[str, ...]]] = [
            ("run_summary", "Run summary", "Run-level cutnrun2tracks result counts.",
             read_tsv(root / "10_reports/run_summary.tsv"), ["value"], ("metric",)),
            ("cohort_membership", "Cohort membership and controls",
             "Auditable target/control membership and shared-control reuse for every derived cohort.",
             read_tsv(root / "00_metadata/cohort_membership.tsv"),
             ["role", "condition", "replicate", "control_key", "control_reused_by_targets"],
             ("cohort_id", "sample_key")),
            ("resource_budget", "Configured CPU budget",
             "Maximum jobs x threads requested by each long-running stage.",
             read_tsv(root / "00_metadata/resource_budget.tsv"),
             ["parallel_jobs", "threads_per_job", "requested_cpus", "budget_cpus", "status"],
             ("stage",)),
            ("warnings", "Workflow warnings and skips", "Recorded sample, cohort, and module warnings.",
             read_tsv(root / "10_reports/warning_summary.tsv"), ["severity", "item"], ("severity", "item")),
            ("observations", "Retained analysis observations", "Post-filter fragments for paired-end libraries or reads for single-end libraries.",
             read_tsv(root / "06_qc/alignment_and_complexity/observation_counts.tsv"),
             ["layout", "signal_unit", "analysis_observations"], ("sample_key",)),
            ("complexity", "Library complexity", "Duplicate-retained q30 library complexity metrics.",
             read_tsv(root / "06_qc/alignment_and_complexity/library_complexity.tsv"),
             ["layout", "total_observations", "distinct_observations", "NRF", "PBC1", "PBC2"],
             ("sample_key",)),
            ("alignment", "Host alignment records", "Host-alignment records before downstream analysis filtering.",
             aggregate("03_alignment/metrics/*.alignment.tsv", root),
             ["layout", "genome", "alignment_records", "spikein_mode"], ("sample_key",)),
            ("filtering", "Filtering and duplicate sensitivity", "Signal-unit counts across the four filtering policies.",
             aggregate("03_alignment/metrics/*.filter_counts.tsv", root),
             ["q0_dup_retained", "q0_dup_removed", "q30_dup_retained", "q30_dup_removed", "analysis_policy"], ("sample_key",)),
            ("peakcalls", "Per-sample peak calling", "Primary peak counts and caller warnings; empty and failed calls remain explicit.",
             read_tsv(root / "05_peaks/per_sample/peakcall_status.tsv"),
             ["primary_caller", "primary_class", "status", "primary_peak_count", "caller_warnings", "reason"], ("sample_key",)),
            ("consensus", "Consensus peak sets", "Successful contributions, exclusions, and retained consensus-region counts by cohort.",
             read_tsv(root / "05_peaks/consensus/consensus_status.tsv"),
             ["status", "total_samples", "successful_peak_samples", "excluded_samples", "regions", "reason"], ("cohort_id",)),
            ("normalization", "Normalized-track families", "Per-cohort normalization availability for each filtering policy.",
             read_tsv(root / "04_tracks/normalized_track_family_status.tsv"),
             ["status", "reason"], ("cohort_id", "policy")),
            ("differential", "Differential enrichment stage", "Run-level differential module failures and cohort skips.",
             read_tsv(root / "08_differential/stage_status.tsv"), None, ("status",)),
            ("frip", "Fraction of signal in consensus peaks", "Descriptive FRiP against each target cohort's consensus set.",
             aggregate("06_qc/frip_and_peak_reproducibility/*.frip.tsv", root),
             ["signal_unit", "total", "in_consensus", "frip"], ("sample_key",)),
            ("peak_features", "Peak genomic-feature composition",
             "Mutually exclusive peak categories; fractions sum to one within each sample/caller peak set.",
             read_tsv(root / "07_annotation/feature_summary/peak_feature_summary.tsv"),
             ["caller", "peak_class", "category", "count", "fraction", "bp_fraction", "total_peaks"],
             ("entity_type", "entity_id", "caller", "peak_class", "category")),
            ("spikein", "Spike-in calibration", "Host/spike observations, scale factors, and calibration status.",
             read_tsv(root / "06_qc/spikein/spikein_scaling.tsv"), None, ("sample_key",)),
            ("comparisons", "Differential comparison summary", "Tested and significant consensus regions for completed DESeq2 comparisons.",
             aggregate("08_differential/**/comparison_summary.tsv", root, "source_table"), None,
             ("source_table", "comparison_id")),
        ]
        manifest_rows: list[tuple[str, str, str]] = []
        for identifier, title, description, rows, columns, keys in tables:
            destination = custom / f"cutnrun2tracks_{identifier}_mqc.tsv"
            if write_custom_table(destination, identifier=f"cutnrun2tracks_{identifier}",
                                  section_name=title, description=description, rows=rows,
                                  columns=columns, key_columns=keys):
                manifest_rows.append(("table", identifier, destination.name))

        for source, destination in stage_images(root, custom):
            manifest_rows.append(("image", source, destination))

        with (custom / "custom_content_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["type", "source", "staged_name"])
            writer.writerows(manifest_rows)
        print(f"Prepared {len(manifest_rows)} MultiQC custom-content items in {custom}")
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"MULTIQC CONTENT ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
