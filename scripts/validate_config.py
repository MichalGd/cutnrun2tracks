#!/usr/bin/env python3
"""Parse and validate the non-executable cutnrun2tracks configuration format."""

from __future__ import annotations

import argparse
import math
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path


BOOLEAN_KEYS = {
    "ALLOW_MIXED_LAYOUTS", "ALLOW_MIXED_GENOMES", "TRIM_ADAPTERS",
    "RUN_FASTQC", "RUN_FASTQC_PER_TECHNICAL_UNIT", "RUN_MULTIQC", "BOWTIE2_DOVETAIL", "BOWTIE2_MIXED",
    "BOWTIE2_DISCORDANT", "CANONICAL_CHROMS_ONLY", "REMOVE_MITO",
    "ALLOW_EMPTY_FILTERED_BAM", "MARK_DUPLICATES",
    "GENERATE_DUPLICATE_RETAINED_BAMS", "GENERATE_DUPLICATE_REMOVED_BAMS",
    "REQUIRE_MATCHED_CONTROL", "ALLOW_SHARED_CONTROLS",
    "ALLOW_CONTROL_FREE_PEAKCALL", "ALLOW_EMPTY_PEAKS", "SEACR_ALLOW_SE",
    "MACS3_CALL_SUMMITS", "ALLOW_SINGLE_SAMPLE_CONSENSUS",
    "CONSENSUS_USE_PRIMARY_CALLER", "REQUIRE_ALL_ENABLED_TRACKS",
    "GENERATE_CPM_TRACKS", "GENERATE_DESEQ2_CONSENSUS_TRACKS",
    "GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS",
    "GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS",
    "GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS",
    "GENERATE_COVERAGE_BIGWIGS", "GENERATE_COVERAGE_BEDGRAPHS",
    "TRACK_STANDARD_CHROMS_ONLY", "ALLOW_FAILED_SPIKEIN",
    "GENERATE_SPIKEIN_CONTROL_TRACKS", "RUN_FRAGMENT_QC",
    "RUN_TSS_SIGNAL_PROFILE", "RUN_ATAQV_QC", "GENERATE_ATAQV_VIEWER",
    "RUN_PRESEQ", "RUN_LIBRARY_COMPLEXITY", "RUN_CROSS_CORRELATION",
    "RUN_REPLICATE_CORRELATION", "RUN_DIFFBIND", "RUN_DESEQ2_ENRICHMENT",
    "RUN_METAGENE", "METAGENE_ALLOW_CPM_FALLBACK",
    "METAGENE_INCLUDE_CONTROLS", "METAGENE_SKIP_ZERO_REGIONS",
    "DIFFERENTIAL_SUBTRACT_CONTROL", "RUN_CONTROL_SUBTRACTED_SENSITIVITY",
    "RUN_TARGET_CONTROL_INTERACTION", "REQUIRE_CONDITION_MATCHED_CONTROLS",
    "RUN_SIMPLE_PEAK_ANNOTATION", "RUN_CCRE_ANNOTATION",
    "RUN_FEATURE_ANNOTATION_SUMMARY", "PEAK_ANNOTATION_INCLUDE_CONSENSUS",
    "RUN_MOTIF_ENRICHMENT", "ENABLE_AUTOMATIC_CLEANUP",
    "KEEP_TRIMMED_FASTQ", "KEEP_RAW_ALIGNMENT_BAMS", "KEEP_MARKED_BAMS",
    "KEEP_FILTERED_BAMS", "KEEP_POLICY_BAMS", "KEEP_RAW_BEDGRAPH",
    "KEEP_SPIKEIN_BAMS", "WRITE_IGV_SESSION", "WRITE_COMMAND_LOG",
    "WRITE_STRUCTURED_LOG", "WRITE_CONSOLE_LOG",
    "WRITE_FILE_CHECKSUMS", "REDACT_PATHS_IN_REPORT",
}

POSITIVE_INTEGER_KEYS = {
    "MIN_TRIMMED_LENGTH", "BOWTIE2_MIN_INSERT", "BOWTIE2_MAX_INSERT",
    "THREADS_FASTQC", "THREADS_TRIMGALORE",
    "THREADS_BOWTIE2", "THREADS_SAMTOOLS", "MIN_MAPQ",
    "CONSENSUS_MIN_BIOLOGICAL_SAMPLES", "SEACR_MAX_FRAGMENT",
    "MACS3_CUTRUN_SE_EXTSIZE", "MACS3_CUTTAG_SE_EXTSIZE",
    "EPIC2_BIN_SIZE", "EPIC2_GAP_SIZE", "EPIC2_FRAGMENT_SIZE", "TRACK_BIN_SIZE",
    "THREADS_BAMCOVERAGE", "SPIKEIN_MIN_MAPQ", "SPIKEIN_MIN_OBSERVATIONS_FAIL",
    "SPIKEIN_MIN_OBSERVATIONS_WARN", "FRAGMENT_PLOT_MAX_BP",
    "TSS_PROFILE_UPSTREAM", "TSS_PROFILE_DOWNSTREAM",
    "DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION", "PEAK_ANNOTATION_PROMOTER_UPSTREAM",
    "PEAK_ANNOTATION_PROMOTER_DOWNSTREAM", "THREADS_PARALLEL_JOBS",
    "QC_SAMPLE_PARALLEL_JOBS", "TRACK_PARALLEL_JOBS",
    "PEAKCALL_PARALLEL_JOBS", "MERGE_PARALLEL_JOBS", "SPIKEIN_PARALLEL_JOBS",
    "NORMALIZED_TRACK_PARALLEL_JOBS", "DIFFERENTIAL_PARALLEL_JOBS",
    "ANNOTATION_PARALLEL_JOBS", "CHECKPOINT_PARALLEL_JOBS", "CHECKSUM_PARALLEL_JOBS",
    "ATAQV_PARALLEL_JOBS", "THREADS_ATAQV", "ATAQV_TSS_EXTENSION",
    "METAGENE_BODY_LENGTH_BP", "METAGENE_BIN_SIZE_BP", "METAGENE_DPI",
    "METAGENE_THREADS_COMPUTEMATRIX", "METAGENE_PARALLEL_JOBS",
}

NONNEGATIVE_INTEGER_KEYS = {
    "BOWTIE2_SEED", "PERMISSIVE_MIN_MAPQ", "INTERMEDIATE_MIN_MAPQ",
    "METAGENE_REFERENCE_UPSTREAM_BP", "METAGENE_REFERENCE_DOWNSTREAM_BP",
    "METAGENE_BODY_UPSTREAM_BP", "METAGENE_BODY_DOWNSTREAM_BP",
    "PEAK_ANNOTATION_GENE_END_WINDOW",
}

FLOAT_KEYS = {
    "SEACR_NO_CONTROL_THRESHOLD", "MACS3_QVALUE", "MACS3_BROAD_CUTOFF", "EPIC2_FDR",
    "SPIKEIN_SCALE_TARGET", "SPIKEIN_WARN_LOW_FRACTION",
    "SPIKEIN_WARN_HIGH_FRACTION", "DIFFERENTIAL_ALPHA",
    "DIFFERENTIAL_MIN_ABS_LOG2FC",
}

ENUMS = {
    "BOWTIE2_MODE": {"end-to-end", "local"},
    "BOWTIE2_PRESET": {"very-sensitive", "sensitive", "very-sensitive-local", "sensitive-local"},
    "BOWTIE2_REPORTING": {"best"},
    "TARGET_DEFAULT_DUPLICATE_POLICY": {"retain", "remove"},
    "CONTROL_DEFAULT_DUPLICATE_POLICY": {"retain", "remove"},
    "SPIKEIN_DUPLICATE_POLICY": {"retain", "remove"},
    "PRIMARY_PEAK_CALLER": {"auto", "seacr", "macs3", "epic2"},
    "PEAKCALL_FAILURE_POLICY": {"fail", "continue"},
    "SEACR_MODE": {"stringent", "relaxed"},
    "SEACR_CONTROL_NORMALIZATION": {"norm", "non"},
    "MACS3_KEEP_DUP": {"all", "auto"},
    "CONSENSUS_PEAK_CLASS": {"auto", "narrow", "broad"},
    "SE_SIGNAL_MODE": {"read"},
    "SPIKEIN_MODE": {"none", "dm6", "ecoli", "custom"},
    "QC_THRESHOLDS_MODE": {"descriptive"},
    "METAGENE_TRACK_FAMILY": {"auto", "cpm", "spikein"},
    "METAGENE_MISSING_DATA_POLICY": {"zero", "na"},
    "DIFFERENTIAL_NORMALIZATION": {"deseq2", "spikein"},
    "DIFFERENTIAL_CONTROL_MODE": {"peak_calling_only", "control_subtracted", "interaction"},
    "CHECKPOINT_MODE": {"signature_and_outputs"},
    "RESOURCE_CHECK_MODE": {"warn", "fail"},
}

REFERENCE_KEY = re.compile(
    r"^(INDEX|FASTA|CHROM_SIZES|CANONICAL_CONTIGS|GTF|BLACKLIST|EFFECTIVE_GENOME_SIZE|TSS_BED|CCRE_BED)_[A-Z0-9_]+$"
)
KEY_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")


def template_keys(template: Path) -> set[str]:
    keys: set[str] = set()
    for raw in template.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = re.match(r"^([A-Z][A-Z0-9_]*)\s*=", line)
        if match:
            keys.add(match.group(1))
    return keys


def parse_value(raw: str, line_number: int) -> str:
    if any(token in raw for token in ("$(", "${", "`")):
        raise ValueError(f"line {line_number}: shell expansion is prohibited")
    lexer = shlex.shlex(raw, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = "#"
    parts = list(lexer)
    if not parts:
        return ""
    if len(parts) != 1:
        raise ValueError(
            f"line {line_number}: values containing spaces must be quoted"
        )
    return parts[0]


def parse_config(path: Path, allowed: set[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            raise ValueError(f"line {number}: use KEY=VALUE without export")
        if "=" not in line:
            raise ValueError(f"line {number}: expected KEY=VALUE")
        key, raw_value = line.split("=", 1)
        key = key.strip()
        if not KEY_RE.fullmatch(key):
            raise ValueError(f"line {number}: invalid key {key!r}")
        if key not in allowed and not REFERENCE_KEY.fullmatch(key):
            raise ValueError(f"line {number}: unknown configuration key {key}")
        if key in values:
            raise ValueError(f"line {number}: duplicate configuration key {key}")
        values[key] = parse_value(raw_value.strip(), number)
    return values


def validate(values: dict[str, str], required_keys: set[str]) -> list[str]:
    errors: list[str] = []
    missing = sorted(required_keys - values.keys())
    if missing:
        errors.append("missing template keys: " + ", ".join(missing))

    for key in BOOLEAN_KEYS:
        if key in values and values[key].lower() not in {"true", "false"}:
            errors.append(f"{key} must be true or false")
        elif key in values:
            values[key] = values[key].lower()
    for key in POSITIVE_INTEGER_KEYS:
        if key in values:
            try:
                if int(values[key]) <= 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{key} must be a positive integer")
    for key in NONNEGATIVE_INTEGER_KEYS:
        if key in values:
            try:
                if int(values[key]) < 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{key} must be a non-negative integer")
    for key in FLOAT_KEYS:
        if key in values:
            try:
                number = float(values[key])
                if not math.isfinite(number) or number < 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{key} must be a finite non-negative number")
    for key, choices in ENUMS.items():
        if key in values and values[key].lower() not in choices:
            errors.append(f"{key} must be one of: {', '.join(sorted(choices))}")
        elif key in values:
            values[key] = values[key].lower()

    callers = [item.strip().lower() for item in values.get("PEAK_CALLERS", "").split(",") if item.strip()]
    if not callers or not set(callers) <= {"seacr", "macs3", "epic2"}:
        errors.append("PEAK_CALLERS must be a comma-separated subset of seacr,macs3,epic2")
    values["PEAK_CALLERS"] = ",".join(dict.fromkeys(callers))
    if values.get("PRIMARY_PEAK_CALLER") not in {"auto", *callers}:
        errors.append("PRIMARY_PEAK_CALLER must be auto or an enabled PEAK_CALLERS value")

    if values.get("BOWTIE2_MAX_INSERT", "0").isdigit() and values.get("BOWTIE2_MIN_INSERT", "0").isdigit():
        if int(values["BOWTIE2_MAX_INSERT"]) < int(values["BOWTIE2_MIN_INSERT"]):
            errors.append("BOWTIE2_MAX_INSERT must be >= BOWTIE2_MIN_INSERT")
    if values.get("SPIKEIN_MODE") != "none":
        for key in ("SPIKEIN_REFERENCE_ID", "SPIKEIN_INDEX", "SPIKEIN_FASTA", "SPIKEIN_CHROM_SIZES", "SPIKEIN_ALLOWED_CONTIGS"):
            if not values.get(key):
                errors.append(f"{key} is required when SPIKEIN_MODE is enabled")
    if values.get("DIFFERENTIAL_NORMALIZATION") == "spikein" and values.get("SPIKEIN_MODE") == "none":
        errors.append("DIFFERENTIAL_NORMALIZATION=spikein requires SPIKEIN_MODE")
    for key in ("SEACR_NO_CONTROL_THRESHOLD", "MACS3_QVALUE", "MACS3_BROAD_CUTOFF", "EPIC2_FDR",
                "SPIKEIN_WARN_LOW_FRACTION", "SPIKEIN_WARN_HIGH_FRACTION", "DIFFERENTIAL_ALPHA"):
        if key in values:
            try:
                if not 0 < float(values[key]) <= 1:
                    errors.append(f"{key} must be >0 and <=1")
            except ValueError:
                pass
    if values.get("DIFFERENTIAL_SUBTRACT_CONTROL") == "true":
        errors.append("DIFFERENTIAL_SUBTRACT_CONTROL cannot alter the primary model; use RUN_CONTROL_SUBTRACTED_SENSITIVITY=true")
    total_cpu_budget = values.get("TOTAL_CPU_BUDGET", "")
    if total_cpu_budget != "auto":
        try:
            if int(total_cpu_budget) <= 0:
                raise ValueError
        except ValueError:
            errors.append("TOTAL_CPU_BUDGET must be auto or a positive integer")
    precedence = [item.strip().lower() for item in values.get("PEAK_ANNOTATION_FEATURE_PRECEDENCE", "").split(",") if item.strip()]
    expected_categories = {"promoter", "enhancer", "exon", "intron", "gene_end", "other_regulatory", "intergenic", "unclassified"}
    if len(precedence) != len(expected_categories) or set(precedence) != expected_categories:
        errors.append("PEAK_ANNOTATION_FEATURE_PRECEDENCE must contain each supported category exactly once")
    values["PEAK_ANNOTATION_FEATURE_PRECEDENCE"] = ",".join(precedence)
    annotation_formats = [item.strip().lower() for item in values.get("PEAK_ANNOTATION_PLOT_FORMATS", "").split(",") if item.strip()]
    if not annotation_formats or not set(annotation_formats) <= {"png", "pdf", "svg"}:
        errors.append("PEAK_ANNOTATION_PLOT_FORMATS must be a nonempty subset of png,pdf,svg")
    values["PEAK_ANNOTATION_PLOT_FORMATS"] = ",".join(dict.fromkeys(annotation_formats))
    gene_sets = [item.strip() for item in values.get("METAGENE_GENE_SETS", "").split(",") if item.strip()]
    if not gene_sets or any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", item) for item in gene_sets):
        errors.append("METAGENE_GENE_SETS must be a nonempty comma-separated list of safe IDs")
    values["METAGENE_GENE_SETS"] = ",".join(dict.fromkeys(gene_sets))
    modes = [item.strip().lower() for item in values.get("METAGENE_MODES", "").split(",") if item.strip()]
    if not modes or not set(modes) <= {"tss", "tes", "gene_body"}:
        errors.append("METAGENE_MODES must be a comma-separated subset of tss,tes,gene_body")
    values["METAGENE_MODES"] = ",".join(dict.fromkeys(modes))
    formats = [item.strip().lower() for item in values.get("METAGENE_PLOT_FORMATS", "").split(",") if item.strip()]
    if not {"png", "pdf"} <= set(formats) or not set(formats) <= {"png", "pdf"}:
        errors.append("METAGENE_PLOT_FORMATS must contain png and pdf only")
    values["METAGENE_PLOT_FORMATS"] = ",".join(dict.fromkeys(formats))
    for key in ("METAGENE_ZMIN", "METAGENE_ZMAX"):
        value = values.get(key, "")
        if value != "auto":
            try:
                if not math.isfinite(float(value)):
                    raise ValueError
            except ValueError:
                errors.append(f"{key} must be auto or a finite number")
    if values.get("RUN_METAGENE") == "true":
        if not values.get("METAGENE_GENE_SET_MANIFEST"):
            errors.append("METAGENE_GENE_SET_MANIFEST is required when RUN_METAGENE=true")
        if values.get("GENERATE_COVERAGE_BIGWIGS") != "true":
            errors.append("RUN_METAGENE=true requires GENERATE_COVERAGE_BIGWIGS=true")
        family = values.get("METAGENE_TRACK_FAMILY")
        if family == "cpm" and values.get("GENERATE_CPM_TRACKS") != "true":
            errors.append("METAGENE_TRACK_FAMILY=cpm requires GENERATE_CPM_TRACKS=true")
        if family == "spikein" and values.get("SPIKEIN_MODE") == "none":
            errors.append("METAGENE_TRACK_FAMILY=spikein requires SPIKEIN_MODE")
        if family == "auto" and values.get("SPIKEIN_MODE") == "none" and values.get("GENERATE_CPM_TRACKS") != "true":
            errors.append("RUN_METAGENE auto track selection requires CPM tracks when spike-in is disabled")
    for required in ("SAMPLESHEET", "OUTPUT_DIR"):
        value = values.get(required, "")
        if not value or value.startswith("/absolute/path/to/"):
            errors.append(f"{required} must be set to a real path")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--samplesheet")
    parser.add_argument("--output-dir")
    parser.add_argument("--write-shell", type=Path, required=True)
    parser.add_argument("--write-tsv", type=Path)
    args = parser.parse_args()

    try:
        allowed = template_keys(args.template)
        values = parse_config(args.config, allowed)
        if args.samplesheet:
            values["SAMPLESHEET"] = args.samplesheet
        if args.output_dir:
            values["OUTPUT_DIR"] = args.output_dir
        if values.get("RUN_ID") == "auto":
            values["RUN_ID"] = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        errors = validate(values, allowed)
    except (OSError, ValueError) as exc:
        print(f"CONFIG ERROR: {exc}", file=sys.stderr)
        return 1
    if errors:
        for error in errors:
            print(f"CONFIG ERROR: {error}", file=sys.stderr)
        return 1

    args.write_shell.parent.mkdir(parents=True, exist_ok=True)
    with args.write_shell.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# Generated by validate_config.py; do not edit.\n")
        for key in sorted(values):
            handle.write(f"{key}={shlex.quote(values[key])}\n")
    if args.write_tsv:
        args.write_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.write_tsv.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("key\tvalue\n")
            for key in sorted(values):
                handle.write(f"{key}\t{values[key]}\n")
    print(f"Validated configuration: {args.config}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
