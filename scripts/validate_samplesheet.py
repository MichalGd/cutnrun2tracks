#!/usr/bin/env python3
"""Validate CUT&RUN/CUT&Tag metadata and emit biological/cohort manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path


COLUMNS = [
    "sample_id", "fastq_1", "fastq_2", "layout", "genome",
    "assay_profile", "factor", "antibody_id", "target_class", "condition",
    "treatment", "cell_type", "replicate", "tech_replicate", "is_control",
    "control_type", "control_id", "analysis_duplicate_policy",
    "spikein_to_host_ratio", "spikein_stage", "spikein_lot", "batch",
    "donor", "output_prefix",
]
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
TRUE_VALUES = {"true", "1", "yes"}
FALSE_VALUES = {"false", "0", "no"}
CONTROL_TYPES = {"igg", "input", "mock"}
TARGET_CLASSES = {"narrow", "broad", "mixed"}
SPIKE_STAGES = {"cells", "nuclei", "chromatin", "post_library"}


def boolean(value: str, row_number: int) -> bool:
    normalized = value.strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ValueError(f"row {row_number}: is_control must be TRUE or FALSE")


def positive_integer(value: str, field: str, row_number: int) -> int:
    try:
        result = int(value)
    except ValueError as exc:
        raise ValueError(f"row {row_number}: {field} must be a positive integer") from exc
    if result <= 0:
        raise ValueError(f"row {row_number}: {field} must be a positive integer")
    return result


def checked_id(value: str, field: str, row_number: int) -> str:
    if not SAFE_ID.fullmatch(value):
        raise ValueError(f"row {row_number}: unsafe or empty {field}: {value!r}")
    return value


def cohort_identity(row: dict[str, str], primary_caller: str, primary_class: str,
                    spike_mode: str, spike_reference: str) -> tuple[str, ...]:
    identity = (
        row["genome"], row["assay_profile"], row["factor"], row["antibody_id"],
        row["layout"], row["target_class"], row["analysis_duplicate_policy"],
        primary_caller, primary_class,
    )
    if spike_mode != "none":
        identity += (spike_mode, spike_reference, row["spikein_stage"], row["spikein_lot"])
    return identity


def cohort_slug(identity: tuple[str, ...]) -> str:
    readable = "__".join(identity[:5])
    readable = re.sub(r"[^A-Za-z0-9._-]+", "-", readable).strip("-")[:96]
    digest = hashlib.sha256("\x1f".join(identity).encode()).hexdigest()[:8]
    return f"{readable}__{digest}"


def resolve_primary(row: dict[str, str], configured: str, callers: set[str]) -> tuple[str, str]:
    if configured != "auto":
        caller = configured
    elif row["target_class"] in {"broad", "mixed"} and "epic2" in callers:
        caller = "epic2"
    elif row["layout"] == "PE" and row["target_class"] == "narrow" and "seacr" in callers:
        caller = "seacr"
    else:
        caller = "macs3"
    if caller not in callers:
        raise ValueError(f"primary caller {caller} is not enabled in PEAK_CALLERS")
    if caller == "epic2" and row["target_class"] == "narrow":
        raise ValueError("epic2 cannot be the primary caller for target_class=narrow")
    peak_class = "narrow" if row["target_class"] == "narrow" else "broad"
    return caller, peak_class


def compatible_control(target: dict[str, str], control: dict[str, str], shared: bool) -> list[str]:
    required = [
        "genome", "assay_profile", "layout", "treatment", "cell_type", "batch",
        "spikein_to_host_ratio", "spikein_stage", "spikein_lot",
    ]
    if not shared:
        required.append("condition")
    return [field for field in required if target[field] != control[field]]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: ("." if row.get(key, "") == "" else row.get(key, "."))
                for key in fieldnames
            })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("samplesheet", type=Path)
    parser.add_argument("--spikein-mode", choices=["none", "dm6", "ecoli", "custom"], required=True)
    parser.add_argument("--spikein-reference-id", default="")
    parser.add_argument("--peak-callers", default="seacr,macs3")
    parser.add_argument("--primary-peak-caller", choices=["auto", "seacr", "macs3", "epic2"], default="auto")
    parser.add_argument("--blacklist-map", action="append", default=[], metavar="GENOME=PATH")
    parser.add_argument("--allow-shared-controls", action="store_true")
    parser.add_argument("--allow-mixed-layouts", action="store_true")
    parser.add_argument("--allow-mixed-genomes", action="store_true")
    parser.add_argument("--allow-control-free", action="store_true")
    parser.add_argument("--check-files", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    blacklist_map: dict[str, str] = {}
    try:
        for item in args.blacklist_map:
            genome, separator, path = item.partition("=")
            if not separator or not genome or not path:
                raise ValueError(f"invalid --blacklist-map value: {item!r}; expected GENOME=PATH")
            blacklist_map[genome] = path
    except ValueError as exc:
        print(f"SAMPLESHEET ERROR: {exc}", file=sys.stderr)
        return 1

    errors: list[str] = []
    try:
        with args.samplesheet.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != COLUMNS:
                missing = [name for name in COLUMNS if name not in (reader.fieldnames or [])]
                extra = [name for name in (reader.fieldnames or []) if name not in COLUMNS]
                raise ValueError(f"header must exactly match template; missing={missing}, extra={extra}")
            raw_rows = list(reader)
    except (OSError, ValueError) as exc:
        print(f"SAMPLESHEET ERROR: {exc}", file=sys.stderr)
        return 1
    if not raw_rows:
        print("SAMPLESHEET ERROR: no data rows", file=sys.stderr)
        return 1

    rows: list[dict[str, str]] = []
    seen_units: set[tuple[str, int, int]] = set()
    seen_exact: set[tuple[str, ...]] = set()
    fastq_paths: set[str] = set()
    output_owner: dict[str, tuple[str, int]] = {}
    groups: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    layouts: set[str] = set()
    genomes: set[str] = set()

    for number, raw in enumerate(raw_rows, 2):
        row = {key: (value or "").strip() for key, value in raw.items()}
        try:
            sid = checked_id(row["sample_id"], "sample_id", number)
            prefix = checked_id(row["output_prefix"], "output_prefix", number)
            checked_id(row["factor"], "factor", number)
            checked_id(row["antibody_id"], "antibody_id", number)
            checked_id(row["condition"], "condition", number)
            replicate = positive_integer(row["replicate"], "replicate", number)
            tech = positive_integer(row["tech_replicate"], "tech_replicate", number)
            row["replicate"] = str(replicate)
            row["tech_replicate"] = str(tech)
            layout = row["layout"].upper()
            if layout not in {"PE", "SE"}:
                raise ValueError(f"row {number}: layout must be PE or SE")
            row["layout"] = layout
            row["assay_profile"] = row["assay_profile"].lower()
            if row["assay_profile"] not in {"cutrun", "cuttag"}:
                raise ValueError(f"row {number}: assay_profile must be cutrun or cuttag")
            if not row["genome"]:
                raise ValueError(f"row {number}: genome is empty")
            if not row["fastq_1"]:
                raise ValueError(f"row {number}: fastq_1 is empty")
            if layout == "PE" and not row["fastq_2"]:
                raise ValueError(f"row {number}: PE sample requires fastq_2")
            if layout == "SE" and row["fastq_2"]:
                raise ValueError(f"row {number}: SE sample must have empty fastq_2")
            is_control = boolean(row["is_control"], number)
            row["is_control"] = "TRUE" if is_control else "FALSE"
            row["control_type"] = row["control_type"].lower()
            row["target_class"] = row["target_class"].lower()
            row["analysis_duplicate_policy"] = row["analysis_duplicate_policy"].lower()
            if row["analysis_duplicate_policy"] not in {"retain", "remove"}:
                raise ValueError(f"row {number}: analysis_duplicate_policy must be retain or remove")
            if is_control:
                if row["target_class"] != "control" or row["control_type"] not in CONTROL_TYPES:
                    raise ValueError(f"row {number}: control needs target_class=control and igg/input/mock type")
                if row["control_id"]:
                    raise ValueError(f"row {number}: control_id must be empty on a control row")
            else:
                if row["target_class"] not in TARGET_CLASSES or row["control_type"] != "none":
                    raise ValueError(f"row {number}: target needs narrow/broad/mixed class and control_type=none")
                if not row["control_id"] and not args.allow_control_free:
                    raise ValueError(f"row {number}: target requires control_id")
            row["blacklist"] = blacklist_map.get(row["genome"], "")
            if not row["blacklist"]:
                raise ValueError(
                    f"row {number}: no configured blacklist reference for genome {row['genome']}"
                )
            if any("\t" in value or ";" in value for value in row.values()):
                raise ValueError(f"row {number}: tabs and semicolons are not allowed in fields")
            if args.spikein_mode == "none":
                if row["spikein_to_host_ratio"] or row["spikein_stage"] or row["spikein_lot"]:
                    raise ValueError(f"row {number}: spike metadata must be empty for SPIKEIN_MODE=none")
            else:
                try:
                    ratio = float(row["spikein_to_host_ratio"])
                    if ratio <= 0:
                        raise ValueError
                except ValueError as exc:
                    raise ValueError(f"row {number}: spikein_to_host_ratio must be >0") from exc
                if row["spikein_stage"] not in SPIKE_STAGES or not row["spikein_lot"]:
                    raise ValueError(f"row {number}: valid spikein_stage and spikein_lot are required")
            unit = (sid, replicate, tech)
            exact = tuple(row[column] for column in COLUMNS)
            if unit in seen_units:
                raise ValueError(f"row {number}: duplicate sample/replicate/tech_replicate key")
            if exact in seen_exact:
                raise ValueError(f"row {number}: exact duplicate row")
            seen_units.add(unit)
            seen_exact.add(exact)
            owner = (sid, replicate)
            if prefix in output_owner and output_owner[prefix] != owner:
                raise ValueError(f"row {number}: output_prefix is assigned to multiple biological samples")
            output_owner[prefix] = owner
            for field in ("fastq_1", "fastq_2"):
                path = row[field]
                if path:
                    if path in fastq_paths:
                        raise ValueError(f"row {number}: FASTQ reused: {path}")
                    fastq_paths.add(path)
                    if args.check_files and not Path(path).is_file():
                        raise ValueError(f"row {number}: file not found: {path}")
            if args.check_files and not Path(row["blacklist"]).is_file():
                raise ValueError(f"row {number}: blacklist not found: {row['blacklist']}")
            layouts.add(layout)
            genomes.add(row["genome"])
            row["_row"] = str(number)
            groups[(sid, replicate)].append(row)
            rows.append(row)
        except ValueError as exc:
            errors.append(str(exc))

    if len(layouts) > 1 and not args.allow_mixed_layouts:
        errors.append("mixed PE/SE runs require ALLOW_MIXED_LAYOUTS=true")
    if len(genomes) > 1 and not args.allow_mixed_genomes:
        errors.append("mixed genomes require ALLOW_MIXED_GENOMES=true")
    if "SE" in layouts and "seacr" in args.peak_callers.split(","):
        errors.append("SEACR is PE-only in cutnrun2tracks v1; remove seacr from PEAK_CALLERS")

    consistent = [column for column in COLUMNS if column not in {"fastq_1", "fastq_2", "tech_replicate"}]
    for group, members in groups.items():
        first = members[0]
        for member in members[1:]:
            changed = [column for column in consistent if member[column] != first[column]]
            if changed:
                errors.append(f"rows for biological key {group} disagree in: {', '.join(changed)}")

    controls_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for members in groups.values():
        row = members[0]
        if row["is_control"] == "TRUE":
            controls_by_id[row["sample_id"]].append(row)

    biological_rows: list[dict[str, object]] = []
    callers = {item.strip() for item in args.peak_callers.split(",") if item.strip()}
    for (sid, replicate), members in sorted(groups.items()):
        row = dict(members[0])
        row.pop("_row", None)
        row["sample_key"] = f"{row['output_prefix']}.bioR{replicate}"
        row["fastq_1_list"] = ";".join(member["fastq_1"] for member in sorted(members, key=lambda x: int(x["tech_replicate"])))
        row["fastq_2_list"] = ";".join(member["fastq_2"] for member in sorted(members, key=lambda x: int(x["tech_replicate"])) if member["fastq_2"])
        row["technical_units"] = len(members)
        row["control_key"] = "."
        row["cohort_id"] = ""
        row["cohort_key"] = ""
        row["primary_peak_caller"] = "none"
        row["primary_peak_class"] = "control"
        if row["is_control"] == "FALSE":
            try:
                caller, peak_class = resolve_primary(row, args.primary_peak_caller, callers)
                row["primary_peak_caller"] = caller
                row["primary_peak_class"] = peak_class
                identity = cohort_identity(row, caller, peak_class, args.spikein_mode, args.spikein_reference_id)
                row["cohort_key"] = "|".join(identity)
                row["cohort_id"] = cohort_slug(identity)
            except ValueError as exc:
                errors.append(f"sample {sid} replicate {replicate}: {exc}")
            control_id = row["control_id"]
            if control_id:
                candidates = controls_by_id.get(control_id, [])
                exact = [candidate for candidate in candidates if int(candidate["replicate"]) == replicate]
                shared = False
                selected = exact
                if not selected and args.allow_shared_controls:
                    selected = [candidate for candidate in candidates if candidate["condition"].lower() == "shared"]
                    shared = True
                if len(selected) != 1:
                    errors.append(
                        f"sample {sid} replicate {replicate}: expected exactly one matched control {control_id}; found {len(selected)}"
                    )
                else:
                    mismatches = compatible_control(row, selected[0], shared)
                    if mismatches:
                        errors.append(
                            f"sample {sid} replicate {replicate}: control {control_id} incompatible in {','.join(mismatches)}"
                        )
                    else:
                        row["control_key"] = f"{selected[0]['output_prefix']}.bioR{selected[0]['replicate']}"
        biological_rows.append(row)

    control_users: dict[str, list[str]] = defaultdict(list)
    for row in biological_rows:
        if row["is_control"] == "FALSE" and row["control_key"] not in {"", "."}:
            control_users[str(row["control_key"])].append(str(row["sample_key"]))
    if not args.allow_shared_controls:
        for control_key, users in control_users.items():
            if len(users) > 1:
                errors.append(
                    f"control {control_key} is assigned to multiple targets ({','.join(users)}); "
                    "set ALLOW_SHARED_CONTROLS=true after confirming compatibility"
                )

    if errors:
        for error in errors:
            print(f"SAMPLESHEET ERROR: {error}", file=sys.stderr)
        return 1

    manifest_fields = [
        "sample_key", "sample_id", "replicate", "layout", "genome", "assay_profile",
        "factor", "antibody_id", "target_class", "condition", "treatment", "cell_type",
        "is_control", "control_type", "control_id", "control_key",
        "analysis_duplicate_policy", "blacklist", "spikein_to_host_ratio",
        "spikein_stage", "spikein_lot", "batch", "donor", "output_prefix",
        "technical_units", "fastq_1_list", "fastq_2_list", "cohort_id", "cohort_key",
        "primary_peak_caller", "primary_peak_class",
    ]
    manifest = args.output_dir / "sample_manifest.tsv"
    write_tsv(manifest, manifest_fields, biological_rows)

    target_rows = [row for row in biological_rows if row["is_control"] == "FALSE"]
    by_cohort: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in target_rows:
        by_cohort[str(row["cohort_id"])].append(row)
    cohort_rows: list[dict[str, object]] = []
    for cohort_id, members in sorted(by_cohort.items()):
        first = members[0]
        conditions = list(dict.fromkeys(str(member["condition"]) for member in members))
        cohort_rows.append({
            "cohort_id": cohort_id,
            "cohort_key": first["cohort_key"],
            "genome": first["genome"],
            "assay_profile": first["assay_profile"],
            "factor": first["factor"],
            "antibody_id": first["antibody_id"],
            "layout": first["layout"],
            "target_class": first["target_class"],
            "analysis_duplicate_policy": first["analysis_duplicate_policy"],
            "primary_peak_caller": first["primary_peak_caller"],
            "primary_peak_class": first["primary_peak_class"],
            "biological_samples": len(members),
            "sample_keys": ",".join(str(member["sample_key"]) for member in members),
            "conditions": ",".join(conditions),
        })
    cohort_fields = [
        "cohort_id", "cohort_key", "genome", "assay_profile", "factor",
        "antibody_id", "layout", "target_class", "analysis_duplicate_policy",
        "primary_peak_caller", "primary_peak_class", "biological_samples",
        "sample_keys", "conditions",
    ]
    write_tsv(args.output_dir / "cohort_manifest.tsv", cohort_fields, cohort_rows)
    membership_rows: list[dict[str, object]] = []
    biological_by_key = {str(row["sample_key"]): row for row in biological_rows}
    for row in target_rows:
        membership_rows.append({
            "cohort_id": row["cohort_id"], "sample_key": row["sample_key"],
            "role": "target", "condition": row["condition"], "replicate": row["replicate"],
            "control_key": row["control_key"], "control_reused_by_targets": ".",
        })
    for control_key, users in sorted(control_users.items()):
        control = biological_by_key.get(control_key)
        membership_rows.append({
            "cohort_id": ",".join(sorted({str(biological_by_key[user]["cohort_id"]) for user in users})),
            "sample_key": control_key, "role": "control",
            "condition": control["condition"] if control else ".",
            "replicate": control["replicate"] if control else ".",
            "control_key": ".", "control_reused_by_targets": ",".join(users),
        })
    write_tsv(args.output_dir / "cohort_membership.tsv", [
        "cohort_id", "sample_key", "role", "condition", "replicate", "control_key",
        "control_reused_by_targets",
    ], membership_rows)
    print(
        f"Validated {len(raw_rows)} sequencing units, {len(biological_rows)} biological libraries, "
        f"and {len(cohort_rows)} target cohorts"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
