#!/usr/bin/env python3
"""Validate spike observations and calculate explicit calibration factors."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("counts", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--scale-target", type=float, required=True)
    parser.add_argument("--fail-below", type=int, required=True)
    parser.add_argument("--warn-below", type=int, required=True)
    parser.add_argument("--warn-low-fraction", type=float, required=True)
    parser.add_argument("--warn-high-fraction", type=float, required=True)
    parser.add_argument("--allow-failed", action="store_true")
    args = parser.parse_args()
    with args.counts.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    provisional = []
    for row in rows:
        host, spike, ratio = int(row["host_observations"]), int(row["spike_observations"]), float(row["spikein_to_host_ratio"])
        if host <= 0 or spike <= 0 or not math.isfinite(ratio) or ratio <= 0:
            scale = math.nan
        else:
            scale = args.scale_target * ratio / spike
        provisional.append(scale)
    valid_scales = [value for value in provisional if math.isfinite(value) and value > 0]
    median_scale = statistics.median(valid_scales) if valid_scales else math.nan
    failed = False
    output_rows = []
    for row, scale in zip(rows, provisional):
        host, spike = int(row["host_observations"]), int(row["spike_observations"])
        fraction = spike / (host + spike) if host + spike > 0 else math.nan
        reasons, warnings = [], []
        if host <= 0:
            reasons.append("zero_host_observations")
        if spike < args.fail_below:
            reasons.append("spike_below_hard_threshold")
        elif spike < args.warn_below:
            warnings.append("spike_below_warning_threshold")
        if row["spikein_stage"] == "post_library":
            reasons.append("post_library_not_upstream_calibrator")
        if math.isfinite(fraction) and fraction < args.warn_low_fraction:
            warnings.append("low_spike_fraction")
        if math.isfinite(fraction) and fraction > args.warn_high_fraction:
            warnings.append("high_spike_fraction")
        if math.isfinite(scale) and math.isfinite(median_scale) and (scale > 10 * median_scale or scale < median_scale / 10):
            warnings.append("scale_more_than_tenfold_from_median")
        status = "FAILED" if reasons else ("WARN" if warnings else "PASS")
        failed |= bool(reasons)
        output_rows.append({
            **row, "spike_fraction": f"{fraction:.12g}", "host_scale_factor": f"{scale:.15g}",
            "spike_cpm_scale_factor": f"{(1e6/spike) if spike else math.nan:.15g}",
            "cohort_median_host_scale": f"{median_scale:.15g}", "status": status,
            "failure_reasons": ",".join(reasons) or ".", "warnings": ",".join(warnings) or ".",
        })
    fields = list(output_rows[0]) if output_rows else []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)
    return 1 if failed and not args.allow_failed else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print(f"SPIKE-IN ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
