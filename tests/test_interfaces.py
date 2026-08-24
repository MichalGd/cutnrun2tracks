#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class InterfaceTests(unittest.TestCase):
    def test_filtering_uses_indexed_marked_bam_for_region_selection(self) -> None:
        script = (ROOT / "scripts/mark_filter_batch.sh").read_text(encoding="utf-8")
        self.assertNotIn('"$tmp/flags.bam"', script)
        self.assertIn('"$marked" "${contigs[@]}"', script)
        self.assertIn("Reusing validated marked BAM", script)
        self.assertNotIn("trap 'rm -rf \"$tmp\"' RETURN", script)
        self.assertIn("trap 'rm -rf -- \"$tmp\"' EXIT", script)

    def test_config_template_is_complete_and_safe(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = directory / "config.conf"
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            resolved = (directory / "resolved.conf").read_text(encoding="utf-8")
            for key in (
                "GENERATE_CPM_TRACKS", "GENERATE_DESEQ2_CONSENSUS_TRACKS",
                "GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS",
                "GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS",
                "GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS",
            ):
                self.assertIn(f"{key}=true", resolved)
            self.assertIn("THREADS_FASTQC=10", resolved)
            self.assertIn("THREADS_TRIMGALORE=8", resolved)
            self.assertIn("PEAKCALL_FAILURE_POLICY=continue", resolved)

    def test_config_rejects_unknown_peakcall_failure_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            text = text.replace("PEAKCALL_FAILURE_POLICY=continue", "PEAKCALL_FAILURE_POLICY=ignore")
            config = directory / "config.conf"
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("PEAKCALL_FAILURE_POLICY must be one of", result.stderr)

    def test_config_rejects_nonpositive_preprocessing_threads(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            text = text.replace("THREADS_TRIMGALORE=8", "THREADS_TRIMGALORE=0")
            config = directory / "config.conf"
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("THREADS_TRIMGALORE must be a positive integer", result.stderr)

    def test_config_rejects_unknown_key(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results")) + "\nUNKNOWN_SETTING=yes\n"
            config = directory / "config.conf"; config.write_text(text, encoding="utf-8")
            result = subprocess.run([sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"), "--write-shell", str(directory / "resolved")],
                text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unknown configuration key", result.stderr)

    def test_samplesheet_resolves_control_and_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "metadata"
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_samplesheet.py"),
                str(ROOT / "config/examples/cutrun_pe.csv"), "--assay-profile", "cutrun",
                "--spikein-mode", "none", "--peak-callers", "seacr,macs3",
                "--primary-peak-caller", "auto", "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "sample_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            target = next(row for row in rows if row["is_control"] == "FALSE")
            self.assertEqual(target["control_key"], "WT_IgG.bioR1")
            self.assertEqual(target["primary_peak_caller"], "seacr")
            self.assertNotEqual(target["cohort_id"], ".")

    def test_seacr_is_rejected_for_single_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_samplesheet.py"),
                str(ROOT / "config/examples/cutrun_se.csv"), "--assay-profile", "cutrun",
                "--spikein-mode", "none", "--peak-callers", "seacr,macs3",
                "--primary-peak-caller", "auto", "--output-dir", str(Path(temporary) / "metadata"),
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("SEACR is PE-only", result.stderr)

    def test_different_antibodies_create_separate_cohorts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); source = ROOT / "config/examples/cutrun_pe.csv"
            with source.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            target = dict(rows[0]); second = dict(target)
            second.update(sample_id="WT_H3K27ac", factor="H3K27ac", antibody_id="AB_H3K27AC_01",
                          output_prefix="WT_H3K27ac", fastq_1="/data/WT_H3K27ac_R1.fastq.gz",
                          fastq_2="/data/WT_H3K27ac_R2.fastq.gz")
            sheet = directory / "samples.csv"
            with sheet.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
                writer.writeheader(); writer.writerows([target, second, rows[1]])
            out = directory / "metadata"
            rejected = subprocess.run([sys.executable, str(ROOT / "scripts/validate_samplesheet.py"), str(sheet),
                "--assay-profile", "cutrun", "--spikein-mode", "none", "--peak-callers", "seacr,macs3",
                "--primary-peak-caller", "auto", "--output-dir", str(directory / "rejected")],
                text=True, capture_output=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("assigned to multiple targets", rejected.stderr)
            result = subprocess.run([sys.executable, str(ROOT / "scripts/validate_samplesheet.py"), str(sheet),
                "--assay-profile", "cutrun", "--spikein-mode", "none", "--peak-callers", "seacr,macs3",
                "--primary-peak-caller", "auto", "--allow-shared-controls", "--output-dir", str(out)],
                text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "cohort_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                cohorts = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(cohorts), 2)
            self.assertNotEqual(cohorts[0]["cohort_id"], cohorts[1]["cohort_id"])


if __name__ == "__main__":
    unittest.main()
