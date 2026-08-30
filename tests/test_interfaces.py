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
    def test_final_report_runs_unified_multiqc_and_validates_outputs(self) -> None:
        report = (ROOT / "scripts/report_batch.sh").read_text(encoding="utf-8")
        unified = (ROOT / "scripts/generate_multiqc_report.sh").read_text(encoding="utf-8")
        recovery = (ROOT / "utilities/regenerate_reports.sh").read_text(encoding="utf-8")
        self.assertIn('generate_multiqc_report.sh', report)
        self.assertIn('cutnrun2tracks_multiqc_report.html', report)
        self.assertIn('--exclude deeptools', unified)
        self.assertIn('multiqc_custom_content_manifest.tsv', unified)
        self.assertIn('report_checksums.sha256', recovery)

    def test_filtering_uses_indexed_marked_bam_for_region_selection(self) -> None:
        script = (ROOT / "scripts/mark_filter_batch.sh").read_text(encoding="utf-8")
        self.assertNotIn('"$tmp/flags.bam"', script)
        self.assertIn('"$marked" "${contigs[@]}"', script)
        self.assertIn("Reusing validated marked BAM", script)
        self.assertNotIn("trap 'rm -rf \"$tmp\"' RETURN", script)
        self.assertIn("trap 'rm -rf -- \"$tmp\"' EXIT", script)

    def test_zero_consensus_counts_are_diagnostic_and_cohort_local(self) -> None:
        factors = (ROOT / "scripts/consensus_track_factors.R").read_text(encoding="utf-8")
        normalized = (ROOT / "scripts/normalized_tracks_batch.sh").read_text(encoding="utf-8")
        differential = (ROOT / "scripts/differential_batch.sh").read_text(encoding="utf-8")
        self.assertIn('"consensus_count_sums.tsv"', factors)
        self.assertIn("zero consensus counts for samples:", factors)
        self.assertIn("skip_or_fail_family", normalized)
        self.assertIn("normalized_track_family_status.tsv", normalized)
        self.assertIn("consensus normalization unavailable", differential)

    def test_annotation_uses_configured_genome_order(self) -> None:
        script = (ROOT / "scripts/annotate_browser.sh").read_text(encoding="utf-8")
        self.assertGreaterEqual(script.count('bedtools sort -faidx "$chrom_sizes"'), 2)
        self.assertIn('bedtools closest -a "$sorted_consensus" -b "$genes" -d -g "$chrom_sizes"', script)
        self.assertIn("summarize_peak_annotations.py", script)
        feature_script = (ROOT / "scripts/summarize_peak_annotations.py").read_text(encoding="utf-8")
        self.assertIn("axis.barh", feature_script)
        self.assertIn("peak_feature_counts.tsv", feature_script)
        self.assertIn("peak_feature_fractions.tsv", feature_script)

    def test_launcher_pins_environment_without_interactive_activation(self) -> None:
        launcher = (ROOT / "utilities/cutnrun2tracks_shared_launcher.sh").read_text(encoding="utf-8")
        self.assertIn("CUTNRUN2TRACKS_MAIN_ENV", launcher)
        self.assertIn('export PATH="$MAIN_ENV/bin:/usr/local/bin:/usr/bin:/bin"', launcher)
        self.assertNotIn("source /opt/miniconda", launcher)

    def test_resource_and_structured_logging_interfaces_are_present(self) -> None:
        preflight = (ROOT / "scripts/preflight.sh").read_text(encoding="utf-8")
        driver = (ROOT / "cutnrun2tracks.sh").read_text(encoding="utf-8")
        self.assertIn("resource_budget.tsv", preflight)
        self.assertIn("TOTAL_CPU_BUDGET", preflight)
        self.assertIn("workflow_events.tsv", driver)
        self.assertIn("stage_timing.tsv", driver)
        self.assertIn("cutnrun2tracks.console.log", driver)

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
            self.assertNotIn("ASSAY_PROFILE=", resolved)
            self.assertIn("RUN_FEATURE_ANNOTATION_SUMMARY=true", resolved)
            self.assertIn("TOTAL_CPU_BUDGET=auto", resolved)
            self.assertIn("RUN_FASTQC_PER_TECHNICAL_UNIT=true", resolved)
            self.assertIn("RUN_REPLICATE_CORRELATION=true", resolved)
            self.assertIn("EPIC2_COMMAND=epic2", resolved)

    def test_samplesheet_keeps_sample_metadata_without_reference_paths(self) -> None:
        header = (ROOT / "config/samplesheet_template.csv").read_text(encoding="utf-8").strip().split(",")
        self.assertIn("assay_profile", header)
        self.assertIn("genome", header)
        self.assertNotIn("blacklist", header)

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
                str(ROOT / "config/examples/cutrun_pe.csv"), "--blacklist-map", "hg38=/refs/hg38.blacklist.bed",
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
            self.assertTrue((out / "cohort_membership.tsv").is_file())

    def test_epic2_is_auto_primary_for_broad_targets_when_enabled(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = ROOT / "config/examples/cuttag_pe.csv"
            out = directory / "metadata"
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_samplesheet.py"), str(source),
                "--blacklist-map", "mm39=/refs/mm39.blacklist.bed",
                "--spikein-mode", "none", "--peak-callers", "macs3,epic2",
                "--primary-peak-caller", "auto", "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "sample_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                target = next(row for row in csv.DictReader(handle, delimiter="\t")
                              if row["is_control"] == "FALSE")
            self.assertEqual(target["primary_peak_caller"], "epic2")
            self.assertEqual(target["primary_peak_class"], "broad")

    def test_seacr_is_rejected_for_single_end(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_samplesheet.py"),
                str(ROOT / "config/examples/cutrun_se.csv"), "--blacklist-map", "hg38=/refs/hg38.blacklist.bed",
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
                "--blacklist-map", "hg38=/refs/hg38.blacklist.bed", "--spikein-mode", "none", "--peak-callers", "seacr,macs3",
                "--primary-peak-caller", "auto", "--output-dir", str(directory / "rejected")],
                text=True, capture_output=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("assigned to multiple targets", rejected.stderr)
            result = subprocess.run([sys.executable, str(ROOT / "scripts/validate_samplesheet.py"), str(sheet),
                "--blacklist-map", "hg38=/refs/hg38.blacklist.bed", "--spikein-mode", "none", "--peak-callers", "seacr,macs3",
                "--primary-peak-caller", "auto", "--allow-shared-controls", "--output-dir", str(out)],
                text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "cohort_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                cohorts = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(cohorts), 2)
            self.assertNotEqual(cohorts[0]["cohort_id"], cohorts[1]["cohort_id"])


if __name__ == "__main__":
    unittest.main()
