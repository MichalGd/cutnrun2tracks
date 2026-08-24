#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


class ScientificHelperTests(unittest.TestCase):
    def test_consensus_support_counts_biological_samples(self) -> None:
        module = load("build_consensus", ROOT / "scripts/build_consensus.py")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            a = directory / "a.bed"; b = directory / "b.bed"; c = directory / "c.bed"
            a.write_text("chr1\t10\t30\nchr1\t20\t40\n", encoding="utf-8")
            b.write_text("chr1\t20\t35\n", encoding="utf-8")
            c.write_text("chr1\t25\t50\n", encoding="utf-8")
            result = module.consensus([("A", a), ("B", b), ("C", c)], 2)
            self.assertEqual(result, [("chr1", 20, 40, 3)])

    def test_consensus_excludes_failed_peak_samples_and_continues(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output_root = directory / "05_peaks/consensus"
            sample_manifest = directory / "sample_manifest.tsv"
            cohort_manifest = directory / "cohort_manifest.tsv"
            sample_manifest.write_text(
                "sample_key\tis_control\tcohort_id\n"
                "A\tFALSE\tC1\nB\tFALSE\tC1\nC\tFALSE\tC1\n",
                encoding="utf-8",
            )
            cohort_manifest.write_text(
                "cohort_id\tprimary_peak_caller\tprimary_peak_class\n"
                "C1\tmacs3\tbroad\n",
                encoding="utf-8",
            )
            for sample, status, intervals, reason in (
                ("A", "SUCCESS", "chr1\t10\t30\n", "."),
                ("B", "SUCCESS", "chr1\t20\t40\n", "."),
                ("C", "ERROR", "", "primary_caller_error"),
            ):
                root = directory / "05_peaks/per_sample" / sample
                peak_dir = root / "macs3"
                peak_dir.mkdir(parents=True)
                (peak_dir / f"{sample}.macs3.broad.bed").write_text(intervals, encoding="utf-8")
                (root / "peakcall_metadata.tsv").write_text(
                    "sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\tprimary_peak_count\tcaller_warnings\treason\n"
                    f"{sample}\t.\tmacs3\tbroad\t{status}\t{int(bool(intervals))}\t.\t{reason}\n",
                    encoding="utf-8",
                )
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/build_consensus.py"),
                "--sample-manifest", str(sample_manifest),
                "--cohort-manifest", str(cohort_manifest),
                "--output-root", str(output_root),
                "--minimum-support", "2",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            consensus_files = list(output_root.rglob("*.consensus.bed"))
            self.assertEqual(len(consensus_files), 1)
            self.assertEqual(consensus_files[0].read_text(encoding="utf-8").split("\t")[:3],
                             ["chr1", "20", "30"])
            with (output_root / "consensus_status.tsv").open(encoding="utf-8", newline="") as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["status"], "SUCCESS")
            self.assertEqual(summary["successful_peak_samples"], "2")
            self.assertEqual(summary["excluded_samples"], "1")
            exclusion_files = list(output_root.rglob("excluded_peak_samples.tsv"))
            with exclusion_files[0].open(encoding="utf-8", newline="") as handle:
                exclusions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(exclusions, [{
                "sample_key": "C", "status": "ERROR", "reason": "primary_caller_error"
            }])

    def test_spikein_formula_and_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            counts = directory / "counts.tsv"; output = directory / "scales.tsv"
            counts.write_text(
                "sample_key\tlayout\tcohort_id\tspikein_stage\tspikein_lot\tspikein_to_host_ratio\thost_observations\tspike_observations\n"
                "A\tPE\tC1\tnuclei\tL1\t2\t90000\t10000\n", encoding="utf-8")
            result = subprocess.run([sys.executable, str(ROOT / "scripts/spikein_qc.py"), str(counts), str(output),
                "--scale-target", "1000000", "--fail-below", "1000", "--warn-below", "10000",
                "--warn-low-fraction", "0.001", "--warn-high-fraction", "0.20"], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertAlmostEqual(float(row["host_scale_factor"]), 200.0)
            self.assertAlmostEqual(float(row["spike_cpm_scale_factor"]), 100.0)
            self.assertEqual(row["status"], "PASS")

    def test_checkpoint_detects_changed_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); source = directory / "source"; output = directory / "output"
            source.write_text("a", encoding="utf-8"); output.write_text("result", encoding="utf-8")
            signature = subprocess.check_output([sys.executable, str(ROOT / "scripts/checkpoint.py"), "signature", str(source)], text=True).strip()
            checkpoint = directory / "checkpoint.json"
            subprocess.check_call([sys.executable, str(ROOT / "scripts/checkpoint.py"), "write", "--checkpoint", str(checkpoint),
                "--stage", "test", "--signature", signature, "--outputs", str(output)])
            valid = subprocess.run([sys.executable, str(ROOT / "scripts/checkpoint.py"), "check", "--checkpoint", str(checkpoint),
                "--stage", "test", "--signature", signature])
            self.assertEqual(valid.returncode, 0)
            output.write_text("changed", encoding="utf-8")
            invalid = subprocess.run([sys.executable, str(ROOT / "scripts/checkpoint.py"), "check", "--checkpoint", str(checkpoint),
                "--stage", "test", "--signature", signature])
            self.assertNotEqual(invalid.returncode, 0)

    def test_checkpoint_adopts_valid_prior_stage_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            output = directory / "output.txt"
            output.write_text("validated prior result", encoding="utf-8")
            checkpoint = directory / "checkpoint.json"
            old_signature = "old-signature"
            new_signature = "new-signature"
            subprocess.check_call([
                sys.executable, str(ROOT / "scripts/checkpoint.py"), "write",
                "--checkpoint", str(checkpoint), "--stage", "alignment",
                "--signature", old_signature, "--outputs", str(output),
            ])
            adopted = subprocess.run([
                sys.executable, str(ROOT / "scripts/checkpoint.py"), "adopt",
                "--checkpoint", str(checkpoint), "--stage", "alignment",
                "--signature", new_signature,
            ])
            self.assertEqual(adopted.returncode, 0)
            valid = subprocess.run([
                sys.executable, str(ROOT / "scripts/checkpoint.py"), "check",
                "--checkpoint", str(checkpoint), "--stage", "alignment",
                "--signature", new_signature,
            ])
            self.assertEqual(valid.returncode, 0)
            payload = json.loads(checkpoint.read_text(encoding="utf-8"))
            self.assertEqual(payload["signature"], new_signature)
            self.assertEqual(payload["signature_adoptions"][-1]["from"], old_signature)
            output.write_text("changed", encoding="utf-8")
            rejected = subprocess.run([
                sys.executable, str(ROOT / "scripts/checkpoint.py"), "adopt",
                "--checkpoint", str(checkpoint), "--stage", "alignment",
                "--signature", "another-signature",
            ])
            self.assertNotEqual(rejected.returncode, 0)


if __name__ == "__main__":
    unittest.main()
