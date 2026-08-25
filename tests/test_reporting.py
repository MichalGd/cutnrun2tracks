#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReportingTests(unittest.TestCase):
    def test_lightweight_and_custom_reports_recover_from_retained_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            files = {
                "00_metadata/sample_manifest.tsv": "sample_key\tcohort_id\nTARGET.bioR1\tCOHORT_A\n",
                "00_metadata/cohort_manifest.tsv": "cohort_id\tn_samples\nCOHORT_A\t1\n",
                "06_qc/alignment_and_complexity/observation_counts.tsv": (
                    "sample_key\tlayout\tsignal_unit\tanalysis_observations\n"
                    "TARGET.bioR1\tPE\tfragment\t1234\n"
                ),
                "05_peaks/per_sample/peakcall_status.tsv": (
                    "sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\t"
                    "primary_peak_count\tcaller_warnings\treason\n"
                    "TARGET.bioR1\tCTRL.bioR1\tmacs3\tbroad\tEMPTY\t0\tmacs3:broad=EMPTY\t"
                    "primary_caller_produced_no_peaks\n"
                ),
                "05_peaks/consensus/consensus_status.tsv": (
                    "cohort_id\tstatus\ttotal_samples\tsuccessful_peak_samples\texcluded_samples\tregions\treason\n"
                    "COHORT_A\tSKIPPED\t1\t0\t1\t0\tinsufficient successful samples\n"
                ),
                "04_tracks/normalized_track_family_status.tsv": (
                    "cohort_id\tpolicy\tstatus\treason\tlog\n"
                    "COHORT_A\tanalysis\tSKIPPED\tconsensus unavailable\t.\n"
                ),
                "08_differential/stage_status.tsv": (
                    "status\tfailed_modules\tskipped_cohorts\nCOMPLETED_WITH_WARNINGS\t0\t1\n"
                ),
            }
            for relative, content in files.items():
                path = output / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content, encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/generate_report.py"), str(output)],
                text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report_dir = output / "10_reports"
            for name in ("pipeline_report.html", "run_summary.tsv", "warning_summary.tsv"):
                self.assertTrue((report_dir / name).is_file(), name)
                self.assertGreater((report_dir / name).stat().st_size, 0)
            warnings = (report_dir / "warning_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("peakcalling:TARGET.bioR1:EMPTY", warnings)
            self.assertIn("normalization:COHORT_A:analysis:SKIPPED", warnings)

            custom_dir = Path(temporary) / "custom"
            module = load_module(
                "prepare_multiqc_content",
                ROOT / "scripts/prepare_multiqc_content.py",
            )
            original_argv = sys.argv
            try:
                sys.argv = ["prepare_multiqc_content.py", str(output), str(custom_dir)]
                self.assertEqual(module.main(), 0)
            finally:
                sys.argv = original_argv
            manifest = (custom_dir / "custom_content_manifest.tsv").read_text(encoding="utf-8")
            self.assertIn("cutnrun2tracks_observations_mqc.tsv", manifest)
            self.assertIn("cutnrun2tracks_peakcalls_mqc.tsv", manifest)
            custom_peakcalls = (custom_dir / "cutnrun2tracks_peakcalls_mqc.tsv").read_text(encoding="utf-8")
            self.assertIn("# plot_type: table", custom_peakcalls)
            self.assertIn("TARGET.bioR1", custom_peakcalls)


if __name__ == "__main__":
    unittest.main()
