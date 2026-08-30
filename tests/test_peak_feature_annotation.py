#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class PeakFeatureAnnotationTests(unittest.TestCase):
    def test_exclusive_counts_fractions_and_zero_categories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            metadata = output / "00_metadata"
            peak_root = output / "05_peaks/per_sample/S1.bioR1"
            annotation = output / "07_annotation/feature_summary"
            metadata.mkdir(parents=True)
            peak_root.mkdir(parents=True)
            peak_file = peak_root / "macs3/S1.bioR1.macs3.narrow.bed"
            peak_file.parent.mkdir()
            peak_file.write_text(
                "chr1\t950\t990\tpromoter\n"
                "chr1\t1100\t1150\texon\n"
                "chr1\t1500\t1600\tintron\n"
                "chr1\t2000\t2050\tgene_end\n"
                "chr1\t3000\t3100\tenhancer\n"
                "chr1\t4000\t4100\tintergenic\n"
                "chrUn\t0\t10\tunclassified\n",
                encoding="utf-8",
            )
            (peak_root / "caller_status.tsv").write_text(
                "sample_key\tcaller\tpeak_class\tstatus\tpeak_count\tpeak_file\tlog\treason\n"
                f"S1.bioR1\tmacs3\tnarrow\tSUCCESS\t7\t{peak_file}\t.\t.\n",
                encoding="utf-8",
            )
            (metadata / "sample_manifest.tsv").write_text(
                "sample_key\treplicate\tgenome\tfactor\tcondition\tis_control\tcohort_id\n"
                "S1.bioR1\t1\ttest\tH3K27ac\tWT\tFALSE\tC1\n",
                encoding="utf-8",
            )
            (metadata / "cohort_manifest.tsv").write_text(
                "cohort_id\tgenome\tfactor\tconditions\tprimary_peak_caller\tprimary_peak_class\n"
                "C1\ttest\tH3K27ac\tWT\tmacs3\tnarrow\n",
                encoding="utf-8",
            )
            gtf = root / "genes.gtf"
            gtf.write_text(
                'chr1\ttest\tgene\t1001\t2000\t.\t+\t.\tgene_id "G1"; gene_name "Gene1";\n'
                'chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "G1"; gene_name "Gene1";\n'
                'chr1\ttest\texon\t1801\t2000\t.\t+\t.\tgene_id "G1"; gene_name "Gene1";\n',
                encoding="utf-8",
            )
            chrom = root / "chrom.sizes"
            chrom.write_text("chr1\t5000\n", encoding="utf-8")
            ccre = root / "ccre.bed"
            ccre.write_text("chr1\t3000\t3100\tEH1\tdELS\n", encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/summarize_peak_annotations.py"),
                "--output-dir", str(output), "--annotation-dir", str(annotation),
                "--gtf", f"test={gtf}", "--chrom-sizes", f"test={chrom}",
                "--ccre", f"test={ccre}", "--promoter-upstream", "100",
                "--promoter-downstream", "50", "--gene-end-window", "100",
                "--plot-formats", "png", "--skip-plots",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (annotation / "peak_feature_summary.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            sample_rows = [row for row in rows if row["entity_id"] == "S1.bioR1"]
            self.assertEqual(len(sample_rows), 8)
            counts = {row["category"]: int(row["count"]) for row in sample_rows}
            for category in ("promoter", "enhancer", "exon", "intron", "gene_end", "intergenic", "unclassified"):
                self.assertEqual(counts[category], 1)
            self.assertEqual(counts["other_regulatory"], 0)
            self.assertAlmostEqual(sum(float(row["fraction"]) for row in sample_rows), 1.0)
            with gzip.open(annotation / "peak_feature_assignments.tsv.gz", "rt", encoding="utf-8") as handle:
                assignments = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(assignments), 7)


if __name__ == "__main__":
    unittest.main()
