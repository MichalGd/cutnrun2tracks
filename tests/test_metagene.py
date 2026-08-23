#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "common/metagene"


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


class MetageneTests(unittest.TestCase):
    def test_reference_builder_selects_transcript_and_filters(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = write(root / "genes.gtf", """chr1\ttest\tgene\t101\t500\t.\t+\t.\tgene_id \"ENSGA.1\"; gene_type \"protein_coding\"; gene_name \"A\";
chr1\ttest\ttranscript\t101\t400\t.\t+\t.\tgene_id \"ENSGA.1\"; transcript_id \"ENSTA1.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"A\"; tag \"MANE_Select\";
chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"ENSGA.1\"; transcript_id \"ENSTA1.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"A\";
chr1\ttest\texon\t301\t400\t.\t+\t.\tgene_id \"ENSGA.1\"; transcript_id \"ENSTA1.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"A\";
chr1\ttest\ttranscript\t101\t500\t.\t+\t.\tgene_id \"ENSGA.1\"; transcript_id \"ENSTA2.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"A\";
chr1\ttest\texon\t101\t500\t.\t+\t.\tgene_id \"ENSGA.1\"; transcript_id \"ENSTA2.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"A\";
chr1\ttest\tgene\t601\t800\t.\t-\t.\tgene_id \"ENSGB.2\"; gene_type \"protein_coding\"; gene_name \"B\";
chr1\ttest\ttranscript\t601\t800\t.\t-\t.\tgene_id \"ENSGB.2\"; transcript_id \"ENSTB.2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"B\";
chr1\ttest\texon\t601\t800\t.\t-\t.\tgene_id \"ENSGB.2\"; transcript_id \"ENSTB.2\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"B\";
chr1\ttest\tgene\t901\t1000\t.\t+\t.\tgene_id \"ENSGC.1\"; gene_type \"protein_coding\"; gene_name \"C\";
chr1\ttest\ttranscript\t901\t1000\t.\t+\t.\tgene_id \"ENSGC.1\"; transcript_id \"ENSTC.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"C\";
chr1\ttest\texon\t901\t1000\t.\t+\t.\tgene_id \"ENSGC.1\"; transcript_id \"ENSTC.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"C\";
chrUn\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"ENSGD.1\"; gene_type \"protein_coding\"; gene_name \"D\";
chrUn\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id \"ENSGD.1\"; transcript_id \"ENSTD.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"D\";
chrUn\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"ENSGD.1\"; transcript_id \"ENSTD.1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\"; gene_name \"D\";
""")
            chrom_sizes = write(root / "chrom.sizes", "chr1\t2000\nchrUn\t500\n")
            contigs = write(root / "canonical.txt", "chr1\n")
            blacklist = write(root / "blacklist.bed", "chr1\t950\t960\n")
            out = root / "reference"
            result = subprocess.run([
                sys.executable, str(MODULE / "prepare_metagene_reference.py"),
                "--gtf", str(gtf), "--assembly", "GRCh38", "--annotation-source", "gencode",
                "--annotation-release", "test", "--chrom-sizes", str(chrom_sizes),
                "--canonical-contigs", str(contigs), "--blacklist", str(blacklist),
                "--min-gene-span", "0", "--min-spliced-length", "0", "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            bed_rows = [line.split("\t") for line in (out / "protein_coding.filtered.bed12").read_text().splitlines()]
            self.assertEqual([row[3] for row in bed_rows], ["ENSGA", "ENSGB"])
            self.assertEqual(bed_rows[0][1:3], ["100", "400"])
            self.assertEqual(bed_rows[0][9:12], ["2", "100,100,", "0,200,"])
            with (out / "gene_model_metadata.tsv").open(encoding="utf-8", newline="") as handle:
                metadata = {row["stable_gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(metadata["ENSGA"]["stable_transcript_id"], "ENSTA1")
            self.assertEqual(metadata["ENSGC"]["exclusion_reasons"], "blacklist_overlap")
            self.assertEqual(metadata["ENSGD"]["exclusion_reasons"], "noncanonical_contig")

    def test_hpa_and_mouse_ortholog_subset(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            human = write(root / "human.bed12", "chr1\t0\t100\tENSGA\t0\t+\t0\t100\t0\t1\t100,\t0,\nchr1\t200\t300\tENSGB\t0\t+\t200\t300\t0\t1\t100,\t0,\n")
            mouse = write(root / "mouse.bed12", "chr1\t0\t100\tENSMUSG1\t0\t+\t0\t100\t0\t1\t100,\t0,\n")
            hpa = write(root / "hpa.tsv", "Gene\tRNA tissue specificity\tRNA tissue distribution\nENSGA.9\tLow tissue specificity\tDetected in all\nENSGB\tLow tissue specificity\tDetected in many\n")
            orthologs = write(root / "orthologs.tsv", "Gene stable ID\tMouse gene stable ID\tMouse homology type\tMouse orthology confidence [0 low, 1 high]\nENSGA\tENSMUSG1\tortholog_one2one\t1\n")
            out = root / "subset"
            result = subprocess.run([
                sys.executable, str(MODULE / "build_hpa_reference_subset.py"),
                "--hpa-table", str(hpa), "--hpa-release", "test",
                "--human-gene-model", str(human), "--ortholog-table", str(orthologs),
                "--mouse-gene-model", str(mouse), "--ensembl-release", "test", "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("ENSGA", (out / "broadly_expressed.GRCh38.bed12").read_text())
            self.assertNotIn("ENSGB", (out / "broadly_expressed.GRCh38.bed12").read_text())
            self.assertIn("ENSMUSG1", (out / "broadly_expressed.GRCm39.bed12").read_text())

    def test_manifest_validator_resolves_all_modes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = write(root / "genes.bed12", "chr1\t0\t100\tENSGA\t0\t+\t0\t100\t0\t1\t100,\t0,\n")
            bigwig = write(root / "sample.bw", "placeholder")
            blacklist = write(root / "blacklist.bed", "chr1\t500\t600\n")
            chrom_sizes = write(root / "chrom.sizes", "chr1\t1000\n")
            tracks = write(root / "tracks.tsv", f"sample_id\tassay\tgenome\tbigwig\tnormalization\tnormalization_detail\tblacklist\tchrom_sizes\nS1\tcutrun\tGRCh38\t{bigwig}\tCPM\ttest\t{blacklist}\t{chrom_sizes}\n")
            genes = write(root / "gene_sets.tsv", f"gene_set_id\tgenome\tbed12\tlabel\nprotein_coding\tGRCh38\t{bed}\tProtein coding\n")
            tasks = root / "tasks.tsv"
            result = subprocess.run([
                sys.executable, str(MODULE / "validate_metagene_inputs.py"),
                "--track-manifest", str(tracks), "--gene-set-manifest", str(genes),
                "--gene-sets", "protein_coding", "--modes", "tss,tes,gene_body",
                "--output", str(tasks), "--skip-bigwig-validation",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with tasks.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["mode"] for row in rows}, {"tss", "tes", "gene_body"})
            self.assertEqual(len(rows), 3)

    def test_cut_adapter_selects_cpm_track(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "run"; (output / "04_tracks/cpm").mkdir(parents=True)
            bigwig = write(output / "04_tracks/cpm/S1.bioR1.CPM.bw", "placeholder")
            blacklist = write(root / "blacklist.bed", "chr1\t1\t2\n")
            chrom_sizes = write(root / "chrom.sizes", "chr1\t1000\n")
            samples = write(root / "samples.tsv", f"sample_key\tis_control\tgenome\tassay_profile\tblacklist\tcohort_id\nS1.bioR1\tFALSE\tGRCh38\tcutrun\t{blacklist}\tC1\n")
            config = write(root / "resolved.conf", f"SPIKEIN_MODE=none\nCHROM_SIZES_GRCH38={shlex_quote(str(chrom_sizes))}\n")
            manifest = root / "tracks.tsv"
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/build_metagene_track_manifest.py"),
                "--sample-manifest", str(samples), "--resolved-config", str(config),
                "--output-dir", str(output), "--track-family", "auto", "--output", str(manifest),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with manifest.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["normalization"], "CPM")
            self.assertEqual(Path(row["bigwig"]), bigwig.resolve())


def shlex_quote(value: str) -> str:
    import shlex
    return shlex.quote(value)


if __name__ == "__main__":
    unittest.main()
