# cutnrun2tracks

`cutnrun2tracks` is a Linux Bash workflow for paired-end or single-end
CUT&RUN and CUT&Tag data. It is an independent scientific refactor of the
ATACseq2tracks v4.2.0 engineering base; it does not modify or call the ATAC
installation.

Version 0.2.0 implements:

- explicit run-wide `cutrun` and `cuttag` profiles;
- strict, non-executable configuration parsing and CSV validation;
- exact replicate-matched IgG/input/mock control mapping;
- Bowtie2 host or competitive host+spike alignment;
- duplicate marking once and four inspectable MAPQ/duplicate BAM branches;
- pair-safe canonical/blacklist filtering;
- CPM, DESeq2-consensus, and three DESeq2 robust-CPM track families;
- MACS3 `BAMPE` and CUT-specific SE presets, plus PE-only SEACR;
- target/antibody/caller-specific consensus peaks;
- fragment-level PE counting, FRiP, fragment-length, fingerprint, and
  descriptive TSS QC;
- reusable TSS-, TES-, and scaled-gene-body aggregate-signal profiles and
  heatmaps from upstream-normalized bigWigs;
- DESeq2Enrichment, DiffBind target-only analysis, control-subtracted
  sensitivity analysis, and an optional target-control interaction model;
- gene/cCRE reference overlaps, UCSC/IGV assets, provenance, checksummed JSON
  checkpoints, HTML reporting, and guarded cleanup.

Different factors and antibody IDs are separate cohorts. The workflow refuses
to normalize different targets together.

## Quick start

Execution requires Linux and Bash 4.4 or newer.

```bash
conda env create -f environment.yml
conda activate cutnrun2tracks-0.2.0

cp config/config.conf.template config/config.conf
cp config/examples/cutrun_pe.csv samplesheet.csv
# Edit all input/reference/output paths.

bash cutnrun2tracks.sh --config config/config.conf --plan
bash cutnrun2tracks.sh --config config/config.conf --preflight-only
bash cutnrun2tracks.sh --config config/config.conf
```

SEACR 1.3 must be installed separately from the
[upstream release](https://github.com/FredHutch/SEACR) and its immutable path
set in `SEACR_COMMAND`. Record the downloaded commit and checksum in the local
environment record.

## Important defaults

- Targets retain marked duplicates by default; controls remove them.
- The primary BAM is MAPQ 30. MAPQ 0 appears only in sensitivity tracks.
- PE alignment uses `--no-mixed --no-discordant --dovetail -I 10 -X 1000`.
- SEACR is rejected for SE libraries. MACS3 SE defaults are provisional:
  CUT&RUN 0/150 and CUT&Tag -75/150 shift/extension.
- IgG controls affect background/peak calling, not standard DESeq2 size
  factors. Spike-in calibration is a separate operation.
- Primary differential testing uses raw target fragments. Control subtraction
  and target-control interaction outputs are explicitly labelled sensitivity
  analyses.
- TSS and periodicity metrics are descriptive; no ATAC pass/fail thresholds are
  applied.
- Metagene plots never renormalize bigWig values. They use a versioned BED12
  gene-set manifest and retain zero-signal genes by default.

## Documentation

- [Inputs and configuration](docs/01_inputs_and_configuration.md)
- [Methods and track formulas](docs/02_methods.md)
- [Controls and differential enrichment](docs/03_differential_enrichment.md)
- [Outputs and recovery](docs/04_outputs_and_recovery.md)
- [Limitations and pilot decisions](docs/05_limitations.md)
- [Metagene aggregate-signal module](docs/06_metagene.md)

## Release status

This is a development release. Synthetic interface/scientific tests are
included, but release against small real CUT&RUN/CUT&Tag datasets and a locked
Linux environment is still required before shared production deployment.
