# Documentation coverage compared with ATACseq2tracks

[Documentation index](README.md)

This audit compares the current topic coverage of `cutnrun2tracks/docs` with
the numbered pages in `ATACseq2tracks/docs`. It is a planning aid, not a claim
that CUT&RUN/CUT&Tag and ATAC-seq should use identical scientific methods.

## Topic mapping

| ATACseq2tracks page | Closest cutnrun2tracks coverage | Assessment |
|---|---|---|
| `01_overview.md` | root `README.md`, `02_methods.md` | covered |
| `02_quickstart.md` | `19_server_quick_start.md`, root `README.md` | covered |
| `03_installation.md` | `08_server_installation.md` | covered; CUT page is server-deployment focused |
| `04_inputs.md` | `01_inputs_and_configuration.md`, `config/README.md` | covered |
| `05_running.md` | root `README.md`, `04_outputs_and_recovery.md` | covered, but split |
| `06_pipeline_steps.md` | `07_pipeline_stages.md` | covered |
| `07_outputs.md` | `04_outputs_and_recovery.md` | covered |
| `08_diffbind.md` | `03_differential_enrichment.md` | covered within the complete differential-method page |
| `09_troubleshooting.md` | `16_troubleshooting.md` | covered with status-table and checkpoint-aware CUT diagnostics |
| `10_reference_files.md` | `17_reference_preparation_and_provenance.md` | covered; specialist details remain linked from blacklist, annotation, metagene, and installation pages |
| `11_blacklist_filtering.md` | `10_references_blacklist_and_filtering.md` | covered |
| `12_post_alignment_qc.md` | `11_quality_control.md` | covered |
| `13_differential_accessibility.md` | `03_differential_enrichment.md` | covered with CUT-specific models |
| `14_replicates_and_experimental_design.md` | `18_replicates_and_experimental_design.md` | covered with CUT-specific control and cohort rules |
| `15_deseq2atac_legacy_method_review.md` | no direct equivalent | ATAC-specific legacy review; no CUT page needed unless a CUT method is retired |
| `16_peak_annotation.md` | `13_genomic_annotation.md` | covered with CUT-specific schemas and reference provenance |

ATACseq2tracks also retains versioned update and hotfix pages. `cutnrun2tracks`
keeps release history in `CHANGELOG.md` and the focused migration page instead
of adding release-event documents to the current user manual.

## Pages specific to cutnrun2tracks

These CUT pages have no direct standalone ATACseq2tracks counterpart:

| cutnrun2tracks page | Subject |
|---|---|
| `05_limitations.md` | consolidated scientific and implementation limitations |
| `06_metagene.md` | TSS, TES, and scaled-gene-body aggregate signal |
| `06_test_matrix.md` | automated and real-data validation coverage |
| `09_peak_calling.md` | caller-specific CUT&RUN/CUT&Tag parameters and fault isolation |
| `12_tracks_and_normalization.md` | unified CPM, DESeq2, robust-CPM, and spike-in track definitions |
| `14_migration_from_0.2.md` | metadata and checkpoint migration guidance |
| `16_troubleshooting.md` | stage/status-aware diagnosis and safe recovery |
| `17_reference_preparation_and_provenance.md` | complete host, spike-in, annotation, and metagene reference contract |
| `18_replicates_and_experimental_design.md` | technical/biological replication, controls, cohorts, blocking, and confounding |
| `19_server_quick_start.md` | one-page launch and monitoring path for the shared deployment |

## Remaining page-level differences

The four previously identified user-facing gaps are now covered. Remaining
differences are primarily organizational or assay-specific:

- ATACseq2tracks has a DESeq2ATAC legacy-method review; no CUT analogue is
  needed unless a CUT method is retired and requires a migration audit.
- ATACseq2tracks retains individual historical update/hotfix pages;
  `cutnrun2tracks` deliberately keeps release events in `CHANGELOG.md` and
  current-use migration guidance in one page.
- CUT has standalone metagene, limitation, test-matrix, peak-calling, and
  normalization pages because these need CUT-specific interpretation.

Future pages should be added for a concrete user need rather than to make the
two repositories' filenames symmetrical.
