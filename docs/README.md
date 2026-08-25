# cutnrun2tracks documentation

This index groups the documentation by user task. Start with the root
[README](../README.md) for the workflow map, installation, quick start, stage
summary, and principal output tree.

## Prepare a run

| Topic | Document |
|---|---|
| Samplesheet rows, technical and biological replicates, matched controls, configuration safety, and plan/preflight validation | [Inputs and configuration](01_inputs_and_configuration.md) |
| Complete key list and example samplesheets | [`config/README.md`](../config/README.md) and [`config/config.conf.template`](../config/config.conf.template) |
| Exact stage order, enable/disable behavior, checkpoints, and declared outputs | [Pipeline stages](07_pipeline_stages.md) |
| Installation beside an existing ATACseq2tracks server environment and shared reference collection | [Server installation and shared-reference reuse](08_server_installation.md) |

## Understand the analysis

| Topic | Document |
|---|---|
| PE/SE signal units, filtering branches, CPM and DESeq2 track formulas, cohort isolation, and caller behavior | [Methods and normalization](02_methods.md) |
| The distinct roles of IgG/input/mock controls, primary raw-count models, DiffBind, and sensitivity analyses | [Controls and differential enrichment](03_differential_enrichment.md) |
| TSS/TES/gene-body aggregate plots, BED12 reference preparation, HPA subsets, manifests, and standalone reuse | [Metagene aggregate-signal module](06_metagene.md) |

## Operate and validate the workflow

| Topic | Document |
|---|---|
| Output organization, checkpoint recovery, partial reruns, and guarded cleanup | [Outputs and recovery](04_outputs_and_recovery.md) |
| Scientific limitations and real-data pilot decisions | [Limitations and pilot decisions](05_limitations.md) |
| Synthetic coverage and outstanding Linux/real-data fixtures | [Test matrix](06_test_matrix.md) |

## Project information

- [Release history](../CHANGELOG.md)
- [Contribution requirements](../CONTRIBUTING.md)
- [Metagene shared-module interface](../common/metagene/README.md)

Documentation describes version 0.2.5 unless a page explicitly says otherwise.
The executable behavior is defined by `cutnrun2tracks.sh`, the scripts under
`scripts/`, and the validated configuration template.
