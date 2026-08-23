# Outputs and recovery

Results use numbered directories from `00_metadata` through `10_reports`.
Sample files use `<output_prefix>.bioR<replicate>`. Cohort IDs contain a readable
prefix and an eight-character SHA-256 suffix; the complete key is retained in
`cohort_manifest.tsv`.

Each stage checkpoint is JSON containing the complete run signature and hashes
of declared outputs. The signature covers the sanitized CSV, resolved config,
workflow version, scripts, reusable `common/` modules, and reference manifest.
A changed or missing output invalidates the checkpoint.

Metagene outputs are retained under `06_qc/metagene`. Each sample/gene-set/mode
directory contains the deepTools matrix, exported profile values, sorted BED,
PNG/PDF profile, PNG/PDF heatmap, and task metadata. The run-wide
`artifacts.tsv` is the stable reporting interface.

Use `--from-stage NAME` to force that and later stages. Use `--stop-after NAME`
for controlled partial runs. Cleanup begins only after the final report exists,
targets only explicit children of the resolved output root, and records every
deletion in `00_metadata/cleanup_manifest.tsv`.
