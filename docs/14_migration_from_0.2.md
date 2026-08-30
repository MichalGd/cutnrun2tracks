# Migrating a 0.2.x project

[Documentation index](README.md) | [Inputs and configuration](01_inputs_and_configuration.md)

Version 0.3 introduced a stricter metadata contract, additional resource and QC
controls, optional epic2 support, expanded annotation, and a pinned shared
launcher. Do not reuse a 0.2.x configuration unchanged.

## Safe migration procedure

1. Copy `config/config.conf.template` from the installed 0.3.1 release into the
   project configuration directory.
2. Transfer only project-specific reference paths, output path, analysis
   choices, and resource limits.
3. Remove the obsolete run-wide `ASSAY_PROFILE` setting. Keep `assay_profile`
   on each samplesheet row.
4. Remove the old samplesheet `blacklist` column. Configure
   `BLACKLIST_<GENOME>` once in `config.conf`.
5. Keep row-level `genome`, `layout`, FASTQ paths, biological metadata, and
   exact control mappings in the samplesheet; do not duplicate them in config.
6. Review the new QC, annotation, logging, and resource-budget keys rather than
   relying on an earlier default.
7. Validate the metadata plan and inspect cohort membership:

   ```bash
   cutnrun2tracks --config /absolute/path/config.conf --plan
   column -t -s $'\t' /absolute/path/output/00_metadata/cohort_membership.tsv | less -S
   ```

8. Run complete preflight in the pinned runtime environment:

   ```bash
   cutnrun2tracks --config /absolute/path/config.conf --preflight-only
   ```

## Output and checkpoint compatibility

The 0.3 metadata contract changes run signatures and cohort definitions.
Existing 0.2.x checkpoints must not be adopted as if they were 0.3 outputs.
Use a new output directory for a new analysis. If recovery from retained 0.2
files is scientifically necessary, audit every upstream artifact and begin at
the earliest affected stage; copying a checkpoint JSON is not sufficient.

For historical implementation details, see [CHANGELOG.md](../CHANGELOG.md).
Current users should follow the 0.3.1 configuration template and current
method pages, not old release notes.
