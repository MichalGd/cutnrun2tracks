# Contributing

Contributions are welcome. Reproducibility, explicit scientific assumptions,
and separation of assay targets are the primary review criteria.

## Before proposing a change

- Open or reference an issue describing the problem, intended behavior, and
  affected assay/layout/profile.
- Keep CUT&RUN/CUT&Tag assumptions separate from ATAC-seq assumptions.
- Do not update shared references or an ATACseq2tracks installation from this
  repository.
- Do not commit FASTQ, BAM, bedGraph, bigWig, peak, or run-output directories.

## Scientific invariants

Changes must preserve:

- complete target-cohort isolation, including factor and antibody ID;
- exact, explicit target-to-control mapping;
- primary differential analysis from raw target fragment/read counts;
- clear separation of primary and control-aware sensitivity results;
- independence of CPM, DESeq2-consensus, robust-CPM, and spike-in track
  families;
- fragment-aware paired-end and read-aware single-end signal units;
- immutable reference provenance and checksummed run metadata.

Any intentional change to one of these invariants needs a documented scientific
rationale, a version bump, and a synthetic regression test.

## Validation

Run the complete repository checks in the locked Linux environment:

```bash
bash tests/check_bash_syntax.sh
bash tests/run_tests.sh
```

Workflow or scientific changes also require a small representative end-to-end
dataset for every affected combination of assay profile, PE/SE layout, peak
caller, control policy, and spike-in mode. Record tool versions, commands,
reference releases, checksums, and expected-versus-observed results.

Documentation-only changes should at least verify local Markdown links and
confirm that code, configuration defaults, environment files, and tests were
not modified unintentionally.

## Configuration and documentation

Every new configuration parameter must be:

1. added to `config/config.conf.template` with a safe default;
2. parsed and validated explicitly;
3. covered by an interface or scientific regression test;
4. described in the relevant page under `docs/`;
5. recorded in `CHANGELOG.md` with the associated version.

Update the root workflow diagram and stage/output documentation whenever a
stage boundary or data dependency changes. Release manifests must record script
and reference checksums.

## Reporting a problem

Include:

- the workflow version from `VERSION`;
- the assay profile, genome, PE/SE layout, and failing stage name;
- the exact launch command with sensitive paths redacted;
- the relevant stage log and error message;
- `00_metadata/preflight_status.tsv`, the resolved metadata manifests, and the
  relevant checkpoint JSON when safe to share;
- the output of `conda list` restricted to the tools involved in the failure.

Do not attach raw sequencing data, credentials, private sample identifiers, or
unredacted absolute paths to a public issue.
