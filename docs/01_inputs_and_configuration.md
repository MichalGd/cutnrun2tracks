# Inputs and configuration

The samplesheet header must exactly match `config/samplesheet_template.csv`.
One row is a sequencing unit; identical `sample_id + replicate` rows with
different `tech_replicate` values are merged before trimming. Biological
replicates remain independent.

The ownership boundary is deliberate:

- samplesheet: FASTQs, genome, PE/SE layout, CUT&RUN/CUT&Tag profile, target,
  antibody, biological condition, replicate identities, control mapping,
  duplicate policy, batch/donor, and optional spike metadata;
- `config.conf`: output directory, genome reference paths, caller/QC policy,
  resource limits, optional modules, logging, and retention.

Do not add blacklist/reference paths to the CSV. The validator resolves the
row's `genome` against `BLACKLIST_<GENOME>` and other assembly-specific config
keys. Mixed assay profiles are allowed because assay is sample metadata; mixed
genomes or layouts still require their explicit expert overrides.

Targets reference a control `sample_id` through `control_id`. The validator
first requires the same replicate. A condition=`shared` control is eligible
only when `ALLOW_SHARED_CONTROLS=true` and all other biological context matches.
There is no replicate-1 fallback.

Configuration is plain `KEY=VALUE`, not a sourced user shell script. Unknown
keys, duplicate keys, shell expansions, invalid booleans, and incomplete spike
metadata fail before tools run. The resolved configuration is written to
`00_metadata/resolved_config.conf` and `resolved_config.tsv`.

Run `--plan` for metadata/cohort validation without checking FASTQ/reference
existence or installed tools. Run `--preflight-only` for the full environment
and reference audit.

Both validation modes emit `sample_manifest.tsv`, `cohort_manifest.tsv`, and
`cohort_membership.tsv`. The membership table states which biological samples
and controls feed each cohort and lists every shared-control reuse, making the
otherwise automatic cohort derivation directly auditable.

The installed launcher needs only the config path:

```bash
cutnrun2tracks --config /absolute/path/to/config.conf
```

It pins the release and Conda environment internally; interactive
`conda activate` and manual PATH exports are not required.

The optional metagene stage uses a separate tab-delimited gene-set manifest so
annotation releases and HPA-derived subsets can be replaced without modifying
workflow code. Filtering thresholds belong to the offline reference-build
command and are recorded in its reference manifest; plotting-window and
rendering settings belong to `config.conf`.
