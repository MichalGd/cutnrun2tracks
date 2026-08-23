# Inputs and configuration

The samplesheet header must exactly match `config/samplesheet_template.csv`.
One row is a sequencing unit; identical `sample_id + replicate` rows with
different `tech_replicate` values are merged before trimming. Biological
replicates remain independent.

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

The optional metagene stage uses a separate tab-delimited gene-set manifest so
annotation releases and HPA-derived subsets can be replaced without modifying
workflow code. Filtering thresholds belong to the offline reference-build
command and are recorded in its reference manifest; plotting-window and
rendering settings belong to `config.conf`.
