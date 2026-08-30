# Server installation and shared-reference reuse

This runbook installs `cutnrun2tracks` beside an already validated
ATACseq2tracks installation. It does not modify the production ATAC workflow,
its environment, or its shared genomic references.

## Recommended strategy

Reuse the ATAC environment as a tested package baseline, but do **not** extend
the live environment in place.

1. Clone `/opt/miniconda/envs/ATACseq2tracks` to a new CUT-specific prefix.
2. Add only the missing CUT dependencies to the clone with installed packages
   frozen where possible.
3. Install SEACR in an immutable, checksummed tool directory.
4. Reuse the existing host FASTA, Bowtie2 indexes, chromosome sizes, GTFs,
   blacklists, cCREs, and validated composite indexes through absolute paths.
5. Create only small CUT-specific derived lists when the existing resources do
   not have the coordinate names required by this workflow.
6. Validate a candidate release with synthetic tests and small real CUT data
   before creating the shared `current` symlink.

The final installation is intended for every server user. Workflow code,
Conda packages, SEACR, and shared references must therefore be readable and
traversable by ordinary users, while remaining writable only by the deployment
administrator. Individual users write configurations, logs, and results to
their own project directories, never into the shared installation.

Do not run `conda install` against
`/opt/miniconda/envs/ATACseq2tracks`. A solver transaction there could upgrade
R, Bioconductor, Python, or command-line tools and invalidate the already tested
ATACseq2tracks deployment.

## All-user access model

The instructions use four shared, read-only runtime locations:

```text
/opt/bioinformatics/workflows/cutnrun2tracks/
/opt/miniconda/envs/cutnrun2tracks-0.3.0/
/opt/bioinformatics/tools/SEACR/releases/1.3/
/opt/bioinformatics/references/
```

For ordinary server users:

- every parent directory must have directory traversal permission;
- workflow, environment, SEACR, and reference files must be readable;
- executable programs and scripts must be executable;
- no ordinary user should be able to modify releases, environments, tools,
  reference files, manifests, or the `current` symlink;
- each user must have a separate writable project configuration and output
  directory.

Use an installation umask of `022` so newly created shared directories are
normally `0755` and files are normally `0644` or `0755` when executable. A
site-managed read-only group/ACL is also acceptable, but it must contain every
intended workflow user. The validation in Step 13 must be run as a real
non-administrator account, not only as root or the deployment owner.

## Why cloning is preferable for the first server pilot

The installed ATAC environment already provides nearly all CUT dependencies:
Bowtie2, samtools, bedtools, Picard, Trim Galore, FastQC, MultiQC,
MACS3, deepTools, R, DESeq2, DiffBind, GenomicRanges, GenomicAlignments,
Rsamtools, rtracklayer, BiocParallel, ggplot2, and matplotlib.

The CUT specification adds or makes explicit:

| Component | Role | Expected action after cloning |
|---|---|---|
| Bash 5.1+ | workflow runtime | Verify; install the pinned Conda Bash if needed |
| `preseq` | optional/default complexity extrapolation | Install or set `RUN_PRESEQ=false` |
| `pyBigWig` | metagene bigWig validation | Install if `RUN_METAGENE=true` |
| `ucsc-bedClip` | environment-specification parity | Install; not a current preflight command |
| `cutadapt` | Trim Galore backend | Install the pinned standalone command |
| `pandas` | metagene reference and summary utilities | Install into the CUT clone |
| SEACR 1.3 | CUT peak calling | Install separately and pin the exact revision |
| epic2 sidecar | optional broad-domain peak calling | Create from `environment.epic2.yml`; expose only a version-pinned launcher |

The repository's `environment.lock.yml` is a platform-neutral version
constraint file, not a resolved `conda-lock` artifact with exact Linux package
URLs/build strings. It also proposes newer R/Bioconductor versions than the
reported validated ATAC environment. Do not treat it as the production lock
until that combination has been solved and tested on `linux-64`.

## 1. Define installation locations

Run the following as an account allowed to write the selected workflow,
environment, tool, and shared-reference directories. Adapt ownership/group
handling to the server's local policy.

```bash
set -euo pipefail
umask 022

CUT_VERSION=0.3.1
CUT_COMMIT=REPLACE_WITH_RELEASE_COMMIT
CUT_STAGE="/home/micgdu/Analysis/workflows/cutnrun2tracks/install_sources/${CUT_VERSION}"
CUT_WORKFLOW_ROOT=/opt/bioinformatics/workflows/cutnrun2tracks
CUT_RELEASE_DIR="${CUT_WORKFLOW_ROOT}/releases/${CUT_VERSION}"
CUT_DEPLOYMENT_DIR="${CUT_WORKFLOW_ROOT}/deployment/${CUT_VERSION}"
CUT_ENV=/opt/miniconda/envs/cutnrun2tracks-0.3.0
ATAC_ENV=/opt/miniconda/envs/ATACseq2tracks
ATAC_CURRENT=/opt/bioinformatics/workflows/ATACseq2tracks/current
SEACR_RELEASE=/opt/bioinformatics/tools/SEACR/releases/1.3
CUT_DERIVED_REFS=/opt/bioinformatics/references/cutnrun2tracks/0.2.0
```

Do not continue if any intended new release or environment destination already
contains an unrelated installation:

```bash
test ! -e "$CUT_RELEASE_DIR" && test ! -L "$CUT_RELEASE_DIR"
test ! -e "$CUT_ENV" && test ! -L "$CUT_ENV"
test -d "$ATAC_ENV"
test -e "$ATAC_CURRENT"
```

## 2. Audit the working ATAC installation and shared references

These commands are read-only:

```bash
readlink -f "$ATAC_CURRENT"
cat "$ATAC_CURRENT/VERSION"

source /opt/miniconda/etc/profile.d/conda.sh
conda activate "$ATAC_ENV"

for command_name in \
    bash python3 Rscript java bowtie2 bowtie2-inspect samtools bedtools \
    picard trim_galore cutadapt fastqc multiqc macs3 bamCoverage \
    computeMatrix plotProfile plotHeatmap plotFingerprint bedGraphToBigWig
do
    command -v "$command_name" || printf 'MISSING\t%s\n' "$command_name"
done

python3 --version
Rscript --version
bowtie2 --version | head -n 1
samtools --version | head -n 1
bedtools --version
macs3 --version
bamCoverage --version
```

Verify the R namespaces used by the CUT differential modules:

```bash
Rscript -e '
packages <- c(
  "DESeq2", "DiffBind", "BiocParallel", "GenomicAlignments",
  "GenomicRanges", "IRanges", "Rsamtools", "rtracklayer", "ggplot2"
)
for (package in packages) {
  cat(package, "\t", requireNamespace(package, quietly=TRUE), "\n", sep="")
}'
```

Verify the known shared host resources:

```bash
shared_files=(
  /opt/bioinformatics/references/hg38/hg38.fa
  /opt/bioinformatics/references/hg38/hg38.chrom.sizes
  /opt/bioinformatics/references/mm39/mm39.fa
  /opt/bioinformatics/references/mm39/mm39.chrom.sizes
  /opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf
  /opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.blacklist.bed
  /opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz
  /opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf
  /opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.blacklist.bed
  /opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed
  /opt/bioinformatics/ATACseq2tracks_shared/config/reference_manifest.tsv
)

for shared_file in "${shared_files[@]}"; do
    if [[ -s "$shared_file" && -r "$shared_file" ]]; then
        stat -c 'OK  %A  %U:%G  %s bytes  %n' "$shared_file"
    else
        printf 'FAIL\t%s\n' "$shared_file"
    fi
done

for index_prefix in \
    /opt/bioinformatics/references/hg38/bowtie2/hg38 \
    /opt/bioinformatics/references/mm39/bowtie2/mm39
do
    bowtie2-inspect -n "$index_prefix" >/dev/null
    printf 'OK\t%s\n' "$index_prefix"
done
```

Stop and resolve every `FAIL` before deployment. Do not copy large references
into the workflow release directory.

## 3. Deploy the immutable offline release

Set `CUT_COMMIT` to the exact commit referenced by the signed/tagged
`v${CUT_VERSION}` release. Transfer the archive, its portable SHA-256 sidecar,
and the commit record from the administrator workstation into
`$CUT_STAGE`. The server does not need GitHub access. Keep all CUT-specific
staging material under `/home/micgdu/Analysis/workflows/`; only the shared,
root-owned runtime installation belongs under `/opt`.

Validate the staged release before using administrator privileges:

```bash
cd "$CUT_STAGE"
sha256sum -c "cutnrun2tracks-${CUT_VERSION}.tar.gz.sha256"

test "$(tr -d '\r\n' < "cutnrun2tracks-${CUT_VERSION}.commit.txt")" = \
    "$CUT_COMMIT"
test "$(tar -xOzf "cutnrun2tracks-${CUT_VERSION}.tar.gz" \
    "cutnrun2tracks-${CUT_VERSION}/VERSION" | tr -d '\r\n')" = \
    "$CUT_VERSION"
```

Install the candidate release without creating `current`:

```bash
sudo install -d -o root -g root -m 0755 \
    "$CUT_WORKFLOW_ROOT" \
    "$CUT_WORKFLOW_ROOT/releases" \
    "$CUT_WORKFLOW_ROOT/deployment" \
    "$CUT_RELEASE_DIR" \
    "$CUT_DEPLOYMENT_DIR"

sudo tar --extract --gzip \
    --file "$CUT_STAGE/cutnrun2tracks-${CUT_VERSION}.tar.gz" \
    --directory "$CUT_RELEASE_DIR" \
    --strip-components=1 --no-same-owner

sudo install -o root -g root -m 0644 \
    "$CUT_STAGE/cutnrun2tracks-${CUT_VERSION}.commit.txt" \
    "$CUT_DEPLOYMENT_DIR/workflow_commit.txt"
sudo install -o root -g root -m 0644 \
    "$CUT_STAGE/cutnrun2tracks-${CUT_VERSION}.tar.gz.sha256" \
    "$CUT_DEPLOYMENT_DIR/workflow_archive.sha256"

sudo chown -R root:root "$CUT_RELEASE_DIR" "$CUT_DEPLOYMENT_DIR"
sudo chmod -R a+rX "$CUT_RELEASE_DIR" "$CUT_DEPLOYMENT_DIR"
sudo chmod -R go-w "$CUT_RELEASE_DIR" "$CUT_DEPLOYMENT_DIR"

test "$(tr -d '\r\n' < "$CUT_RELEASE_DIR/VERSION")" = "$CUT_VERSION"
test "$(tr -d '\r\n' < "$CUT_DEPLOYMENT_DIR/workflow_commit.txt")" = \
    "$CUT_COMMIT"
test -r "$CUT_RELEASE_DIR/cutnrun2tracks.sh"
```

The entry point is intentionally invoked through `bash`; repository mode
`100644` is therefore valid as long as the file is readable.

## 4. Clone and extend the ATAC environment safely

Clone the validated environment to the new prefix:

```bash
sudo -H nice -n 10 ionice -c 2 -n 7 \
    /opt/miniconda/bin/conda create --yes --copy \
    --prefix "$CUT_ENV" --clone "$ATAC_ENV"
```

`--copy` avoids hard-linking CUT environment files to the shared Conda package
cache or the validated ATAC installation. It uses more disk than a default
clone, but makes subsequent CUT permission and package management independent.

Record the clone before adding packages:

```bash
/opt/miniconda/bin/conda list --prefix "$CUT_ENV" --explicit | \
    sudo tee \
    "$CUT_DEPLOYMENT_DIR/conda-before-cut-additions.explicit.txt" \
    >/dev/null
```

Check which additions are actually missing:

```bash
for command_name in bash cutadapt preseq bedClip; do
    candidate="$CUT_ENV/bin/$command_name"
    if [[ -x "$candidate" ]]; then
        printf 'ENV_OK\t%s\n' "$candidate"
    else
        printf 'ENV_MISSING\t%s\n' "$command_name"
    fi
done

"$CUT_ENV/bin/python3" -c '
import importlib.util
for module in ("pyBigWig", "pandas", "matplotlib"):
    print(f"{module}\t{importlib.util.find_spec(module) is not None}")
'
```

The audited ATAC clone lacks Conda Bash, standalone cutadapt, preseq, bedClip,
and pandas. First perform a dry run. `--freeze-installed` minimizes changes to
the tested clone; abort if the solver proposes major R/Bioconductor/Python,
deepTools, pyBigWig, or MACS3 changes.

```bash
sudo -H /opt/miniconda/bin/conda install --dry-run --copy \
    --prefix "$CUT_ENV" --freeze-installed --override-channels \
    --channel conda-forge --channel bioconda \
    bash=5.2 cutadapt=5.1 preseq=3.2.0 ucsc-bedclip pandas
```

After reviewing and accepting that transaction, repeat it with `--yes` in place
of `--dry-run`. Do not reinstall deepTools or pyBigWig merely because they are
not listed among the additions; the clone already contains the validated
versions.

Record the completed environment:

```bash
/opt/miniconda/bin/conda list --prefix "$CUT_ENV" --explicit | \
    sudo tee "$CUT_DEPLOYMENT_DIR/conda-linux-64.explicit.txt" >/dev/null
/opt/miniconda/bin/conda env export --prefix "$CUT_ENV" | \
    sudo tee "$CUT_DEPLOYMENT_DIR/conda-environment.yml" >/dev/null
sha256sum \
    "$CUT_DEPLOYMENT_DIR/conda-before-cut-additions.explicit.txt" \
    "$CUT_DEPLOYMENT_DIR/conda-linux-64.explicit.txt" \
    "$CUT_DEPLOYMENT_DIR/conda-environment.yml" | \
    sudo tee "$CUT_DEPLOYMENT_DIR/conda-provenance.sha256" >/dev/null
```

The explicit package file is the reliable same-platform recreation artifact.
Once the pilot passes, it can be converted into the maintained Linux lock for
future releases.

### Clean-environment alternative

A clean environment can instead be solved from `environment.yml`:

```bash
mamba env create --prefix "$CUT_ENV" \
    --file "$CUT_RELEASE_DIR/environment.yml"
```

This is cleaner conceptually but may select versions different from the tested
ATAC environment. Use it after an end-to-end validation, or when cloning is not
permitted. Do not create both alternatives at the same prefix.

## 5. Install pinned SEACR 1.3

SEACR is not included in the Conda specification. The official
[SEACR repository](https://github.com/FredHutch/SEACR) provides tag `v1.3` at
commit `5a0efe59f06fb17cf9d34d415bb0c1a1f7a77a3c`. The shell wrapper requires
`SEACR_1.3.R` in the same directory and uses R and bedtools from `PATH`.

```bash
test ! -e "$SEACR_RELEASE" && test ! -L "$SEACR_RELEASE"
mkdir -p "$(dirname "$SEACR_RELEASE")"

git clone --branch v1.3 --depth 1 \
    https://github.com/FredHutch/SEACR.git "$SEACR_RELEASE"

test "$(git -C "$SEACR_RELEASE" rev-parse HEAD)" = \
    5a0efe59f06fb17cf9d34d415bb0c1a1f7a77a3c
test -s "$SEACR_RELEASE/SEACR_1.3.sh"
test -s "$SEACR_RELEASE/SEACR_1.3.R"
chmod 0755 "$SEACR_RELEASE/SEACR_1.3.sh"

sha256sum \
    "$SEACR_RELEASE/SEACR_1.3.sh" \
    "$SEACR_RELEASE/SEACR_1.3.R" \
    > "$CUT_DEPLOYMENT_DIR/SEACR_1.3.sha256"
```

Use the immutable absolute wrapper path in every CUT configuration:

```text
SEACR_COMMAND=/opt/bioinformatics/tools/SEACR/releases/1.3/SEACR_1.3.sh
```

Do not point production configurations at a moving SEACR branch or a user-home
copy.

## 6. Reuse host references and create canonical-contig lists

The host FASTA, Bowtie2 indexes, chromosome sizes, GTFs, blacklists, and cCRE
BEDs are assay-independent and should be reused unchanged. The CUT workflow
also requires a one-column canonical-contig file. First search for an existing
validated file:

```bash
find /opt/bioinformatics/references \
     /opt/bioinformatics/ATACseq2tracks_shared/references \
     -type f -iname '*canonical*contig*' -print
```

If suitable files do not exist, create small derived lists without modifying
the shared source files:

```bash
mkdir -p "$CUT_DERIVED_REFS/hg38" "$CUT_DERIVED_REFS/mm39"

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1}' \
    /opt/bioinformatics/references/hg38/hg38.chrom.sizes \
    > "$CUT_DERIVED_REFS/hg38/hg38.canonical_contigs.txt"

awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y)$/ {print $1}' \
    /opt/bioinformatics/references/mm39/mm39.chrom.sizes \
    > "$CUT_DERIVED_REFS/mm39/mm39.canonical_contigs.txt"

test "$(wc -l < "$CUT_DERIVED_REFS/hg38/hg38.canonical_contigs.txt")" -eq 24
test "$(wc -l < "$CUT_DERIVED_REFS/mm39/mm39.canonical_contigs.txt")" -eq 21

sha256sum "$CUT_DERIVED_REFS"/*/*.canonical_contigs.txt > \
    "$CUT_DERIVED_REFS/canonical_contigs.sha256"
```

The expected counts assume UCSC-style `chr` names and include autosomes, X, and
Y. If either count check fails, stop: inspect the chromosome-size and Bowtie2
index names rather than silently translating contigs.

## 7. Create the server configuration

Keep run configuration outside the installed release:

```bash
CUT_PROJECT=/path/to/project
mkdir -p "$CUT_PROJECT/config"

cp "$CUT_RELEASE_DIR/config/config.conf.template" \
    "$CUT_PROJECT/config/config.conf"
cp "$CUT_RELEASE_DIR/config/examples/cutrun_pe.csv" \
    "$CUT_PROJECT/config/samplesheet.csv"
```

For a human run, use these existing server resources:

```text
INDEX_HG38=/opt/bioinformatics/references/hg38/bowtie2/hg38
FASTA_HG38=/opt/bioinformatics/references/hg38/hg38.fa
CHROM_SIZES_HG38=/opt/bioinformatics/references/hg38/hg38.chrom.sizes
CANONICAL_CONTIGS_HG38=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/hg38/hg38.canonical_contigs.txt
GTF_HG38=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf
BLACKLIST_HG38=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.blacklist.bed
EFFECTIVE_GENOME_SIZE_HG38=2913022398
TSS_BED_HG38=
CCRE_BED_HG38=/opt/bioinformatics/ATACseq2tracks_shared/references/hg38/hg38.ccre.bed.gz
```

For a mouse run:

```text
INDEX_MM39=/opt/bioinformatics/references/mm39/bowtie2/mm39
FASTA_MM39=/opt/bioinformatics/references/mm39/mm39.fa
CHROM_SIZES_MM39=/opt/bioinformatics/references/mm39/mm39.chrom.sizes
CANONICAL_CONTIGS_MM39=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/mm39/mm39.canonical_contigs.txt
GTF_MM39=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/annotation.gtf
BLACKLIST_MM39=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.blacklist.bed
EFFECTIVE_GENOME_SIZE_MM39=2654621783
TSS_BED_MM39=
CCRE_BED_MM39=/opt/bioinformatics/ATACseq2tracks_shared/references/mm39/mm39.ccre.bed
```

Use environment commands for Picard, MACS3, and bigWig conversion. The CUT
configuration accepts a Picard executable, not the ATAC `PICARD_JAR` variable:

```text
PICARD_COMMAND=picard
MACS3_COMMAND=macs3
SEACR_COMMAND=/opt/bioinformatics/tools/SEACR/releases/1.3/SEACR_1.3.sh
```

Recommended first-pilot settings are:

```text
SPIKEIN_MODE=none
RUN_METAGENE=false
RUN_ATAQV_QC=false
GENERATE_ATAQV_VIEWER=false
RUN_MOTIF_ENRICHMENT=false
ENABLE_AUTOMATIC_CLEANUP=false
```

Set each row's `assay_profile` to `cutrun` or `cuttag`; there is no duplicate
run-wide assay setting. Do not mix genomes or PE/SE layouts in the first
validation run. For SE data, set `PEAK_CALLERS=macs3`; SEACR is PE-only.
Reference/blacklist paths belong only in config and are resolved from the row's
genome. Review control IDs carefully: a target must map to the correct
biological-replicate IgG/input/mock control, with no implicit replicate fallback.

## 8. Validate the installed software without processing reads

```bash
source /opt/miniconda/etc/profile.d/conda.sh
conda activate "$CUT_ENV"

bash "$CUT_RELEASE_DIR/tests/check_bash_syntax.sh"
python3 -m unittest discover \
    -s "$CUT_RELEASE_DIR/tests" -p 'test_*.py' -v

bash "$CUT_RELEASE_DIR/cutnrun2tracks.sh" --help
```

Verify the complete command surface used by the default first pilot:

```bash
for command_name in \
    python3 Rscript bowtie2 bowtie2-inspect samtools bedtools picard \
    trim_galore cutadapt fastqc multiqc macs3 bamCoverage \
    bedGraphToBigWig preseq computeMatrix plotProfile plotHeatmap plotFingerprint
do
    command -v "$command_name"
done

python3 -c 'import pandas, matplotlib, pyBigWig'
```

## 9. Validate one project configuration

The plan check validates configuration syntax, samplesheet structure, cohort
isolation, and control mapping without requiring the input/reference files or
bioinformatics tools:

```bash
bash "$CUT_RELEASE_DIR/cutnrun2tracks.sh" \
    --config "$CUT_PROJECT/config/config.conf" \
    --plan
```

Then perform the full tool, FASTQ, and reference audit:

```bash
bash "$CUT_RELEASE_DIR/cutnrun2tracks.sh" \
    --config "$CUT_PROJECT/config/config.conf" \
    --preflight-only
```

Inspect at least:

```text
<OUTPUT_DIR>/00_metadata/preflight_status.tsv
<OUTPUT_DIR>/00_metadata/software_versions.tsv
<OUTPUT_DIR>/00_metadata/resolved_config.tsv
<OUTPUT_DIR>/00_metadata/sample_manifest.tsv
<OUTPUT_DIR>/00_metadata/cohort_manifest.tsv
<OUTPUT_DIR>/00_metadata/reference_manifest.tsv
```

Do not continue if references from different assemblies are mixed, target and
control rows resolve incorrectly, or the resolved paths differ from the shared
server resources audited above.

## 10. Run staged real-data pilots

Use small representative data with at least two biological target replicates
and their correctly matched controls.

1. Run a PE CUT&RUN candidate with `SPIKEIN_MODE=none`,
   `RUN_METAGENE=false`, and cleanup disabled.
2. Confirm BAM integrity, fragment counts, duplicate policies, CPM tracks,
   MACS3/SEACR commands, consensus support, FRiP, fingerprints, differential
   raw-count inputs, report generation, and checkpoint resume.
3. Run CUT&Tag separately with its own assay profile.
4. Run SE separately with `PEAK_CALLERS=macs3` and confirm read—not inferred
   fragment—signal units.
5. Enable metagene only after its BED12 gene-set manifest has been built from
   the reused GTF/blacklist/chromosome resources and validated as described in
   [Metagene aggregate-signal module](06_metagene.md).
6. Enable automatic cleanup only after a deliberate failed-stage/resume and
   retention-policy test.

Keep each candidate's configuration, sanitized samplesheet, metadata manifests,
logs, environment export, workflow commit, and expected-versus-observed test
record.

## 11. Optional reuse of the existing dm6 composite references

Keep `SPIKEIN_MODE=none` for the first installation pilot. The existing ATAC
composite indexes can be reused only when the CUT samples actually contain dm6
spike material and the experimental spike-in stage/ratio metadata are valid.

Verify the existing resources:

```bash
for composite_prefix in \
    /opt/bioinformatics/references/composite/hg38_dm6/bowtie2/hg38_dm6 \
    /opt/bioinformatics/references/composite/mm39_dm6/bowtie2/mm39_dm6
do
    bowtie2-inspect -n "$composite_prefix" | grep -Fx 'dm6__chr2L'
    bowtie2-inspect -n "$composite_prefix" | grep -Fx 'dm6__chrX'
done

test -s /opt/bioinformatics/references/composite/hg38_dm6/reference_manifest.tsv
test -s /opt/bioinformatics/references/composite/mm39_dm6/reference_manifest.tsv
```

The shared dm6 chromosome sizes and blacklist use unprefixed dm6 coordinates,
whereas the composite indexes use `dm6__chr*`. The current CUT split BAM retains
those namespaced contigs. Create small derived coordinate files; do not pass the
unprefixed dm6 chromosome sizes directly to `SPIKEIN_CHROM_SIZES`:

```bash
mkdir -p "$CUT_DERIVED_REFS/dm6_namespaced"

awk 'BEGIN{OFS="\t"} $1 ~ /^chr(2L|2R|3L|3R|4|X)$/ {
  print "dm6__" $1, $2
}' /opt/bioinformatics/references/dm6/dm6.chrom.sizes \
  > "$CUT_DERIVED_REFS/dm6_namespaced/dm6.allowed.chrom.sizes"

cut -f 1 "$CUT_DERIVED_REFS/dm6_namespaced/dm6.allowed.chrom.sizes" \
  > "$CUT_DERIVED_REFS/dm6_namespaced/dm6.allowed_contigs.txt"

awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
  contig=$1
  if (contig !~ /^chr/) contig="chr" contig
  $1="dm6__" contig
  print
}' /opt/bioinformatics/references/dm6/dm6-blacklist.v2.bed \
  > "$CUT_DERIVED_REFS/dm6_namespaced/dm6.blacklist.bed"

test "$(wc -l < "$CUT_DERIVED_REFS/dm6_namespaced/dm6.allowed_contigs.txt")" -eq 6
sha256sum \
  "$CUT_DERIVED_REFS/dm6_namespaced/dm6.allowed.chrom.sizes" \
  "$CUT_DERIVED_REFS/dm6_namespaced/dm6.allowed_contigs.txt" \
  "$CUT_DERIVED_REFS/dm6_namespaced/dm6.blacklist.bed" > \
  "$CUT_DERIVED_REFS/dm6_namespaced/SHA256SUMS.txt"
```

For human plus dm6:

```text
SPIKEIN_MODE=dm6
SPIKEIN_REFERENCE_ID=dm6_release6
SPIKEIN_INDEX=/opt/bioinformatics/references/composite/hg38_dm6/bowtie2/hg38_dm6
SPIKEIN_FASTA=/opt/bioinformatics/references/composite/hg38_dm6/hg38_dm6.fa
SPIKEIN_CHROM_SIZES=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/dm6_namespaced/dm6.allowed.chrom.sizes
SPIKEIN_ALLOWED_CONTIGS=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/dm6_namespaced/dm6.allowed_contigs.txt
SPIKEIN_BLACKLIST=/opt/bioinformatics/references/cutnrun2tracks/0.2.0/dm6_namespaced/dm6.blacklist.bed
```

Use the `mm39_dm6` index and composite FASTA for mouse.

Current limitation: `SPIKEIN_BLACKLIST` is recorded in reference
provenance, but the current `spikein_batch.sh` does not apply it to the spike
BAM. Treat spike-in outputs as development/pilot results until that behavior is
explicitly corrected or accepted and tested. Also verify that the namespaced
spike control bigWigs open correctly before interpretation.

No shared *E. coli* composite reference was established by the prior audit. Do
not relabel dm6 resources for *E. coli*; build a new species-specific composite
only if the libraries require it.

## 12. Promote the validated release

Do not create the shared `current` link until the relevant CUT&RUN/CUT&Tag pilot
matrix passes. For a first installation, fail if a current link unexpectedly
exists:

```bash
test ! -e "$CUT_WORKFLOW_ROOT/current"
test ! -L "$CUT_WORKFLOW_ROOT/current"
ln -s "$CUT_RELEASE_DIR" "$CUT_WORKFLOW_ROOT/current"
readlink -f "$CUT_WORKFLOW_ROOT/current"
cat "$CUT_WORKFLOW_ROOT/current/VERSION"
```

### Required system-wide launcher

Install the repository-provided root-owned launcher. It pins the release and
main environment, exports a controlled PATH, and does not require interactive
Conda activation or inherit an ATAC/base environment.

```bash
sudo install -o root -g root -m 0755 \
  "$CUT_RELEASE_DIR/utilities/cutnrun2tracks_shared_launcher.sh" \
  /usr/local/bin/cutnrun2tracks

CUTNRUN2TRACKS_MAIN_ENV=/opt/miniconda/envs/cutnrun2tracks-0.3.0 \
  /usr/local/bin/cutnrun2tracks --version
```

When epic2 is enabled, create its sidecar and install the companion launcher:

```bash
mamba env create --prefix /opt/miniconda/envs/cutnrun2tracks-epic2-0.3.0 \
  -f "$CUT_RELEASE_DIR/environment.epic2.yml"
sudo install -o root -g root -m 0755 \
  "$CUT_RELEASE_DIR/utilities/epic2_shared_launcher.sh" /usr/local/bin/epic2
epic2 --version
```

The launcher does not grant access by itself: the workflow, environment,
references, SEACR release, and every parent directory must still satisfy the
read/execute policy above.

Users launch a complete run with one command:

```bash
nohup cutnrun2tracks --config /absolute/path/to/project/config/config.conf \
  > /absolute/path/to/project/cutnrun2tracks.nohup.log 2>&1 &
```

The external nohup file and internal raw/structured logs are both retained.

Retain the previous immutable release and environment during future upgrades.
Promote a new version by validating a new release directory and environment,
then deliberately switching `current`; never modify an already validated
release in place.

## 13. Apply and verify all-user access

First inspect traversal and permissions without changing them:

```bash
namei -l /opt/bioinformatics/workflows/cutnrun2tracks/current/cutnrun2tracks.sh
namei -l /opt/miniconda/envs/cutnrun2tracks-0.3.0/bin/python3
namei -l /opt/bioinformatics/tools/SEACR/releases/1.3/SEACR_1.3.sh
namei -l /opt/bioinformatics/references/hg38/bowtie2/hg38.1.bt2
namei -l /usr/local/bin/cutnrun2tracks
```

If the installation was created with a restrictive umask, correct only the
exact CUT-owned workflow, deployment-provenance, SEACR, and derived-reference
trees. Resolve and guard the paths before any recursive permission change:

```bash
cut_owned_runtime_roots=(
  "$CUT_RELEASE_DIR"
  "$CUT_DEPLOYMENT_DIR"
  "$SEACR_RELEASE"
  "$CUT_DERIVED_REFS"
)

for shared_root in "${cut_owned_runtime_roots[@]}"; do
    resolved_root="$(realpath -e -- "$shared_root")"
    case "$resolved_root" in
      /opt/bioinformatics/workflows/cutnrun2tracks/releases/* | \
      /opt/bioinformatics/workflows/cutnrun2tracks/deployment/* | \
      /opt/bioinformatics/tools/SEACR/releases/* | \
      /opt/bioinformatics/references/cutnrun2tracks/*)
        ;;
      *)
        printf 'REFUSING_UNEXPECTED_PATH\t%s\n' "$resolved_root" >&2
        exit 1
        ;;
    esac
done

chmod -R a+rX -- "${cut_owned_runtime_roots[@]}"
chmod -R go-w -- "${cut_owned_runtime_roots[@]}"
find "$CUT_RELEASE_DIR" -type f \
    \( -name '*.sh' -o -name '*.py' \) -exec chmod a+rx -- {} +

chmod a+rx -- "$CUT_WORKFLOW_ROOT" \
    "$CUT_WORKFLOW_ROOT/releases" "$CUT_WORKFLOW_ROOT/deployment"
chmod go-w -- "$CUT_WORKFLOW_ROOT" \
    "$CUT_WORKFLOW_ROOT/releases" "$CUT_WORKFLOW_ROOT/deployment"
```

The CUT Conda environment was created independently with `umask 022` and
`--copy`; verify it rather than applying a blanket recursive chmod. If its
ordinary-user access test fails, have the Conda administrator repair or recreate
that one CUT prefix according to the site's Conda permission policy.

Do not recursively chmod the complete `/opt/bioinformatics/references`,
`/opt/miniconda`, or ATACseq2tracks installation. Those are existing shared
assets outside the CUT deployment scope. If a parent directory lacks traversal
permission, have the server administrator grant the minimum required execute
permission or site-approved read-only ACL on that exact parent.

Verify access as a representative ordinary server user. Replace the placeholder
with a real username that does not own the installation and has no administrative
privileges:

```bash
CUT_TEST_USER=REPLACE_WITH_ORDINARY_SERVER_USERNAME
id "$CUT_TEST_USER"

sudo -u "$CUT_TEST_USER" test -r \
    /opt/bioinformatics/workflows/cutnrun2tracks/current/README.md
sudo -u "$CUT_TEST_USER" test -r \
    /opt/bioinformatics/workflows/cutnrun2tracks/current/cutnrun2tracks.sh
sudo -u "$CUT_TEST_USER" test -x \
    /opt/miniconda/envs/cutnrun2tracks-0.3.0/bin/python3
sudo -u "$CUT_TEST_USER" test -x \
    /opt/bioinformatics/tools/SEACR/releases/1.3/SEACR_1.3.sh
sudo -u "$CUT_TEST_USER" test -r \
    /opt/bioinformatics/references/hg38/bowtie2/hg38.1.bt2
sudo -u "$CUT_TEST_USER" test -r \
    /opt/bioinformatics/ATACseq2tracks_shared/references/hg38/annotation.gtf

sudo -u "$CUT_TEST_USER" /usr/local/bin/cutnrun2tracks --help

sudo -u "$CUT_TEST_USER" \
  /opt/miniconda/bin/conda run --no-capture-output \
  --prefix /opt/miniconda/envs/cutnrun2tracks-0.3.0 \
  python3 -c 'import pyBigWig; print("CUT environment access: OK")'

sudo -u "$CUT_TEST_USER" \
  /opt/miniconda/bin/conda run --no-capture-output \
  --prefix /opt/miniconda/envs/cutnrun2tracks-0.3.0 \
  bowtie2-inspect -n \
  /opt/bioinformatics/references/hg38/bowtie2/hg38 >/dev/null
```

Finally confirm the representative user cannot change the shared deployment:

```bash
if sudo -u "$CUT_TEST_USER" test -w "$CUT_WORKFLOW_ROOT" || \
   sudo -u "$CUT_TEST_USER" test -w "$CUT_RELEASE_DIR" || \
   sudo -u "$CUT_TEST_USER" test -w "$CUT_ENV" || \
   sudo -u "$CUT_TEST_USER" test -w "$SEACR_RELEASE" || \
   sudo -u "$CUT_TEST_USER" test -w "$CUT_DERIVED_REFS"
then
    echo 'FAIL: ordinary user can modify a shared CUT runtime tree' >&2
    exit 1
fi

echo 'PASS: ordinary user can execute but cannot modify cutnrun2tracks'
```

Repeat the read/execute checks for mm39 and dm6 resources before making those
branches available. If the server does not use `sudo`, perform the equivalent
checks with the site's supported account-switching mechanism.

## Installation acceptance checklist

- [ ] Production ATAC workflow and environment remain unchanged.
- [ ] Workflow, CUT environment, SEACR, references, and all parent directories
      are readable/traversable by every intended server user.
- [ ] Shared runtime trees and `current` are not writable by ordinary users.
- [ ] `/usr/local/bin/cutnrun2tracks --help` succeeds as a representative
      non-administrator user.
- [ ] CUT workflow is pinned to a full reviewed Git commit.
- [ ] CUT environment has its own prefix and explicit Linux package export.
- [ ] SEACR v1.3 commit and both script checksums are recorded.
- [ ] Host Bowtie2 indexes pass `bowtie2-inspect`.
- [ ] FASTA, chromosome sizes, canonical lists, GTF, blacklist, and cCRE are
      readable and assembly-matched.
- [ ] Repository syntax and synthetic tests pass inside the CUT environment.
- [ ] `--plan` and `--preflight-only` pass for the first project.
- [ ] Small real PE CUT&RUN and CUT&Tag pilots pass before shared promotion.
- [ ] SE and spike-in remain separately validated branches.
- [ ] Cleanup stays disabled until recovery and retention tests pass.
