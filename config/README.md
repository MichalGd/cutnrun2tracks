# Configuration files

Copy `config.conf.template` to an untracked `config.conf`, then set every input,
output, and reference path used by the selected genome. The parser accepts only
the keys present in the template plus assembly-specific reference keys matching
`INDEX_*`, `FASTA_*`, `CHROM_SIZES_*`, `CANONICAL_CONTIGS_*`, `GTF_*`, `BLACKLIST_*`, `TSS_BED_*`, or
`CCRE_BED_*`.

`RUN_METAGENE=true` additionally requires `METAGENE_GENE_SET_MANIFEST` to name
a tab-separated gene-set manifest prepared with the utilities under
`common/metagene/`. The plotting stage consumes existing CPM or spike-in
bigWigs and does not perform normalization.

The examples contain illustrative `/data` and `/refs` paths and are not directly
executable. `cutrun_se.csv` requires `PEAK_CALLERS=macs3`.
