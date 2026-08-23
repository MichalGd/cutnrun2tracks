# Reusable metagene module

This directory is a pipeline-neutral Bash/Python module for aggregate-signal
profiles and heatmaps. It accepts normalized bigWigs through a track manifest
and versioned BED12 gene sets through a gene-set manifest. It never normalizes
signal values.

The public entry point is `run_metagene.sh`. It validates manifests, constructs
one task per sample, gene set, and mode, then runs deepTools `computeMatrix`,
`plotProfile`, and `plotHeatmap` headlessly. See `../../docs/06_metagene.md` for
the schemas, reference-building commands, configuration, and output contract.

`prepare_metagene_reference.py` and `build_hpa_reference_subset.py` are offline
reference-preparation utilities. Downloaded GTF, HPA, and BioMart inputs must be
cached and versioned by the caller; routine plotting does not access the
network.
