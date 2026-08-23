# Contributing

Changes must preserve target-cohort isolation, raw-count differential inputs,
explicit control mapping, and the five independent track switches. Add a
synthetic regression test for every scientific branch change. Run
`tests/run_tests.sh` in the locked Linux environment before review.

Do not update shared references or the ATACseq2tracks installation from this
repository. Release manifests must record script and reference checksums.
