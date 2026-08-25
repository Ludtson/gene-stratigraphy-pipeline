---
name: Bug report
about: Report unexpected behavior in the classification or cross-check
title: ""
labels: bug
---

**Command used**

The exact command that was run, including flags (especially `--hog-filename`).

**Expected behavior**

What you expected to happen.

**Actual behavior**

What happened instead. If a candidate gene was misclassified, include the
species and gene ID, and whether `CLSGClassificationVerifier.py` reports it
as a violation.

**Environment**

- OS:
- Python version (`python3 --version`):
- OrthoFinder version used to generate the input:

**Input data**

Which `N?.tsv` was used and why (root-level run, or a specific clade of
interest — see the README's "Why `--hog-filename` might not be `N0.tsv`"
section). Whether `Orthogroups.tsv` and `Orthologues/` were available for
the cross-check, or `--skip-cross-check` was used.
