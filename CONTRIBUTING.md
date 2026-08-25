# Contributing

This project is maintained by the Casola Lab (Ecology & Conservation Biology
Program, Texas A&M University) as part of the LINGUA comparative genomics
framework. It is Stage 2 (gene stratigraphy) — see the README for how it
relates to Stage 1 (protein preprocessing) and Stages 3–4 (homolog exclusion,
synteny validation).

## Getting set up

Requires bash, Python 3.7+, and OrthoFinder (only needed if you're running
the generation step; the classification step only reads OrthoFinder's
output). No third-party Python dependencies. See the README's Requirements
section for the recommended conda install of OrthoFinder.

```bash
git clone https://github.com/Ludtson/gene-stratigraphy-pipeline.git
cd gene-stratigraphy-pipeline
```

## Verifying a change

Before proposing a change to the classification logic
(`OrthoFinderDataProcessor.py`, `OrthoFinderN0GeneSorter.py`,
`OrthoFinderTreeParser.py`, `OrthologyCrossChecker.py`), run
`CLSGClassificationVerifier.py` against the output to confirm it still
reports zero violations:

```bash
python bin/CLSGClassificationVerifier.py \
  --hog-tsv <path to the N?.tsv used> \
  --results-dir <your output dir> \
  --report violations_report.tsv
```

A change that introduces violations here needs to be understood and fixed,
not just noted, before it's proposed. If a change intentionally alters what
counts as internally consistent, explain why in the pull request.

## Code style

- One class per file, matching the file's name (`CSVFileMerger.py` →
  `class CSVFileMerger`, etc.) — see the README's Components section for the
  existing pattern. New scripts should follow it, not introduce a mix of
  styles.
- Locate expected input files by searching under a results directory
  (`_find_file`/`_find_dir`'s pattern) rather than requiring exact paths,
  matching the existing modularity.
- Standard library only; a new third-party dependency needs a strong reason.
- Match the existing logging conventions rather than introducing a new
  style.
- Missing/not-applicable values in output tables are written as the string
  `"NA"`, not an empty string — keep this consistent across any new output.

## Submitting changes

1. Open an issue first for anything beyond a small fix.
2. Keep pull requests focused on one change.
3. Include the verifier's output (or an explanation of why it's expected to
   change) in the PR description.
