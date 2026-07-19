# Simulation Verification

This directory contains the local simulation-verification workspace.
The full simulation corpus is large, so it is stored as a release artifact at
`FBartos/bfpwr@sim-corpus-v1` rather than committed to the package repository.

## Quick Start

From the repository root, run:

```sh
Rscript simulations/scripts/build_verification.R
```

This is the normal one-command workflow. It downloads the prepared fixture
bundle if needed, recomputes package-side references with the current checkout,
runs the validation checks, and writes:

```text
simulations/simulation-verification.pdf
```

Use this before release or before reviewing simulation-verification changes.
The raw Monte Carlo simulations are not rebuilt by this command; those are
cluster jobs and are only rerun when the simulation corpus itself changes.
The checked package cases are read from the committed verification manifest by
default, so a normal report build should not rewrite tracked files.

## What The Build Does

The build script performs these steps:

1. Downloads and unpacks the release fixture bundle.
2. Reads the committed package-verification case manifest.
3. Recomputes package probabilities, sample-size searches, and timing
   diagnostics against the current package code.
4. Validates the recomputed results against the fixture summaries.
5. Renders the PDF report.

Each recomputation writes `package-reference-provenance.csv` and records the
same package Git revision in the search-comparison bundle. A release result
should name a clean package checkout; a dirty package tree means the result is
diagnostic and cannot be attributed to the recorded commit alone.

If a prebuilt search-validation bundle is present, the build refreshes it. If it
is missing, the build creates the search-validation bundle from the package
manifest so the default command still completes end to end.

## Useful Commands

Validate an already materialized corpus without downloading or rendering:

```sh
Rscript simulations/scripts/run_verification.R --corpus-root simulations/corpus/v1
```

Regenerate the PDF from existing validation outputs:

```sh
Rscript simulations/scripts/build_verification.R --skip-fetch --skip-recompute
```

Refresh the committed verification manifest after intentionally changing the
case-selection rules:

```sh
Rscript simulations/scripts/build_verification.R --refresh-manifests --skip-report
```

Compare timings from two fixture-reference recomputations performed in the
same machine session:

```sh
Rscript simulations/scripts/compare_verification_timings.R \
  --baseline-fixture-dir <baseline-output> \
  --candidate-fixture-dir <candidate-output> \
  --max-total-ratio 1.10
```

The comparison aligns explicit verification-case keys and ignores functions
whose complete baseline workload took less than one second. It is intended for
matched development runs; timings from different machines or unrelated
sessions are descriptive only.

Download only the fixture-test assets used by package fixture tests:

```sh
Rscript simulations/scripts/download_corpus_release.R \
  --corpus-root simulations/corpus/v1 \
  --asset-set fixture-tests
```

## Local Requirements

The one-command build needs:

- R packages used by the simulation tooling, including `knitr`.
- A `tar` executable that can extract `.tar.zst` files.
- A LaTeX installation available to `knitr::knit2pdf()`.

## Repository Hygiene

Commit the scripts, manifests, documentation, and rendered verification PDF.
Do not commit the downloaded corpus, release archives, logs, generated figures,
coverage summaries, generated registry CSVs, or intermediate TeX files. These
local outputs are ignored by `.gitignore`:

```text
simulations/corpus/
simulations/downloads/
simulations/logs/
simulations/figure/
simulations/registry/manifests/*.csv except package-verification-cases.csv
```

## Report Contents

The report is meant to answer practical release questions:

- Do simulated fixture summaries still agree with the package functions?
- Do sequential probability and expected-sample-size checks still pass?
- Do fixed and sequential sample-size searches still find the expected cases?
- Which package calls are slow enough to deserve attention?

Some rows are diagnostic rather than pass/fail checks. In particular, fixed
`t` references include known approximation diagnostics, and sequential binomial
fixtures are checked for fixture integrity because the package does not expose a
sequential binomial probability API.
