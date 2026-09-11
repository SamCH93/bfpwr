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
To refresh the additional one-sided normal scenarios as well, use the complete
prediction refresh described below. The report refuses to combine bundles with
different or unrecorded numerical accuracy settings.

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

Compute the one-sided normal z-test validation on the original 10,000
trajectories per generating condition, then refresh the PDF and preview plots:

```sh
Rscript simulations/scripts/download_corpus_release.R --asset-set designs
Rscript simulations/scripts/verify_one_sided_z.R
```

This adds 70 realistic prior/design/direction combinations and selected long
schedules. It does not simulate new observations or increase replication counts.
The independent Bayes factors, probability and stopping-time summaries,
sample-size search checks, and provenance are saved under
`simulations/corpus/v1/one-sided-z/`. Preview images are written to
`simulations/figure/one-sided-z-curves.png` and `one-sided-z-agreement.png`.
The report combines this result bundle with the existing z-test diagnostic
plots and prior tables when present. Stopping curves remain a separate preview
and are not included in the PDF. Use `--skip-report` to
compute without rendering, or `--package-root <clean-checkout>` to attribute
package predictions to a clean revision while preserving local work.
The current validation exposes accuracy limitations in the default sequential
integration for some dense schedules. Selected finer-grid/backend calculations
are saved in `convergence.csv`; they do not replace the default predictions.
`probability_failures.csv` identifies comparisons failing the Holm-adjusted
binomial diagnostic, and `moments.csv` records stopping-time moment checks.

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

## Refreshing package predictions

Recompute the report's package predictions using the package's numerical defaults,
while retaining the existing 10,000 simulated replications:

```sh
Rscript simulations/scripts/build_package_verification_manifest.R
Rscript simulations/scripts/refresh_simulation_grid.R --workers 6
```

The refresh writes separate results under
`simulations/corpus/v1/predictions-default-ngrid-10000/`. It recomputes the fixed
and sequential references, the required sample-size searches, and all one-sided
z-test predictions, searches, and Bayes factor spot checks. Sequential t-test
references include H1, H0, and inconclusive probabilities at every selected look,
expected sample sizes, and the four 20-look cases shown in the original report.
Exploratory search rows remain explicitly marked as unevaluated.

The validation calls omit numerical tuning arguments, including `ngrid`,
`tail.eps`, and `tail.nquad`, so they use the shared defaults in
`package/R/numerical-defaults.R`. The current settings are 10,000 integration
points, root tolerance `1e-8`, BF integration tolerance `1e-8` with up to 1,000
subdivisions, predictive tail mass `1e-6`, and 512 fallback quadrature nodes.
These values are recorded in the result provenance. The old simulation-only
`bfpwr.sim.ngrid` option and `--ngrid` override no longer control verification.
Explicit accuracy comparisons in the separate convergence diagnostic remain
labelled with their settings and do not replace the default predictions.

Use `--results-root` for a different output directory.
`--skip-report` runs the calculations without rendering, and
`--assemble-only` assembles completed per-case results without recomputing them.
Assembly requires every recorded accuracy setting to match the package defaults;
otherwise, rerun without `--assemble-only`.
The report reads only results directly in the selected directory, so older
subdirectory caches cannot silently fill missing comparisons.
Validation failures are retained in the diagnostic files and report; the refresh
finishes reporting them before returning a failing exit status.

To render the completed results again:

```sh
Rscript simulations/scripts/render_simulation_verification_background.R \
  --results-root simulations/corpus/v1/predictions-default-ngrid-10000
```

The grid refresh fixes the integration RNG kind to
`Mersenne-Twister/Inversion/Rejection`. Recorded elapsed times come from concurrent
workers and are descriptive, rather than isolated performance benchmarks.
