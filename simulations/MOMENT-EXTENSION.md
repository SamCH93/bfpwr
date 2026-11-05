# PR 8 moment simulation extension

The extension adds `less` and `greater` to the existing 35 moment-prior cases.
It preserves the 15 generating conditions, 10,000 original trajectories per
condition, look grids, thresholds and schedules. The complete grid contains
146,400 fixed and 222,570 sequential comparisons, 4,332 stopping distributions
and 4,200 sample-size searches. No new trajectories are generated.

## Calculation and provenance

Use `ANALYSIS_RSCRIPT` locally. The PBS wrapper retains R 4.5.1 and the existing
cluster library, and accepts `SOURCE`, `INPUT`, `OUTPUT`, `VERSION`, `RELEASE`,
`WORKERS` and `PHASE`, plus `REUSE_FROM`, `BF_SOURCE_ROOT` and optional
`OFFICIAL_METADATA` for import. Its defaults write to a new `moment-v2` destination.
Do not change the frozen calculation sources during a run.
Run plan, import, summaries and assembly with the same R version, platform and
numerical libraries. In particular, MetaCentrum R 4.5.1 needs its own fresh output
directory and a new plan/import there; it cannot resume the Windows R 4.6 preflight.

```sh
Rscript simulations/scripts/check_moment_extension.R
Rscript simulations/scripts/run_moment_extension.R --phase plan --input-corpus INPUT --output OUTPUT
Rscript simulations/scripts/run_moment_extension.R --phase chunks --input-corpus INPUT --output OUTPUT --workers 40 --worker 1
Rscript simulations/scripts/run_moment_extension.R --phase summaries --input-corpus INPUT --output OUTPUT --workers 40 --worker 1
Rscript simulations/scripts/finalize_moment_extension_chunks.R --input-corpus INPUT --output OUTPUT
Rscript simulations/scripts/run_moment_extension.R --phase assemble --input-corpus INPUT --output OUTPUT
Rscript simulations/scripts/validate_moment_extension.R --output OUTPUT
Rscript simulations/scripts/package_moment_extension.R --output OUTPUT --release RELEASE
```

Planning writes separate BF-dependency and reference identities, dynamically
records the Git revision, and freezes cases, designs, original input indexes,
schedules, thresholds, targets, numerical defaults, RNG kind, runtime and source
snapshots. Source manifests record physical SHA256 hashes. Computational
compatibility normalizes line endings and trailing empty lines; raw input,
output, snapshot and archive integrity checks always use physical SHA256.
The reference snapshot contains the package source and simulation/report tools.
Its cache key includes package R code, DESCRIPTION/NAMESPACE, simulation R and
registry modules, and the numerical runner, support, finalization and validation
scripts. Report text, rendering/preparation layout, PDF assembly, README, NEWS
and help pages can change without staling numerical results. Every saved source
file, including presentation files, still has a strict physical integrity check.
On PBS, the source revision is exported before copying files to scratch.

Workers use atomic lock directories and retain error files. Completed work is
reused only under compatible identities. Fresh BF identities also record R,
platform and stats versions, so fresh chunks from different runtimes cannot mix.
Inspect errors and confirm the worker
has finished before removing a failed lock. Input paths in chunk metadata are
relative to `provenance$source_corpus`. Finalization verifies the actual original
input files against their SHA256 indexes and writes standard BF indexes. Imported
statuses retain the original input verbatim in `historical_input`; its path base
labels an absolute historical path or the historical source corpus as appropriate.

## Verified reuse of v1 BF chunks

BF evaluation and probability references have different dependencies.
The unchanged `nmbf01.R` and BF-forwarding source, exact cases, original input
hashes and BF output hashes must be proven before importing old BF chunks.
Unknown BF source changes require regeneration. Import never relabels historical
BF work as a new computation; it retains the original identity and provenance,
stores the historical BF source and a separate compatibility proof, and leaves
all candidate reference calculations to be run afresh.
An explicit verified import may cross runtimes: the historical identity and
provenance describe the old computation, and the new compatibility proof binds
the verified data to the current run. This does not authorize fresh BF evaluation
in that imported output directory.

```sh
Rscript simulations/scripts/download_moment_extension.R --version moment-v1 --output OLD --with-bf-chunks
Rscript simulations/scripts/run_moment_extension.R --phase plan --input-corpus INPUT --output NEW
Rscript simulations/scripts/run_moment_extension.R --phase import --input-corpus INPUT --output NEW --reuse-from OLD --bf-source-root OLD/sources/bf-evaluation-historical --official-metadata DOWNLOADS/bfpwr-sim-corpus-v1-moment-metadata.tar.zst
Rscript simulations/scripts/run_moment_extension.R --phase summaries --input-corpus INPUT --output NEW --workers 40 --worker 1
```

Then finalize, assemble, validate and package `NEW` as above. Assembly compares
every empirical column with v1 for exact equality. It also loads the independent
`moment-reference-check.R` module, freezes `bundle$designs`, and saves
`numerical-diagnostics.rds` and `numerical-diagnostics.csv`.
Original v1 results and release assets remain immutable.
`--official-metadata` is optional. It verifies the published metadata archive's
raw SHA256 (`2c5906d4dd598074a75f9d641b1e8183a1b8b99418782d58efbdfd3a1557e4a9`),
checks archive members before private extraction, and requires the complete
imported case/chunk keys, original input hashes and actual BF output hashes to
match its authoritative indexes. Only this proof authorizes the original
published BF archive links. A generic compatible v1-shaped import remains valid
without this option and packages its copied chunks into new v2 BF archives.

## Comparisons and timings

For each case, threshold and schedule family, only the longest compatible
schedule is integrated. Its cumulative prefixes supply every existing
comparison schedule. Prefix EN and VarN are recomputed with the existing package
moment helper, assigning unresolved mass to that prefix's own final look.
The dense 50/100-look curves also supply the sequential search reference scan.
All 4,200 actual `nnmbf01()`/`nbf01seq()` target checks remain, and the fixed
integer-grid reference scan remains independent of the fast solver.

The v1 MC screening rules are retained. Moment tables record MCSE, observed
maximum N, the actual tolerance, signed error (package minus empirical), relative
error and standardized error. Tolerances are the larger of five MCSEs and
0.002 times observed max N (EN), or observed max N squared (VarN).
Relative errors with a zero empirical denominator are NA. EN has sample-size
units; VarN has squared sample-size units. A missed sequential first candidate
fails validation; the lower-bound exception applies only to fixed searches.
The [manuscript check](one-sided-moment-manuscript-check.md) records a fixed
search example with a later crossing and checks the exhaustive first crossing.

`bundle$timing_schema_version = 2` identifies actual package-call timings,
including solver calls, with operation, phase, schedule, threshold, looks and
rows covered. Reused rows are separate `bundle$reuse` records.
`bundle$case_timings` records full case elapsed time, including preflight,
empirical summaries, prediction and search phases. Worker startup and final
serialization are outside that interval. Concurrent case times must not be
summed to claim run wall time; use their observed timestamp span.

MC/search screening and independent-reference convergence are separate checks
and CSV files. Convergence of the 64/128-node reference is not a pass for the
package's default numerical accuracy. Package errors and reference
nonconvergence remain visible in the diagnostic results.

## Release and reports

The v2 packager creates new results, metadata, source, checksum and CSV-manifest
assets. Every run checks the complete planned BF chunk set, per-case indexes,
and chunk hashes against their indexes and statuses. Imports also have to
match their saved import proof. For imports proven against the official metadata,
an additional check against that proof authorizes references to the original v1
BF archives by their published SHA256 hashes, avoiding duplicate 1.5 GB uploads.
Other imports and fresh calculations get their own v2 BF archives. It never replaces v1
assets and rejects incompatible payload reuse in an existing release directory.
Fresh packaging needs no manually copied identities or source tar commands.
The release checksum list includes the manifest; the downloader verifies the
manifest and every selected archive before extraction, and accepts only regular
files and directories with safe relative paths. Archive file lists and payload
identities share radix ordering and exclude `.error`, `.tmp` and `.lock` evidence;
those files remain available locally.

Download and extract the verified baseline archive as shown in
[One-sided moment extension reports](README.md#one-sided-moment-extension-reports).
Use its extraction directory as `BASELINE`.

```sh
Rscript simulations/scripts/download_moment_extension.R --version moment-v2 --output OUTPUT
Rscript simulations/scripts/render_moment_extension.R --baseline-root BASELINE --moment-root OUTPUT
```

Add `--with-bf-chunks` for full BF data. Result-only downloads reproduce the
reports. `--release-base` and `--downloads` support release mirrors and private
download destinations. Source reconstruction is under `OUTPUT/sources/reference`.
The renderer reads the saved verification bundle directly. The optional
fixture export reads from `--results`, which defaults to
`simulations/corpus/moment-v2`. It writes to a sibling directory with
`-report` appended to the results directory name by default. To choose a
different destination, pass `--output`. The export rejects in-place use.

```sh
Rscript simulations/scripts/prepare_moment_extension_report.R --results OUTPUT --output PREPARED
```

PBS assembly saves `known-search-limitation.csv` with the other computational
diagnostics in `OUTPUT` for packaging.

The full PDF pools existing and new priors in the z-test tables, settings,
and fixed and sequential results. Its fixed z-test search plot pools the
existing Normal searches with the new Moment searches. The [incremental review
excerpt](https://github.com/FBartos/bfpwr/releases/download/sim-corpus-v1/simulation-verification-moment-v2-incremental.pdf)
selects the affected sections and shows only the new Moment+ and Moment-
results, searches, and timings. Both scopes use the approved table
and section structure. Additional diagnostics stay in machine-readable
artifacts only. Both PDFs render natively from complete cached
baseline and moment inputs without reference calculations. Independent cohort
provenance and input/output hashes are saved alongside the PDFs. Rendering
requires knitr and LaTeX.

Each invocation writes `report-input-provenance.rds` for the requested scope.
Use `--scope incremental` to reproduce the single-output sidecar published as
`incremental-input-provenance.rds`. The default `--scope both` sidecar records
both output hashes. The historical full report retains its original release
provenance. The renderer's RDS contains the following fields.

- `scope` and the absolute `baseline_root` and `moment_root`.
- `baseline`, `normal`, and `moment` provenance from the saved inputs.
- `baseline_input_sha256`, `moment_input_sha256`, and `source_sha256` as
  named hash vectors keyed by paths relative to their respective roots.
- `renderer` with the R version, platform, library paths, and source root.
- `output_sha256` keyed by the PDF filenames rendered in this invocation.

Focused engineering checks are `check_moment_runner_regressions.R` and
`check_moment_pipeline_regressions.R`. The latter uses a private small corpus
and the shipped reference module.
To check the export with a saved v2 bundle, run
`Rscript simulations/tests/test-moment-report-export.R --results OUTPUT`.
The check copies the required RDS inputs to a temporary directory.
To check incremental report selection with the saved inputs, run
`Rscript simulations/tests/test-moment-report-scope.R --baseline-root BASELINE --moment-root OUTPUT`.
