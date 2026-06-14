# In-depth review: bfpwr package + PR #6 ("Add sequential sample-size search")

- **Date:** 2026-06-10
- **Reviewed at:** branch `seq-n-find`, commit `e670fcf` (= PR #6 head; base branch `gsd`)
- **Method:** Seven parallel review agents, each owning one domain (search engine, sequential helpers, seq wrappers, sequential probability functions, t-test family, z-test/normal-moment family, binomial family, cross-cutting tests/packaging). Agents re-derived the mathematics independently, cross-checked against Monte Carlo simulation, `BayesFactor`, brute-force enumeration, and 256-bit arithmetic (`Rmpfr`), and ran the full test suite (default + extended). The three highest-impact findings were additionally re-verified independently of the agents. Findings are labeled **confirmed** (reproduced by running code), **high**, or **medium** confidence, and **NEW-IN-PR** vs **PRE-EXISTING** (on `gsd`/`main` before this PR).
- **Original review scope:** findings only; the 2026-06-13/14 follow-up below
  records later implementation changes.

## 2026-06-13/14 follow-up: one-sided t adaptive search update

The PR6 follow-up replaces the old one-sided adaptive t-search terminal cap
(`|t| <= 256`) with a predictive-tail cutoff controlled by `tail.eps`
(default `1e-3`) in the fixed-n and sequential t power functions. This makes
the adaptive stopping rule a probability guarantee under the relevant
predictive distribution, rather than an arbitrary t-statistic magnitude that
behaves differently across sample sizes and designs.

Implementation comment for PR6:

- `tcrit()` remains distribution-agnostic: it receives a finite
  `search_limit`, not `tail.eps`. The fixed-n and sequential callers know the
  predictive mean/variance of the t-statistic, so they convert `tail.eps` into
  the appropriate positive and negative t cutoffs before calling `tcrit()`.
- For fixed-n power calculations, if no boundary is found before the searched
  tail has remaining predictive mass at most `tail.eps`, the function returns
  the boundary-free 0/1 approximation implied by the searched side and warns
  with the omitted-tail bound. This approximation is now tied directly to a
  user-visible probability tolerance.
- For sequential one-sided t designs, the same `tail.eps` is used for each
  adaptive boundary search. It is not split across looks or boundaries. If a
  boundary is unresolved at the cutoff, the omitted marginal mass for that
  boundary search is bounded by `tail.eps`; users can lower `tail.eps` or pass a
  numeric `trange` to force exact boundary searches over that interval.
- One-sided searches now scan only the mathematically expected direction. This
  removes the PR5/PR6 opposite-side fallback while preserving valid wrong-tail
  H0 roots when the sign at the origin implies that side.
- The adaptive search now marks a tail cutoff only after reaching a finite,
  certifying endpoint. Non-finite endpoint behavior remains a numerical search
  failure instead of being converted into a valid-looking 0/1 approximation.
- Sequential region construction now treats an unresolved earlier H1 boundary
  as unbounded continuation, rather than allowing a bare `NaN` to zero all later
  paths. A missing final stopping boundary remains an empty stopping event.
- Sequential boundary diagnostics now carry explicit statuses such as
  `ok`, `impossible`, `tail_cutoff`, and `search_failed`, so mathematically
  empty regions and unresolved numerical searches are no longer represented by
  the same bare `NaN` state.
- Benchmark check against the previous PR6 head (`e670fcf`): full extended
  tinytest runtime was effectively unchanged (current 51.7 s vs baseline
  51.2 s). Regular adaptive outputs matched. In an extreme sequential case,
  current adaptive differed from explicit `trange = c(-50, 50)` by `2.95e-5`
  for final H1 probability and `3.52e-4` for final H0 probability, within
  `tail.eps = 0.001`; the explicit range took about 10-18 s versus about
  0.03 s for adaptive search.

## 2026-06-14 follow-up: search policy and final hardening

Current PR6 search-risk status is tracked in
`PR5_PR6_SEARCH_RISK_REVIEW.md`. The search-risk follow-up addressed the
original review's active search-policy issues:

- Failed sequential t boundary searches no longer enter integration as
  ordinary empty `NaN` regions; boundary results carry explicit status before
  region construction.
- Sequential sample-size search catches only typed candidate-level numerical
  invalidity. Structural evaluator errors and malformed evaluator returns now
  propagate.
- `search = "exhaustive"` now scans the full candidate domain before relying
  on adaptive bracketing, skips transient invalid candidates, and stops at
  terminal invalid candidates.
- Fixed-`n` sequential wrapper schedules now round-trip searched increment
  schedules by respecting `nrange[1]` as the default first look when `minN` is
  missing. The t wrapper also allows non-decreasing group-2 schedules generated
  by small allocation ratios.

Do not read the historical "Suggested merge gate" section below as the current
status for C2/M1/M2/M4. It remains useful as the original review record. The
remaining non-search and pre-existing findings in that review still need a
separate decision if they are in scope for this PR.

Final hardening checks run on 2026-06-14:

- Full extended tinytest suite: 660 results passed.
- `R CMD check` on a built source tarball with `--no-manual --ignore-vignettes`:
  `Status: OK`.
- `R CMD check` with `--no-build-vignettes` reported only the existing vignette
  profile warnings about `vignettes/` files without generated `inst/doc`
  outputs.

## Executive summary

The mathematical core of the package — including everything PR #6 adds — is in very good shape. The sequential covariance structure, first-passage stopping-region decomposition, BF↔boundary inversions, the new log-scale tail computations, the new exact wrong-tail t-BF path, and the binomial/dirbf01 reworks were all re-derived and verified numerically (details in §5). The search engine in `seqsearch.R` is carefully built, with deterministic caching and honest certification semantics.

The defects that matter cluster in **failure-mode handling, not formulas**: two paths return **silently, drastically wrong probabilities** (one pre-existing in `ptbf01`, one substantially new in `ptbf01seq`), the new search layer **swallows errors** and misattributes them to `nrange`, and the bracketing phase assumes power is eventually increasing in `n`, which is false for `target = "h0"` designs. Separately, there is a **process risk**: virtually all numerical-correctness tests are gated behind `BFPWR_RUN_EXTENDED_TESTS` (off by default, no CI exists to set it), so CRAN and default checks exercise almost none of the math validated here.

Suggested merge gate: fix **C1, C2, M1–M4** before merge; M5–M9 and the minor items can follow in this PR or a fast-follow.

---

## 1. Critical — silent wrong results in realistic use

### C1. `ptbf01` two-sided guard assumes max BF01 is at `t = null`; false for informed (asymmetric) priors → silently returns power 0/1
- **File:** `package/R/ptbf01.R:110-123` — **PRE-EXISTING** — **confirmed**
- **Status after 2026-06-14 change:** **Fixed.** The early-exit guard keeps
  the null-point shortcut only for centered two-sided priors. For shifted
  informed priors where `BF01(null) < k`, `ptbf01()` now performs a bounded
  maximum over the searched interval and uses that maximum as the split point
  for the two root searches. The regression is pinned against explicit
  bracketing roots in `test-ptbf01.R`.
- **Problem:** The early-exit guard for "k above the maximum attainable BF01" evaluates BF01 only at the null. That is the argmax only for priors symmetric around the null. For an informed prior (`plocation ≠ 0`, two-sided), the max of BF01(t) is at some t ≠ 0 and can far exceed BF01(0).
- **Reproduction (re-verified):** `ptbf01(k = 30, n = 25, plocation = 0.6, pscale = 0.15, pdf = 10, type = "one.sample", alternative = "two.sided", dpm = -0.5, dpsd = 0, lower.tail = FALSE)` returns **0**; Monte Carlo truth ≈ **0.95**. BF01(0) ≈ 18.6 while max BF01 ≈ 132 (argmax ≈ −0.46). No warning is emitted.
- **Why it matters:** Informed-prior two-sided designs with k > BF01(0) silently get power 0 or 1. PR #6 fixed exactly this flaw in `tcrit()` (`seqhelpers.R:928-944`, bounded `optimize` over the search interval) but did not port the fix to `ptbf01`.
- **Fix:** Replace the point evaluation at `null` with the same bounded maximization of `rootFun` used in `tcrit()`, or drop the guard and rely on the existing `uniroot(extendInt = ...)` fallbacks which already handle the zero/one-root cases.

### C2. `ptbf01seq`: failed t-boundary searches are silenced and NaN boundaries silently zero stage probabilities → drastically wrong sequential power with no diagnostic
- **Files:** `package/R/ptbf01seq.R:134-146` (warning handler — **NEW-IN-PR**), `package/R/ptbf01seq.R:147-168` + `package/R/seqhelpers.R:160-162, 280-288, 390-397` (NaN-region semantics — **PRE-EXISTING**) — **confirmed**
- **Problem (two interacting defects):**
  1. The new `evalTcrit` handler counts only "Adaptive t critical-value search reached…" warnings but calls `invokeRestart("muffleWarning")` for **every** warning, so `tcrit`'s "Numerical problems finding critical value" / "maximum BF is less than k" warnings vanish entirely.
  2. Downstream, `.bfseq_intstage_sum` maps any region containing NaN to probability 0. That is the intended encoding for a *mathematically nonexistent* H0 boundary (`zk0 = NaN`, handled via `H0nan` → bound extends to ∓Inf), but a NaN **H1** boundary is always a numerical search failure (for one-sided BFs the H1 boundary provably exists), and it zeroes the stopping probability of that look — and, when it sits in a continuation column, of **all later looks** (both H1 and H0).
- **Reproduction (re-verified):** `ptbf01seq(k1 = 1/10, k0 = 10, n = c(10, 3000), alternative = "greater", dpm = 0.5, dpsd = 0, trange = c(-2, 3))` — the look-2 H1 boundary (≈3.02) lies outside `trange`, `tcrit` returns NaN, and the function reports `cumpH1 = (0.040, 0.040)` with **zero warnings**. The fixed-design probability at n = 3000 is **1.0**. Calling `tcrit` directly on the same inputs *does* warn — the wrapper eats it.
- **Fix:** (a) Move `invokeRestart("muffleWarning")` inside the `if` that matches the adaptive-limit message (or aggregate and re-emit all distinct messages); (b) in `ptbf01seq`/`pbf01seq`, after computing `zk1`, `stop()` (or warn loudly and extend the bound conservatively) when `any(is.nan(zk1))` instead of feeding NaN into region construction; mirror the `H0nan` treatment if a non-error path is preferred.

---

## 2. Major

### M1. `seqsearch.R` bracketing assumes power is eventually increasing in n; unimodal curves (`target = "h0"`) yield spurious `NaN` with a misleading warning — and `search = "exhaustive"` never actually scans in that case
- **File:** `package/R/seqsearch.R:149-151, 286-327, 253-276, 478-489` — **NEW-IN-PR** — **confirmed**
- **Problem:** The bracket grows geometrically (2, 4, 8, …) and gives up if no probed point meets the target, warning that `nrange[2]` "leads to lower power than specified". But for `target = "h0"` with a design prior putting mass on nonzero effects, P(BF01 ≥ k0) is **unimodal** in n (rises to a peak < 1, then → 0), so an above-target window between doubling points is missed entirely. Reproduced: `k1 = 1/6, k0 = 6, psd = 1, dpm = 0.3, dpsd = 0, target = "h0", looks = 1, power = 0.055` — true solution window n ≈ 83–110 (peak power 0.060 at n = 96), but the bracket probes 64 (power 0) and 128 (0.051), returns `NaN`, and the warning misdirects the user toward raising `nrange[2]`. Crucially, the **exhaustive first-crossing scan only runs after a successful bracket**, so `search = "exhaustive"` returns the identical `NaN` — violating its documented contract exactly when it matters.
- **Fix:** When the bracket fails under `search = "exhaustive"`, fall back to a true scan of the candidate range (or a golden-section peak search — O(log) for unimodal curves). For `search = "adaptive"`, detect a decreasing tail after a positive maximum in the bracket trace and warn "power appears non-monotone; maximum observed power X at n = Y" instead of the upper-bound message. Document the eventually-increasing assumption in `nbf01seq`/`ntbf01seq`.

### M2. `seqsearch.R` evaluation loop swallows all errors (including `...` typos) and misreports them as an `nrange` problem
- **File:** `package/R/seqsearch.R:188-199` (manifests in all four new wrappers); related abort-path issue at `seqsearch.R:306-356` — **NEW-IN-PR** — **confirmed**
- **Problem:** `evalN` wraps every design evaluation in `try()` and converts any error to `criterion = NA`. The search then returns `NaN` with the warning "upper bound of sample size search range ('nrange') leads to Power = NaN". Reproduced: `nbf01seq(..., bogusarg = 1)` → `NaN` with that warning, while `pbf01seq(..., bogusarg = 1)` correctly errors "unused argument". Same for argument collisions like `M = 5000`. This also makes the two modes of `powerbf01seq`/`powertbf01seq` inconsistent (fixed-`n` mode propagates errors; power mode returns `NaN`). The real message is recoverable only via `details = TRUE` (`$error`). Relatedly, when a *transient* evaluation failure occurs mid-range during bracketing, `.bfseq_search_before_invalid` only bisects below the failed point and never explores above it, then blames `nrange[2]` (which was never evaluated).
- **Fix:** Rethrow errors raised at the first evaluated candidate (structural errors fail every n; numerical failures usually don't); append the captured `$error` message to the NaN warning; after the invalid-boundary bisection finds nothing, resume doubling past the invalid candidate, treating invalid points as failures (the binary search already does).

### M3. `nbf01seq(type = "moment")` without `dpm` silently returns an empty `list()`
- **File:** `package/R/nbf01seq.R:215-226` — **NEW-IN-PR** — **confirmed (re-verified)**
- **Problem:** The inner function correctly stops with "argument 'dpm' must be specified when type = \"moment\"", but the exported wrapper always includes `"dpm"` in `vectorizeArgs`; when `dpm` is NULL, `mapply` receives a zero-length argument and returns `list()` with no error or warning. `powerbf01seq` has its own guard; `ntbf01seq` is immune (`dpm` defaults to `plocation`).
- **Fix:** Validate `dpm` in the wrapper before vectorizing (mirror the inner check), or exclude NULL `dpm`/`dpsd` from `vectorizeArgs` as already done for `pm`. Same pattern check for `pm` missing with type "normal", which currently fails with the cryptic `length(pm) == 1 is not TRUE` (see m12).

### M4. `powerbf01seq`/`powertbf01seq`: fixed-`n` mode does not reproduce the design found in power mode
- **Files:** `package/R/powerbf01seq.R:126-128`, `package/R/powertbf01seq.R:102-105, 112` — **NEW-IN-PR** — **confirmed**
- **Problem (two parts):**
  1. In fixed-`n` mode the schedule is rebuilt with `nrange = c(2, max(2, n))`, so the default `minN` silently becomes 2 instead of the user's `nrange[1]`. Round-trip reproduced: power mode with `by = 20, nrange = c(40, 400)` finds n = 60 with looks (40, 60), power 0.699; fixed-`n` mode with identical settings and `n = 60` builds looks (2, 22, 42, 60), power 0.634 — different design, different power, no warning.
  2. Fixed-`n` mode validates `n2` as strictly increasing, while the search evaluator (and `ptbf01seq`) accept non-decreasing `n2` (valid when `ratio < 1`); a design found in power mode can be un-reproducible in fixed-`n` mode (errors "schedule is not strictly increasing").
- **Fix:** Default `minN` from the user-supplied `nrange[1]` (clamped to ≤ n) in fixed-`n` mode, or require explicit `minN` when `by` is combined with `n`. Relax the fixed-mode `n2` check to non-decreasing (keep strict on `n1`).

### M5. Integer-coercion overflow in `.bfseq_minimum_max_n` makes the verification loop spin forever (process hang)
- **File:** `package/R/seqsearch.R:100-122` — **NEW-IN-PR** — **confirmed**
- **Problem:** For extreme `timing[1]` (< ~5e-10, also reachable via large `lookMinN` from tiny `ratio`), `lower` exceeds `.Machine$integer.max`; `as.integer(lower)` yields NA; inside `while (TRUE)`, the schedule call errors (caught by `try`) and `lower <- NA + 1L` stays NA — infinite loop. The swallowing `try()` even defeats `setTimeLimit`, so a scripted session cannot recover. A hang is the worst failure mode for pathological input.
- **Fix:** After computing `lower`, validate `is.finite(lower) && lower <= .Machine$integer.max` with a clear error ("first information fraction too small relative to the look minimum"), and bound the loop by `nrange[2]`.

### M6. Misleading "pass a wider 'trange'" warning when BF01 = k is mathematically unattainable (one-sided t)
- **File:** `package/R/seqhelpers.R:744-800, 979-991` — **NEW-IN-PR** — **confirmed**
- **Problem:** For one-sided t-tests, BF01(t) approaches a *finite supremum* in the wrong tail (e.g. `tcrit(k = 20, n1 = 15, n2 = 15, alternative = "greater")`: BF01 plateaus at ≈16.8, so BF01 = 20 has no root at any t). The code returns NaN (correct) but advises widening `trange` — advice that can never succeed. The two-sided branch has an "impossible" check; the one-sided branch has none. Note: combined with C2, this NaN currently propagates into silently-zeroed probabilities.
- **Fix:** In the one-sided branch, when the limit probe finds `f` finite, same-signed, and nearly flat, warn "BF01 = k appears unattainable for any t (supremum ≈ …); no critical value exists" instead.

### M7. `nbf01` falsely reports achievable targets as unachievable when the power curve overshoots its limit (skeptical design prior)
- **File:** `package/R/nbf01.R:64-88`; numerical fallback also affected via `helpers.R:33-35` — **PRE-EXISTING** — **confirmed**
- **Problem:** For `psd = 0`, the analytical branch compares the target against the n→∞ *limit* power and returns NaN if target > powlim. But the power curve is not monotone: when `(null+pm)/2 − dpm > 2L·dpsd²/usd²` (L = usd²·log(k)/(null−pm)) the curve overshoots its limit at finite n. Reproduced: `nbf01(k = 1/10, power = 0.05, usd = 1, null = 0, pm = 1, psd = 0, dpm = 0.3, dpsd = 0.1)` warns "power (0.05) higher than limiting power (0.02)" and returns NaN, yet `pbf01` reaches 0.101 at n = 15. `analytical = FALSE` fails identically (endpoint-only bracketing in `searchN`).
- **Fix:** The interior extremum is available in closed form (t* = (s·usd² − 2L·dpsd²)/(L·usd²), t = 1/n); compute attained max power and solve the quadratic for the smaller root when the target is attainable in a finite window (with a documented caveat), or at minimum reword the warning to "power may only be attainable in a finite sample-size window".

### M8. `pbinbf01` Brent *minimization* of a concave function can pick the wrong boundary → power exactly 0 (or 1) when the truth is large
- **File:** `package/R/pbinbf01.R:157-167` — **PRE-EXISTING** — **confirmed**
- **Problem:** The point-type log BF01(x) is strictly concave in x, so its minimum over [0, n] is at an endpoint; `optim(method = "Brent")` assumes a unimodal *minimum* and can converge to the wrong endpoint, triggering a false early return of 0/1. Reproduced: `pbinbf01(k = exp(-150), n = 324, p0 = 0.431, type = "point", a = 0.348, b = 47.2, dp = 0.02)` returns **0**; exhaustive enumeration gives **0.5289**. Mitigating: an 8,000-draw scan found no misfire at conventional k (1/10, 1/3, 3, 10); triggering needs extreme k plus a null-conflicting analysis prior.
- **Fix:** One line: `xminval <- min(logbf(0), logbf(n))` — exact for both types (concave or monotone). The Brent *maximization* nearby is sound.

### M9. Test/process: all numerical-correctness tests are gated off by default, no CI exists to ever run them, and the `pbf01seq` oracle has zero active tests
- **Files:** `package/inst/tinytest/helper-extended-tests.R:4-14` + 9 gated files; `package/inst/tinytest/test-pbf01seq.R` (100% commented out); no `.github/` workflows — **gating NEW-IN-PR** — **confirmed**
- **Problem:** Default run = 154 API/smoke assertions; the 578-assertion extended suite (paper values, the 256-case `nbf01` grid, t-boundary regressions, the only external `BayesFactor` cross-check, sequential-region paper fixtures) requires `BFPWR_RUN_EXTENDED_TESTS=true`, which neither CRAN, `package/tests/tinytest.R`, nor any CI sets. Worse, the new `nbf01seq` tests use `pbf01seq` as their oracle, and `pbf01seq` itself has no active test in the default suite — the default suite validates the search against an engine it never checks. The extended suite takes only ~43 s on this machine, so the gate is far broader than necessary. (Both suites currently pass: 154/154 and 578/578.)
- **Fix:** Gate only the wall-clock performance file and the `BayesFactor` dependency file (or use tinytest's standard `at_home()`); let correctness/paper-value tests run by default; add a CI workflow setting the env var; add fast active assertions to `test-pbf01seq.R` (one-look equivalence to `pbf01`, tail complementarity, one pinned paper value).

---

## 3. Minor

### Sequential layer (new in PR unless noted)

- **m1. Single-look `search = "exhaustive"` silently downgraded to "none" yet reports `firstCrossingCertified = TRUE`** (`seqsearch.R:891-896, 497`). The certification then rests on an unproven at-most-one-upcrossing assumption. Honor the user's exhaustive request regardless of look count, or report `firstCrossingCertified = NA` for one-look schedules.
- **m2. t-evaluator re-emits the trange search-limit warning on every candidate-n evaluation reusing a cached boundary** (`seqsearch.R:793-807, 842-852`) — warning spam with inflated counts during a single search. Warn once per distinct boundary.
- **m3. Increase-mode candidates include the off-grid `nrange[1]`**, so the returned "scheduled maximum" may not be of the documented form `minN + j*by`; and `nextend` certification steps in grid units for `by` schedules but unit steps elsewhere, while the shared doc implies unit steps (`seqsearch.R:358-365, 411-445`). Snap `lowerN` to the grid or document; clarify `nextend` units.
- **m4. One-sided boundary "certification" uses `rel.tol = 1e-2` integration while the two-sided path uses full precision** (`seqhelpers.R:871-884, 962-966`; same pattern in `ptbf01.R:79-107, 199-204`). Measured boundary error up to ~2.4e-5 in t — harmless in practice, but "certify" overstates it and precision is silently asymmetric between alternatives. Certify the accepted root with one extra full-precision evaluation.
- **m5. Hardcoded `ngrid = 1000` QMC points; `...` cannot raise integration precision** (`seqhelpers.R:97-101, 137`; docs in `pbf01seq.R:41`, `ptbf01seq.R:15` say "... passed to mvtnorm::lpmvnorm" — **PRE-EXISTING**). Passing `M = 5000` errors with "matched by multiple actual arguments"; `seed` silently does nothing. ~2e-4 per-stage integration error is fixed and not user-controllable. Expose `ngrid`/`M`, document fixed `w`/`seed`.
- **m6. No validation that looks are ordered/increasing in `pbf01seq`/`ptbf01seq`** (**PRE-EXISTING**; `pbf01seq.R:86-98`, `ptbf01seq.R:71-80`). Decreasing `se`/`n` orderings run silently and return meaningless stage probabilities; duplicated looks crash deep in `chol` with an opaque "leading minor" error. The new search layer validates its own schedules, but the direct entry points don't. Add `all(diff(se) < 0)` / `all(diff(n1) > 0)` checks.
- **m7. `ptbf01seq`: `all(n1 != n2)` should be `any(...)`** (**PRE-EXISTING**; `ptbf01seq.R:111-117`) — the warn-and-reset of `n2` for non-two-sample types fires only if *every* element differs; partial overlap leaves a wrong `EN2` in the returned object.
- **m8. `method = "pmvnorm"` errors on empty (crossed-boundary) regions** that the default `lpmvnorm` path evaluates as ≈0 (**PRE-EXISTING**; `seqhelpers.R:176-184`). Short-circuit `region[1,] >= region[2,]` → 0.
- **m9. Silent `na.rm = TRUE` over region probabilities** hides genuinely failed backend integrations (**PRE-EXISTING**; `seqhelpers.R:188`). Distinguish intentional NaN-boundary encoding from backend NA and warn on the latter.
- **m10. Known-variance normal approximation for sequential/fixed t power is undocumented** (**PRE-EXISTING**; `ptbf01.R:82-83`, `ptbf01seq` via `predpars`). Measured error ~1–2 percentage points per look below n ≈ 30/group, vanishing by n ≈ 50 (consistent between fixed and sequential, so internal consistency holds). Document the magnitude; longer-term, when `dpsd = 0` the exact noncentral-t computation is free for fixed designs (`pt(tcrit, df, ncp)`).
- **m11. `integer = FALSE` is a no-op in `nbf01seq`/`ntbf01seq`** but docs (inherited from `nbf01`) imply a fractional solution; one-look equivalence with fixed functions breaks under `integer = FALSE` (`nbf01seq.R:99-102`, `ntbf01seq.R:106-109`). Document or drop.
- **m12. Wrapper UX gaps:** no partial matching for `type`/`alternative`/`target` in seq wrappers, diverging from `match.arg` behavior of fixed counterparts (`seqsearch.R:942-954`); non-scalar `details` silently treated as FALSE (`nbf01seq.R:201`); missing `pm` gives cryptic `length(pm) == 1 is not TRUE` (`nbf01seq.R:70-75`); `minN` silently ignored without `by`, `timing` silently overrides `looks` (`seqsearch.R:37-49`); non-integer `by` silently `ceiling`-ed; failure branches report `nextend = 0` regardless of user setting; `maxEvaluations` meaningless for adaptive timing searches.

### Fixed-design layer (all PRE-EXISTING)

- **m13. `ntbf01` forwards `...` to both `ptbf01` and `uniroot`** (`ntbf01.R:68-77` + `helpers.R:20-50`): passing the documented `ptbf01` argument `drange` makes every `uniroot` call error → NaN with a misleading warning. Split the dots.
- **m14. `pnmbf01` Lambert-W argument overflow** (`pnmbf01.R:46`): for k ≲ 1e-295 or absurd n, `lambertW0(Inf) = Inf` → power 0 instead of ≈1. Compute the argument in log scale and use the asymptotic expansion `W₀(e^la) ≈ la − log(la) + log(la)/la` for large `la`. (`pbf01` is already fully log-scale and immune.)
- **m15. `pbf01` degenerate prior `psd = 0, pm = null`, `k = 1` → NaN** (`pbf01.R:52`); should be 1 (BF01 ≡ 1). Early-return `as.numeric(k >= 1)` or error informatively.
- **m16. `searchNoscil` crashes on unreachable targets** (`helpers.R:84`, affects `nbinbf01`): `searchN` returns NaN by design, then `seq(NaN, ...)` errors with "'from' must be a finite number" instead of returning NaN. Add `if (is.nan(ni)) return(NaN)`.
- **m17. `pbinbf01` uniroot default tolerance (~1.2e-4) can misclassify the boundary integer** when k is within ~1e-4 relative of an attained BF value (`pbinbf01.R:175-178, 209`); worst observed error 0.0057 in fuzzing. Snap the integer boundary by direct `logbf` evaluation, or pass `tol = 1e-9`.
- **m18. `powerbf01`/`powertbf01` accept both `n` and `power` simultaneously** (silently ignoring `power`) while their error message claims exactly one must be NULL (`powerbf01.R:65`, `powertbf01.R:41`). The new seq wrappers get this right (`is.null(n) == is.null(power)`) — port back.
- **m19. Efficiency:** duplicate evaluation of the certify integral at the origin (`ptbf01.R:199` + `seqhelpers.R:640` — one wasted numerical integration per one-sided boundary search, NEW-IN-PR); exact wrong-tail path forced for all one-sided |t| ≥ 4 even where the fast path agrees to 1e-6 (~300× slowdown per call, `tbf01.R:161-175` — correctness-first is defensible; a cheap consistency check would recover most of the cost).

---

## 4. Nits / documentation

- **n1.** `pbf01seq.R:12-14` `@param k0`: "implies evidence for H1" → should be **H0**. Propagates into new man pages via `@inheritParams`.
- **n2.** `pbf01seq.R:20-21` `@param psd`: 'standard deviation (`type = "moment"` and `type = "directional"`)' — first "moment" should be **"normal"**. Propagates into `nbf01seq.Rd:51-52`.
- **n3.** `plot.bfseqdesign` `@param x` says class `"power.bftest"` → should be `"bfseqdesign"` (`pbf01seq.R:451`); orphaned roxygen block at `pbf01seq.R:659-670`; `nullplot` validated twice / `zplot` never in the `stopifnot` (`pbf01seq.R:495-507`); "Logcal" typos.
- **n4.** `powerbf01.R:433-437`: H0-panel threshold label wrong for k > 1 (prints "1/0.333" instead of "1/3").
- **n5.** `nbf01.R:76-78`: warning rounds both powers to 2 digits, producing self-contradictory "power (0.8) higher than limiting power (0.8)". Use `signif(powlim, 4)`.
- **n6.** `NEWS.md` has no entries for the four new exports nor for the **breaking `drange` → `trange` rename** in `ptbf01seq`.
- **n7.** `package/man/ntbf01.Rd` is stale vs roxygen source (regeneration churns); `DESCRIPTION` lacks `RoxygenNote` and carries a nonstandard `Config/roxygen2/version: 8.0.0`.
- **n8.** `package/out/` tracks ~12 built tarballs/PDFs (PR adds two more; repo pack already ~105 MB). `.Rbuildignore` excludes it anyway — stop tracking, use releases.
- **n9.** Stray untracked artifacts at repo root: `.review-tinytest-results.rds` (delete/ignore before merge).
- **n10.** `numerical-helpers.R`: `.bfpwr_lpnorm_interval` doesn't validate `lower <= upper` (NaN/-Inf + base-R warning if crossed); `.bfpwr_logspace_sum` silently drops NaN (not just −Inf) terms — keep NaN propagating so future bugs surface. `pbinbf01.R:231`'s `min(0, logpow)` clamp would likewise hide a future double-count.
- **n11.** `ptbf01.R:317-320`: `"type"` listed twice in `Vectorize(vectorize.args = ...)`. Harmless.
- **n12.** ~1e-4-level QMC integration error means the selected n can differ by ±1 from a higher-precision evaluation; worth one docs sentence (ties into m5).
- **n13.** Document the monotone-BF01-in-t assumption at `.bfpwr_one_sided_adaptive_root` / `genregions1` (valid for all shipped priors; future priors could break it). The two-sided "impossible" check optimizes only over the z-guess interval — fine in practice, optionally probe ±2 prior SDs around `plocation`.
- **n14.** Test-suite polish: conditional assertions that silently skip when boundary cases occur (`test-nbf01seq.R:102,150`); golden-value pins brittle to integrator changes (each is at least paired with independent invariant checks); missing coverage for `null != 0` through the seq wrappers (the centering transforms are exactly the wrapper-owned logic — manually verified correct, but should be pinned), `nbf01seq(type = "moment")` via the exported wrapper (would have caught M3), `integer = FALSE`, vectorized `type`/`alternative`, power < 0.5 and `lower.tail = FALSE` in the `nbf01` analytical branch (would have caught M7), one-sided/informed-prior `BayesFactor` cross-checks, two-sided informed-prior k > 1 (would have caught C1). The active `test-pbinbf01.R` complementarity check is tautological (computed from the same `logpow`); the genuinely validating simulation blocks are commented out.

---

## 5. Verified correct (what was checked and passed)

For calibration, the following were re-derived and/or numerically validated and found **correct**:

- **Sequential covariance/predictive** (`predpars`): mean `dpm/se_i`, cov `sqrt(inf_min/inf_max) + dpsd²·sqrt(inf_i·inf_j)` — matches the hand-derived joint distribution; multi-look designs match 50k-replicate Monte Carlo within MC error (probabilities, expected N, sd(N)).
- **First-passage decomposition** (`genregions1/2`, strict and non-strict): no double counting; region counts match exhaustive enumeration; H0-interval nesting and the `strict = FALSE` collapse logic verified by case analysis.
- **1-look reductions:** `pbf01seq` = `pbf01`/`pnmbf01` to machine precision; `ptbf01seq` = `ptbf01` to ≤1.6e-5; z vs t agreement at large df.
- **BF↔boundary inversions** (`zcrit`, `tcrit`): residuals ~1e-15 (z) and ≤5e-6 (t), including the PR's new log-scale directional path; the new two-sided "impossible" check via bounded `optimize` is sound.
- **Closed-form z-test core:** BF01 and normal-moment BF re-derived symbolically and match exactly; the Lambert-W branch choice in `pnmbf01` (`lambertW0`) is provably the only real branch; `pbf01`'s bounded-above handling (`X < 0`/`Y < 0` → P = 1) is algebraically correct; 2×10⁵-sim MC agreement across 10 configurations; the PR's `2*log(k)` and log-space tail rework is identical to the old arithmetic to 2e-16 while fixing k² underflow and extreme-tail cancellation.
- **t-test BF (`tbf01`):** matches `BayesFactor::ttest.tstat` to 7.6e-7 over a 48-case two-sided grid; the PR's new exact wrong-tail scale-mixture path was confirmed by independent Monte Carlo of the exact identity — in the extreme wrong tail it is *more* correct than `BayesFactor` (which returns impossible values there). Monotonicity of one-sided BF01 in t (the root-search assumption) was proven via MLR and verified numerically. df/neff conventions and one-sided prior renormalization correct. Fast/exact paths agree to 1e-6 in their overlap with no discontinuity at the |t| = 4 cutoff.
- **`ptbf01` MC-validated in 6 configurations** (one-sided both directions, informed and default priors, dpsd ∈ {0, >0}, two-sided k < 1 and symmetric k > 1) — all within MC error except the C1 case.
- **Binomial family:** point-type log BF01 strict concavity and direction-type monotonicity proven (so the two-root/one-root region architecture is sound — the feared non-contiguous region case cannot occur); `pbinbf01` matches brute-force enumeration to ~1e-15 over a 5,760-case grid + 3,000-case fuzz (outside M8/m17 regimes); the PR's log-scale predictive rework (truncated-beta normalization, O(1) tail-stable point-design path — also an efficiency win) confirmed exact; `nbinbf01` paper values (Kelter & Pawel) pass.
- **`dirbf01` rework:** adjudicated against 256-bit arithmetic — new code accurate to 3e-13 in regimes where the old formulation had 14–197% error or NaN.
- **Search engine invariants:** binary-search exit certifies n succeeds and n−1 fails (cached); backward walk and adaptive probe provably terminate; `nextend` certification fail-closed semantics correct; prefix caching keys are correct; power evaluation is fully deterministic (fixed-seed QMC), so caching and exact comparisons are sound; progress-callback contract matches docs.
- **Wrappers:** analysis vs design prior routing verified end-to-end (no `pm`/`dpm` mix-ups); every wrapper argument reaches the engine; one-look seq results match fixed-design counterparts for H1/H0 targets, one/two-sample, and `null ≠ 0`; NAMESPACE/S3 registrations complete; all new man-page examples run (≤0.11 s each).
- **Test suite:** 154/154 default and 578/578 extended assertions pass; the adversarial files are genuinely adversarial (not tautological); no unseeded-RNG dependence; no R CMD check style issues (no T/F, no browser(), clean `checkUsage`, all imports genuinely used).

---

## 6. Suggested course of action (for discussion)

1. **Before merge (correctness):** C2 (both halves), M1, M2, M3, M4. These are silent-wrong-answer or contract-violation bugs in code paths users of this PR will hit. C1 was fixed in the 2026-06-14 follow-up.
2. **Before merge (cheap & high value):** M5 (overflow guard), M6 (warning text), M8 (one-line fix), m7 (`any` vs `all`), n1/n2 (doc typos newly propagated), n6 (NEWS for new exports + breaking rename).
3. **This PR or fast-follow:** M7 (closed-form fix exists), M9 (un-gate correctness tests + minimal CI + activate `test-pbf01seq.R`), m4–m6, m13–m18, remaining docs/nits.
4. **Backlog:** m10 (exact noncentral-t power for `dpsd = 0`), m19 (exact-path speed), n8 (repo hygiene for `package/out/`).
