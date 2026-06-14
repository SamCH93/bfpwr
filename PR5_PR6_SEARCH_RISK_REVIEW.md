# PR5/PR6 Search Risk Review

This note records undesirable defensive/fallback behavior found while reviewing
the PR5 and PR6 search and numerical-stability changes.

The common pattern is that a change intended to make a numerical search more
robust can also hide a mathematical invariant, turn an unresolved boundary into
a plausible numeric answer, or suppress the diagnostic needed to understand the
failure.

## 2026-06-13/14 Implementation Status For PR6

Implemented for the one-sided t search paths:

- Fixed-n and sequential one-sided adaptive t searches no longer use a terminal
  `|t| <= 256` cap. They derive the finite search limit from the predictive
  distribution and `tail.eps` (default `1e-3`).
- `tail.eps` is used as an omitted-tail bound. In fixed-n power calculations,
  a cutoff without a bracket returns the boundary-free 0/1 approximation with a
  warning that reports the bound. In sequential calculations, an unresolved
  boundary has omitted marginal mass bounded by `tail.eps`.
- Sequential `tail.eps` is not split across looks or boundaries. The same
  user-facing value is used for each adaptive boundary search.
- One-sided searches now scan only the mathematically expected direction
  (`try_opposite = FALSE`). This removes the superficial opposite-side retry
  while preserving valid wrong-tail H0 roots when the sign at the origin points
  there.
- Sequential continuation handling was tightened: an unresolved earlier H1
  boundary is converted to unbounded continuation instead of allowing `NaN` to
  zero all later paths; a missing final stopping boundary remains an empty
  stopping event.
- Sequential t boundary searches now carry an internal status sidecar before
  region generation. `ok`, `impossible`, and adaptive `tail_cutoff` boundaries
  are allowed; `search_failed` and `error` statuses stop before integration
  instead of becoming empty regions.
- Sequential integration still treats generated-region `NaN` as an empty
  stopping event, but non-empty `NA` bounds and backend `NA`/`NaN`
  probabilities now error instead of being silently dropped or failing
  indirectly.
- Sequential sample-size searches now preserve boundary-evaluation failures in
  solver diagnostics, and `ntbf01seq(details = FALSE)` / `powertbf01seq()` now
  surface those messages instead of returning only an unexplained `NaN` or a
  generic failure.
- Follow-up audit fixes on 2026-06-14 tightened the status distinctions:
  adaptive one-sided `tail_cutoff` is reported only after the stable
  certifying function is finite at the predictive cutoff; a non-finite endpoint
  remains a `search_failed` boundary.
- For two-sided numeric `trange`, a "maximum BF below threshold" result is no
  longer blindly accepted as an impossible H0 boundary. The wrapper performs an
  adaptive proof pass: if adaptive search also proves the H0 boundary
  impossible, the empty H0 region is allowed; if adaptive search finds roots,
  the numeric interval is treated as a range miss and errors before region
  generation.
- Cached `ntbf01seq()` / `powertbf01seq()` t-boundary failures now report the
  actual failed look index, and structural boundary failures no longer emit the
  older generic `Power = NaN` warning before the real solver diagnostic.
- Region generation now rejects plain `NA` critical values up front. `NaN`
  remains the intentional internal sentinel for empty generated regions.
- One-sided adaptive t root search now separates the fast scout, relaxed
  search/bracketing function, and full final certification function. Returned
  roots and predictive-tail cutoff endpoints must be certified by the full BF
  path, while adaptive probes still use the cheaper search path.
- `.bfpwr_logspace_sum()` now drops only `-Inf` zero-mass terms. `NA`, `NaN`,
  and positive `Inf` now stop as invalid numerical state instead of being
  silently removed from the sum.
- Benchmark check against the previous PR6 head (`e670fcf`) found no material
  unit-test slowdown (51.7 s vs 51.2 s for the full extended suite). Regular
  adaptive outputs matched; the extreme sequential case stayed within
  `tail.eps = 0.001` relative to explicit `trange = c(-50, 50)` while remaining
  hundreds of times faster.

## 1. One-Sided `tcrit()` Opposite-Side Retry

Status after 2026-06-13 change: **Fixed.** One-sided adaptive searches now use
`try_opposite = FALSE` and scan only the mathematically expected direction.
Valid wrong-tail H0 roots remain supported when the sign at the origin points
there.

Locations:

- `package/R/seqhelpers.R`: `.bfpwr_one_sided_adaptive_root()`
- `package/R/seqhelpers.R`: one-sided `tcrit()` now calls the helper with
  `try_opposite = FALSE`
- `package/R/ptbf01.R`: fixed-design `ptbf01()` calls the same helper with
  `try_opposite = FALSE`

Issue:

PR5 replaced `uniroot(..., extendInt = ...)` with a finite manual scan up to
`|t| <= 256`. As part of that change, one-sided `tcrit()` began searching the
side implied by `BF01(0) - k` and then retrying the opposite side if the first
side did not bracket a root. PR6 preserved this behavior through the shared
helper.

For the current one-sided t Bayes factor, the threshold side is determined by
the sign at the origin and the alternative direction. For `alternative =
"greater"`, BF01 decreases as t moves in the greater/evidence-for-H1 direction;
for `alternative = "less"`, the mirror statement holds. The side selected from
`BF01(0) - k` is therefore the only mathematically valid side for the threshold.

Important distinction:

"Wrong-tail" H0 roots are legitimate and should remain supported. For example,
with `alternative = "less"` and `k > 1`, the H0-evidence boundary can be at
positive t. That is not the issue. The issue is retrying the opposite side after
the side implied by `BF01(0) - k` fails.

Suggested change:

- Change sequential one-sided `tcrit()` to pass `try_opposite = FALSE`.
- Keep the helper argument only if we need generic use elsewhere.
- Add a short comment documenting the one-sided monotonicity assumption and the
  distinction between valid wrong-tail H0 roots and invalid opposite-side
  fallback.
- Add a regression test that verifies no opposite-side retry is needed for H1
  and H0 one-sided boundaries.

Decision:
- agree

## 2. `ptbf01()` Search-Limit Fallback Returns Boundary-Free Probability

Status after 2026-06-13 change: **Fixed by replacing the fixed search limit
with a probability-bound cutoff.** The function still returns the
boundary-free 0/1 approximation for compatibility, but it now does so only
after the searched predictive tail is bounded by `tail.eps`, and the warning
reports that omitted-tail bound. A numeric `drange` remains the exact-search
escape hatch.

Location:

- `package/R/ptbf01.R`: one-sided adaptive branch after
  `bfpwr_ptbf01_search_limit`

Issue:

When the adaptive one-sided power-boundary search reaches `|t| <= 256` without
bracketing a root, `ptbf01()` warns but returns the boundary-free probability
implied by the sign at the null. This is correct for a truly impossible
threshold, but it also covers the case where a root exists beyond the finite
search limit.

That means a search limitation can become a valid-looking probability, usually
0 or 1, instead of an unresolved boundary result.

Suggested change:

- Separate "threshold impossible" from "finite adaptive scan did not find a
  root."
- Return `NaN` for the unresolved finite-scan case unless we can prove the
  threshold is impossible.
- If keeping the current behavior for compatibility, expose the distinction in
  a structured diagnostic field or a more specific warning.
- Add tests for a forced small `search_limit` or equivalent helper-level case
  where a root exists beyond the limit.

## 3. Sequential T Boundary Warning Suppression Is Too Broad

Status after 2026-06-13 change: **Fixed by boundary-status classification.**
The sequential t evaluators now collect `tcrit()` warnings into an internal
status sidecar. Adaptive predictive-tail cutoffs and impossible H0 boundaries
are aggregated and re-emitted; numerical search failures stop before region
generation.

Locations:

- `package/R/ptbf01seq.R`: local `evalTcrit()`
- `package/R/seqsearch.R`: `.bfseq_t_schedule_evaluator()` local `evalTcrit()`

Issue:

Both sequential t evaluators call `tcrit()` under `withCallingHandlers()` and
muffle every warning. They only count warnings whose message contains
`"Adaptive t critical-value search reached"`.

As a result, unrelated warnings such as "Numerical problems finding critical
value" can disappear. The downstream code may then treat the resulting `NaN`
boundary as an empty stopping region, producing a plausible sequential design
instead of surfacing the boundary failure.

Suggested change:

- Only muffle the adaptive-limit warning after recording it.
- Let all other `tcrit()` warnings propagate.
- Alternatively collect non-limit warnings and emit one aggregate warning that
  includes their messages and counts.
- Add a test that injects or triggers a non-limit `tcrit()` warning and verifies
  that `ptbf01seq()` or `ntbf01seq()` does not silently swallow it.

Decision:
- Muffle all search warnings, count them, only ommit the in the aggregate form in the end (the function is otherwise way too verbose)

## 4. `NaN` Boundaries Are Used For Both Empty Regions And Failed Searches

Status after 2026-06-14 follow-up: **Fixed for sequential t boundary failures,
adaptive tail cutoffs, numeric `trange` range misses, and integration `NA`.**
The sequential t paths now classify boundary results before region generation
and reject failed searches. Generated-region `NaN` remains the internal
empty-stopping-event sentinel, but failed t-boundary searches no longer reach
integration as bare `NaN`. Non-empty `NA` bounds and backend `NA`/`NaN` values
now error instead of being silently dropped or failing indirectly.

Remaining caveat: this is not a public boundary-status object. The sidecar is
internal to the sequential t calculation, and z-boundary generation still uses
the established numeric sentinel convention.

Locations:

- `package/R/seqhelpers.R`: `.bfseq_intstage_sum()`
- `package/R/seqhelpers.R`: region probability summation with `na.rm = TRUE`

Issue:

Sequential integration treats any region containing `NaN` as probability zero.
This behavior existed before PR5 and is correct when `NaN` encodes a boundary
that mathematically does not exist, such as an impossible H0 stopping boundary.

The risk is that PR5/PR6 introduced more boundary-search failure paths that also
return `NaN`. Those are not empty regions. They are unresolved numerical
failures. The current integration layer cannot distinguish the two meanings.

Suggested change:

- Stop using bare `NaN` for both "impossible boundary" and "search failed."
- Introduce an internal boundary object or sidecar status with values such as
  `ok`, `impossible`, `search_limit`, and `error`.
- Keep impossible boundaries as empty regions.
- Propagate search failures to the caller instead of converting them to zero
  probability.
- Remove or narrow `na.rm = TRUE` so backend integration `NA` values are not
  silently dropped.

Implemented change:

- Added an internal `tcrit()` result sidecar with statuses `ok`, `impossible`,
  `tail_cutoff`, `search_failed`, and `error`.
- Direct `ptbf01seq()` and the `ntbf01seq()`/`powertbf01seq()` schedule
  evaluator validate those statuses before calling `genregions1()` or
  `genregions2()`.
- H0 `impossible` and adaptive `tail_cutoff` statuses are allowed and reported;
  `search_failed`/`error` statuses stop with a boundary-specific message.
- `.bfseq_intstage_sum()` still converts generated-region `NaN` to probability
  zero, but it now errors if a non-empty region contains `NA` bounds or if
  integration returns `NA`/`NaN`.
- Sample-size search no longer hides an invalid boundary candidate behind the
  last finite below-target candidate when no crossing has been found.
- One-sided adaptive tail cutoffs now require a finite certifying endpoint at
  the predictive cutoff. Non-finite endpoint evaluations are `search_failed`,
  not `tail_cutoff`.
- Numeric two-sided `trange` H0 misses are distinguished from proven impossible
  H0 boundaries by an adaptive proof pass.
- `genregions1()` and `genregions2()` reject plain `NA` critical values; only
  `NaN` is accepted as the empty-boundary sentinel.

## 5. One-Sided Root "Certification" Uses Relaxed Integration Tolerance

Status after 2026-06-14 Batch 1 change: **Fixed.** One-sided adaptive t search
now uses a three-function split: the fast direct-integral scout only proposes
candidate brackets, the relaxed stable path performs adaptive search/bracketing,
and the full BF path certifies any returned root or predictive-tail cutoff
endpoint.

Locations:

- `package/R/seqhelpers.R`: `tcrit()` builds `rootFunSearch` with
  `rel.tol.default = 1e-2`
- `package/R/seqhelpers.R`: one-sided `tcrit()` passes
  `certify_fun = rootFunSearch`
- `package/R/ptbf01.R`: same pattern for fixed-design power boundaries

Issue:

PR6 improved the one-sided search by using the fast direct integral only as a
scout and requiring candidate roots to be certified. However, the certification
function is still the search function with relaxed integration tolerance
(`rel.tol = 1e-2` by default), not the full stable BF path.

That is better than trusting the scout, but it makes "certified" slightly
misleading. A root can be certified against the relaxed search tolerance rather
than the most accurate available calculation.

Suggested change:

- Use the full `rootFun` as `certify_fun` for the final root certification.
- Keep `rootFunSearch` and `rootFunFast` only for bracketing/scouting.
- If the full path is too slow for every certification, certify only the final
  root and perhaps the final bracket endpoints.
- Rename comments if relaxed certification is intentionally retained.

Decision: change so the certification happens with the appropriate function

Implemented change:

- `.bfpwr_one_sided_adaptive_root()` now accepts `certify_fun`,
  `search_fun`, and `scout_fun`.
- `ptbf01()` and `tcrit()` pass the full `rootFun` as `certify_fun`, the
  relaxed `rootFunSearch` as `search_fun`, and the fast direct integral as
  `scout_fun`.
- `.bfpwr_certified_root()` no longer accepts a search-function endpoint root
  unless the full final function also certifies that endpoint. If the relaxed
  root is not certified, it attempts a full-function root over the bracket.
- Predictive-tail cutoff status now requires the full final function to be
  finite and same-sign at the cutoff endpoint.
- Added helper-level tests that force `search_fun` and `certify_fun` to have
  different roots and that force a non-finite full endpoint at the cutoff.

## 6. Log-Space Summation Drops All Non-Finite Values

Status after 2026-06-14 Batch 1 change: **Fixed.** `.bfpwr_logspace_sum()` now
removes only `-Inf`, which represents zero probability/mass. `NA`, `NaN`, and
positive `Inf` now error instead of being silently dropped.

Location:

- `package/R/numerical-helpers.R`: `.bfpwr_logspace_sum()`

Issue:

`.bfpwr_logspace_sum()` filters with `is.finite(logx)`. This correctly removes
`-Inf` terms that represent zero probability, but it also removes `NaN`, `Inf`,
and `NA`. That can hide upstream numerical bugs and turn a partially invalid
sum into a valid-looking finite value.

Suggested change:

- Only drop `-Inf`.
- Treat `NaN`, `NA`, and positive `Inf` as errors or return `NaN`.
- Audit callers that rely on current behavior and handle zero-probability terms
  explicitly.

Decision: agree with the suggestion

Implemented change:

- `.bfpwr_logspace_sum()` rejects `NA`, `NaN`, and `Inf` with an internal
  numerical-state error.
- Existing legitimate zero-probability behavior is preserved: finite terms plus
  `-Inf` sum correctly, and all `-Inf` returns `-Inf`.
- Added direct helper tests for finite-plus-`-Inf`, all-`-Inf`, `NA`, `NaN`,
  and positive `Inf`.

## 7. Sequential Sample-Size Search Converts Evaluation Errors To `NA`

Status after 2026-06-14 Batch 2 change: **Fixed.** `.bfseq_search()` now catches
only the typed internal `bfseq_candidate_invalid` condition. Ordinary evaluator
errors, malformed evaluator returns, invalid arguments, and programming bugs now
propagate as structural errors. Expected candidate-level numerical invalidity
is recorded as search state with `status`, `reason`, `terminal`, and `error`
diagnostics. Sequential t boundary failures and non-finite stage probabilities
use this typed path, and `details = FALSE` wrappers surface the stored
diagnostic instead of silently returning `NaN`.

Location:

- `package/R/seqsearch.R`: `.bfseq_search()` local `evalN()`

Issue:

PR6 introduced the generalized sequential sample-size search. Its evaluator
wraps every design evaluation in `try(..., silent = TRUE)` and converts errors
to `criterion = NA`. Later code often reports this as an upper-bound or
non-finite-power problem.

This makes structural errors, invalid arguments, and true numerical failures
look like ordinary unreachable-power cases.

Suggested change:

- Replace broad `try(evaluate(n), silent = TRUE)` handling with a typed
  internal condition for expected candidate-level invalidity, for example
  `bfseq_candidate_invalid`.
- Catch only that typed condition inside `.bfseq_search()`. Ordinary errors
  from invalid arguments, malformed evaluator returns, or programming bugs
  should propagate.
- Let evaluators explicitly classify expected numerical invalidity by signaling
  the typed condition. Candidate records should preserve structured diagnostics
  such as `status`, `reason`, `terminal`, and the user-facing `error` message.
- Keep malformed evaluator outputs as structural errors. For example, missing
  `power`, non-numeric `power`, or missing `result` should not become an
  invalid candidate.
- Use a `terminal` flag on the typed condition so this change prepares the
  policy decision in item 8. Initially, boundary failures from explicit numeric
  `trange` and non-finite integration probabilities can remain terminal unless
  we later prove a transient-invalid case is meaningful.

Implementation sketch:

- Add a small internal helper to construct/signal the condition, e.g.
  `.bfseq_candidate_invalid(message, reason, terminal = TRUE, ...)`.
- Add an internal validator for successful evaluator output before computing
  `criterion`.
- Change `.bfseq_search()` to use `tryCatch()` with a handler only for
  `bfseq_candidate_invalid`; let all other errors escape.
- Preserve the existing public surface where practical: `ntbf01seq(details =
  TRUE)` can still expose `solver$error`, and `details = FALSE` can still warn
  once with the stored message.

Testing guidance:

- Keep tests minimal and targeted. Do not add a large matrix of artificial
  evaluator cases.
- One helper-level test should verify that an ordinary `stop("bug")` from
  `evaluate()` propagates rather than becoming `n = NaN`.
- One helper-level test should verify that a typed candidate-invalid condition
  is recorded in the solver result with the expected message/reason.
- Existing narrow-`trange` sequential t tests should cover the real package
  path; extend them only if the public diagnostic changes.

Decision:

- Adopt typed invalid candidates as the shared mechanism for items 7, 8, and
  9. Real errors should propagate; only expected candidate-level invalidity
  should enter the search state.


## 8. Search Stops After A Transient Invalid Candidate

Status after 2026-06-14 Batch 2 change: **Fixed.** Invalid candidates now carry
a `terminal` flag. Adaptive search remains conservative: terminal invalidity is
treated as a hard limit, while transient invalidity stops adaptive bracketing
with a diagnostic that recommends `search = "exhaustive"`. Exhaustive search is
the robust path for transient invalidity: it skips transient invalid candidates,
stops at terminal invalid candidates, and returns the first later finite
candidate that reaches the target.

Location:

- `package/R/seqsearch.R`: `.bfseq_find_bracket()` and
  `.bfseq_search_before_invalid()`

Issue:

When bracketing encounters an invalid candidate, the search scans between the
last valid point and the invalid point. If no solution is found there, it
returns the last finite point or invalid limit. It does not continue above the
invalid candidate.

That is conservative if invalidity means all larger sample sizes are invalid.
But if the invalidity is transient, the search can miss a valid solution at a
larger sample size.

Example failure shape:

```text
N = 20 valid, below target
N = 40 invalid
N = 80 valid, reaches target
```

The current adaptive bracketing path scans only between 20 and 40 after seeing
the invalid candidate at 40. It never checks 80, so the solver can report
`NaN`/unreached even though a valid later design exists.

Suggested change:

- Tie this to the typed invalid condition proposed in item 7. Invalid
  candidates should carry a `terminal` flag and a structured `reason`.
- Do not treat every invalid candidate as globally terminal by default. With
  timing schedules, a larger maximum sample size changes all interim looks, so
  failure at one candidate does not prove all larger candidates fail.
- Keep `search = "adaptive"` conservative. If it encounters a transient
  invalid candidate while bracketing, stop with a diagnostic that reliable
  adaptive bracketing was invalidated and suggest `search = "exhaustive"` for
  transient-invalid domains. Continue to stop immediately for
  `terminal = TRUE`.
- Make `search = "exhaustive"` the robust mode for transient invalidity. It
  should scan the candidate domain and treat transient invalid candidates as
  skipped/failed points, returning the first finite candidate that reaches the
  target. If no finite candidate reaches the target, report the most useful
  invalid-candidate diagnostic.
- This also overlaps with item 9: exhaustive search should not depend on first
  finding an adaptive bracket.

Testing guidance:

- Keep coverage minimal. Add one synthetic helper-level test where exhaustive
  search skips an invalid candidate and finds a later valid crossing.
- Add one adaptive-path assertion only if the public diagnostic changes.
- Do not add a large grid of artificial invalidity patterns; real sequential t
  boundary tests should continue to cover the package path.

Decision:

- Treat invalidity as typed state, not as a plain `NA` criterion. Adaptive
  search remains conservative: terminal invalidity stops immediately, and
  transient invalidity stops with a diagnostic that reliable bracketing was
  invalidated and `search = "exhaustive"` should be used for transient-invalid
  domains.

## 9. Exhaustive Search Only Applies After Bracketing

Status after 2026-06-14 Batch 2 change: **Fixed.** `search = "exhaustive"` now
has full-range semantics before bracketing. Timing schedules scan every
candidate maximum sample size in the normalized `nrange`; `by`/`minN` schedules
scan the scheduled candidate grid. Exhaustive mode returns the first finite
candidate that reaches the target, so an isolated early crossing no longer
depends on adaptive bracketing finding a later bracket first.

Location:

- `package/R/seqsearch.R`: `.bfseq_binary_search()` and
  `.bfseq_first_scan()`

Issue:

For multi-look timing schedules, `search = "exhaustive"` scans for the first
success only after an initial bracket has already been found. If bracketing
fails because the target function is nonmonotone or has an early isolated
success window, exhaustive mode may still miss it.

Current meaning:

```text
search = "exhaustive" exhaustively scans inside the bracket found by adaptive
bracketing.
```

It does not currently mean:

```text
scan the full sample-size candidate range.
```

Example failure shape:

```text
N = 20 below target
N = 24 reaches target
N = 30 below target
N = 40 below target
N = 80 below target
```

If adaptive bracketing probes 20, 40, and 80, it never sees the isolated
success at 24. Because no bracket is found, the later exhaustive first-crossing
scan never starts, and the solver can report unreached even though a valid
candidate exists.

Suggested change:

- Give `search = "exhaustive"` full-range semantics. It should scan the
  complete candidate domain before relying on adaptive bracketing.
- For timing schedules, scan every integer maximum sample size in the normalized
  `nrange`.
- For `by`/`minN` increase schedules, scan the scheduled candidate grid rather
  than every integer.
- Return the first finite candidate whose target probability reaches the
  requested power.
- If no finite candidate reaches the target, return unreached with the best
  available diagnostic from the upper bound and any invalid candidates.
- Keep `search = "adaptive"` as the fast bracketing mode. This gives the user a
  clear semantic split:

```text
adaptive   = fast, assumes a mostly searchable/monotone target
exhaustive = slower, checks the whole candidate domain
```

Testing guidance:

- Add one synthetic helper-level test with an early isolated success window
  that adaptive bracketing misses but full-range exhaustive search finds.
- Reuse existing public z/t sequential tests for normal monotone cases.
- Do not add a broad artificial nonmonotonicity matrix; the important contract
  is the full-range exhaustive semantics.

Decision:

- Give `search = "exhaustive"` full-range semantics. Exhaustive search should
  scan the whole candidate domain, skip typed transient-invalid candidates, and
  return the first finite candidate that reaches the target. This is the robust
  alternative when adaptive bracketing is invalidated.

## 10. `pbinbf01()` Brent Extrema Are A Pre-Existing Search Assumption

Status after 2026-06-14 review: **Documented; no planned PR6 action.** This is
outside the one-sided t search changes and predates PR5. We are not changing
this behavior in this PR6 follow-up.

Location:

- `package/R/pbinbf01.R`: `xmax <- optim(..., method = "Brent")`
- `package/R/pbinbf01.R`: `xmin <- optim(..., method = "Brent")`

Issue:

This behavior existed before PR5, so it was not introduced by PR5/PR6. However,
it belongs in the same risk class: the code uses continuous optimization over a
discrete binomial support to decide whether thresholds are possible and where
roots should be searched.

For point-null binomial BFs, relying on Brent extrema can pick the wrong
effective side or miss boundary behavior on the integer support. PR5 made the
probability summation more stable but did not remove this search assumption.

History:

- The original `pbinbf01()` implementation already used continuous Brent
  extrema, continuous `uniroot()` critical values, integer rounding, and
  summation over the inferred success counts.
- PR5 changed the probability summation path from ordinary probability-scale
  summation to log-scale summation and added stability improvements, but it did
  not introduce the continuous boundary-search assumption.

Possible future change, if binomial boundary detection is revisited:

- Prefer discrete support evaluation for binomial boundaries.
- Find success regions by evaluating `logbf(0:n) - log(k)` on the integer
  support, then summing the selected predictive probabilities on the log scale.
- Keep continuous optimization only as an optional fast path if it is validated
  against the discrete support.

Decision:

- No change for PR6. Keep this as a documented pre-existing assumption rather
  than an active follow-up item.

## Priority Fix Order

Fixed by the 2026-06-13/14 search-risk changes:

1. Remove the one-sided `tcrit()` opposite-side retry and document the
   monotonicity invariant.
2. Stop muffling non-limit `tcrit()` warnings in sequential t evaluators.
3. Replace the fixed `|t| <= 256` cap with a predictive-tail cutoff controlled
   by `tail.eps`.
4. Split sequential t boundary search status so failed searches do not enter
   integration as empty `NaN` regions.
5. Stop dropping backend sequential integration `NA`/`NaN` probabilities with
   `na.rm = TRUE`.
6. Use full final certification for one-sided adaptive t roots and predictive
   tail cutoff endpoints.
7. Tighten `.bfpwr_logspace_sum()` so it only drops `-Inf` zero-mass terms and
   errors on invalid non-finite values.
8. Implement typed-invalid sequential sample-size search state so real
   evaluator errors propagate while expected candidate invalidity remains
   structured solver diagnostics.
9. Keep adaptive sample-size search conservative for transient invalidity and
   direct users to exhaustive search when bracketing is invalidated.
10. Give `search = "exhaustive"` full-range semantics before adaptive
    bracketing, with transient invalid candidates skipped and terminal invalid
    candidates stopping the scan.

No active follow-up remains in this document. Issue #10 is documented as a
pre-existing `pbinbf01()` search assumption with no planned PR6 action.
