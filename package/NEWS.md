# bfpwr 0.3

- preserve the H0 and H1 stopping sides when sequential t boundary searches
  reach a predictive tail cutoff
- stabilize normal-prior critical values near the point-alternative limit and
  avoid cancellation in sequential stopping-sample-size variances
- share numerical integration setup between H0 and H1 events and preserve
  custom integration controls when plotting sequential null curves
- add sequential sample-size search helpers `nbf01seq()`, `ntbf01seq()`,
  `powerbf01seq()`, and `powertbf01seq()`
- add sequential t-test stopping probabilities via `ptbf01seq()`
- add `tail.eps` control for one-sided adaptive t critical-value searches;
  the old fixed `|t| <= 256` stopping cap is replaced by a predictive-tail
  probability cutoff
- fix extreme wrong-tail one-sided `tbf01()` underflow with fixed
  Gauss-Legendre quadrature; expose the accuracy/speed tradeoff via
  `tail.nquad`, defaulting to 128 nodes
- make one-sided adaptive t searches scan only the mathematically expected
  direction and report unresolved finite searches with clearer diagnostics
- fix two-sided `ptbf01()` power for shifted informed priors where BF01 is
  maximized away from the null
- fix two-sided informed t critical-value searches for shifted, heavy-tailed
  priors whose two roots lie on the same side of the old heuristic split
- stabilize extreme one-sided informed t Bayes factors by integrating the
  exact latent-chi-square likelihood ratio over truncated prior quantiles
- make `pbinbf01()` locate inclusive critical counts on the integer data grid
  instead of rounding continuous numerical roots
- handle the identical point alternative in `pbf01()` and enforce strictly
  separated sequential thresholds (`k1 < 1 < k0`)
- make sequential sample-size search preserve typed numerical invalidity,
  propagate structural evaluator errors, and give `search = "exhaustive"`
  full-range semantics
- respect both integer endpoints of sequential `nrange`, retain alternating
  feasible rounded schedules, and scan past transient t-boundary failures in
  exhaustive searches
- expose the deterministic `lpmvnorm` grid size as `ngrid` and record the
  integration settings in sequential design objects
- reject direct sequential schedules with non-increasing information or sample
  sizes
- make fixed-`n` sequential wrapper schedules round-trip searched increment
  schedules by respecting `nrange[1]` as the default first look when `minN` is
  missing
- new contributor František Bartoš (<https://orcid.org/0000-0002-0018-5573>)

# bfpwr 0.2

- add function `pbf01seq` to compute characteristics of sequential Bayes factor
  designs

# bfpwr 0.1.6

- fixed a bug in `ptbf01` when adaptively determining the search range 

# bfpwr 0.1.5

- update reference to published paper in The American Statistician
  (<https://doi.org/10.1080/00031305.2025.2467919>)
- fixed a bug in `plot.power.bftest` so that the correct data frame is returned
  for H0 (thanks Riko Kelter!)

# bfpwr 0.1.4

- new functions related to Bayes factors for testing a binomial proportion:
  `binbf01`, `pbinbf01`, `nbinbf01`, `powerbinbf01`
- improvements in documentation based on feedback from Riko Kelter (thanks!)

# bfpwr 0.1.3

- changed citation in DESCRIPTION file to adhere to CRAN policy
- properly reset user's par() in the vignette

# bfpwr 0.1.2

- polished documentation based on feedback from Tsz Keung Wong (thanks!)
- changed name of `sd` argument to `usd` in `nbf01`, `nnmbf01`, `pbf01`, `pbf01`
  to more clearly differentiate between the standard deviation of the data and
  the unit standard deviation related to a parameter estimate

# bfpwr 0.1

- vignette `Using the bfpwr package` created
- new functions: `bf01`, `tbf01`, `nmbf01`, `nbf01`, `ntbf01`, `nnmbf01`,
  `pbf01`, `ptbf01`, `pnmbf01`, `powerbf01`, `powertbf01`, `powernmbf01`
- package development started
