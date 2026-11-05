# bfpwr 0.3

- new contributor František Bartoš (<https://orcid.org/0000-0002-0018-5573>)
- add `dirbf01()` to compute directional z-test Bayes factors
- add point-null one-sided z-tests with truncated normal analysis priors via
  `alternative = "less"` or `"greater"`
- add one-sided normal moment priors via `alternative = "less"` or `"greater"`
- add sequential Bayes factor designs for z-tests and t-tests via `pbf01seq()`
  and `ptbf01seq()`, with stopping probabilities, expected sample sizes,
  summaries, and plots
- add sequential power and sample-size calculations via `nbf01seq()`,
  `ntbf01seq()`, `powerbf01seq()`, and `powertbf01seq()`
- add `bfpwrOptions()` to inspect and change session numerical settings
- improve numerical accuracy of small probabilities in `pbf01()` and
  `pnmbf01()`; fix `pbf01()` for very narrow normal priors and point
  alternatives identical to the null hypothesis
- improve numerical stability of one-sided `tbf01()` when observations strongly
  oppose the alternative hypothesis
- fix two-sided `ptbf01()` power and critical-value searches for shifted
  informed priors, including heavy-tailed priors
- compute `tbf01()` accurately when the t-statistic and effect sizes under the
  prior have opposite signs, for extreme t-statistics, and for very large
  sample sizes, where `stats::dt()` with a noncentrality parameter is
  unreliable and a single integral over the prior can miss the narrow
  likelihood; this corrects informed-prior and large-sample Bayes factors and
  the power and sample-size calculations that use them
- add `tail.eps` to control the tail-probability cutoff in one-sided adaptive
  t critical-value searches, and `tail.nquad` to control the accuracy/speed
  tradeoff in one-sided t Bayes factor calculations (defaults: `1e-6` and 512)
- correct `pbinbf01()` power calculations for discrete binomial outcomes,
  including equality at the Bayes factor threshold

# bfpwr 0.1.7

- added functions `pbinbf01seq` and `powerbinbf01seq` for sequential two-stage
  designs for binomial Bayes factors with futility stopping
  (<https://doi.org/10.48550/arXiv.2511.23144>)

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
  (<https://doi.org/10.48550/arXiv.2502.02914>)
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
