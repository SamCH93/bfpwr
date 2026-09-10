# bfpwr 0.3

- add sequential Bayes factor designs for z-tests and t-tests via `pbf01seq()`
  and `ptbf01seq()`, with stopping probabilities, expected sample sizes,
  summaries, and plots
- add sequential power and sample-size calculations via `nbf01seq()`,
  `ntbf01seq()`, `powerbf01seq()`, and `powertbf01seq()`
- add `dirbf01()` to compute directional z-test Bayes factors
- improve numerical accuracy of small probabilities in `pbf01()` and
  `pnmbf01()`; fix `pbf01()` for very narrow normal priors and point
  alternatives identical to the null hypothesis
- improve numerical stability of one-sided `tbf01()` when observations strongly
  oppose the alternative hypothesis
- fix two-sided `ptbf01()` power and critical-value searches for shifted
  informed priors, including heavy-tailed priors
- add `tail.eps` to control the tail-probability cutoff in one-sided adaptive
  t critical-value searches, and `tail.nquad` to control the accuracy/speed
  tradeoff in one-sided t Bayes factor calculations
- correct `pbinbf01()` power calculations for discrete binomial outcomes,
  including equality at the Bayes factor threshold
- new contributor František Bartoš (<https://orcid.org/0000-0002-0018-5573>)

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
