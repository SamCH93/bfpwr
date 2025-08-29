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
