library(tinytest)
library(bfpwr)

## Legacy unnamed numerical arguments continue through ... to searchN.
legacy <- nnmbf01(.1, .8, 1, 0, .7, .5, 0, c(1, 1e5), TRUE, TRUE, 1e-10)
named <- nnmbf01(.1, .8, 1, 0, .7, .5, 0, c(1, 1e5), TRUE, TRUE, tol = 1e-10)
expect_equal(legacy, named)
expect_equal(legacy, 68)
expect_equal(nnmbf01(.1, .8, 1, 0, .7, .5, 0, c(1, 1e5), TRUE, TRUE,
    alternative = "greater", tol = 1e-10), 60)
expect_equal(nnmbf01(.1, .8, 1, 0, .7, c(-.5, .5), 0, c(1, 1e5), TRUE, TRUE,
    alternative = c("less", "greater"), tol = 1e-10), c(60, 60))
