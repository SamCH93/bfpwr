library(tinytest)
library(bfpwr)

res <- dirbf01(estimate = 0.2, se = 0.2, null = 0, pm = 0, psd = 2)
logres <- dirbf01(estimate = 0.2, se = 0.2, null = 0, pm = 0, psd = 2,
                  log = TRUE)

expect_equal(log(res), logres,
             info = "dirbf01 should return log(dirbf01) when log = TRUE")

extreme <- dirbf01(estimate = 10, se = 1, null = 0, pm = 10, psd = 1,
                   log = TRUE)
expect_true(is.finite(extreme),
            info = "dirbf01 should compute directional odds in extreme tails")
