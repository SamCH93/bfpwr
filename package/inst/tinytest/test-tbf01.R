library(tinytest)
library(bfpwr)

res <- tbf01(t = c(-1, 0, 1), n = 100, plocation = 0, pscale = 1, pdf = 1,
             type = "one.sample", alternative = "two.sided", log = FALSE)
logres <- tbf01(t = c(-1, 0, 1), n = 100, plocation = 0, pscale = 1, pdf = 1,
                type = "one.sample", alternative = "two.sided", log = TRUE)

expect_true(is.numeric(res), info = "tbf01 should return a numeric value")

expect_true(length(res) == 3, info = "tbf01 should handle vector inputs")

expect_equal(log(res), logres,
             info = "tbf01 should return log(tbf01) when log = TRUE")
