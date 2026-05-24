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

expect_equal(
    tbf01(t = -20, n1 = 7880, n2 = 7880, alternative = "greater",
          type = "two.sample", log = TRUE),
    7.2302101, tolerance = 1e-6,
    info = "tbf01 should use a stable wrong-tail one-sided calculation"
)

expect_equal(
    tbf01(t = -4, n1 = 53, n2 = 57, plocation = 0.350, pscale = 0.102,
          pdf = 3, alternative = "greater", type = "two.sample", log = TRUE),
    4.3357339, tolerance = 1e-6,
    info = "tbf01 should handle shifted informed priors in the wrong tail"
)

expect_equal(
    tbf01(t = 1, n1 = 1e6, n2 = 1e6, pscale = 1, type = "two.sample",
          log = TRUE),
    6.2869769, tolerance = 1e-4,
    info = "tbf01 should return finite large-n two-sided log Bayes factors"
)
