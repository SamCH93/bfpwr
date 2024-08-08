library(tinytest)
library(bfpwr)

expect_true(is.numeric(nmbf01(estimate = 0, se = 1, null = 0, psd = 1, log = FALSE)),
            info = "nmbf01 should return a numeric value")

expect_true(length(nmbf01(estimate = c(-1, 0, 1), se = 1, null = 0, psd = 1, log = FALSE)) == 3,
            info = "nmbf01 should handle vector inputs")

expect_equal(nmbf01(estimate = 0, se = 1, null = 0, psd = 1, log = TRUE),
             log(nmbf01(estimate = 0, se = 1, null = 0, psd = 1, log = FALSE)),
             info = "nmbf01 should return log(nmbf01) when log = TRUE")
