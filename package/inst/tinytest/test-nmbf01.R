library(tinytest)
library(bfpwr)

expect_true(is.numeric(nmbf01(estimate = 0, se = 1, null = 0, psd = 1, log = FALSE)),
            info = "nmbf01 should return a numeric value")

expect_true(length(nmbf01(estimate = c(-1, 0, 1), se = 1, null = 0, psd = 1, log = FALSE)) == 3,
            info = "nmbf01 should handle vector inputs")

expect_equal(nmbf01(estimate = 0, se = 1, null = 0, psd = 1, log = TRUE),
             log(nmbf01(estimate = 0, se = 1, null = 0, psd = 1, log = FALSE)),
             info = "nmbf01 should return log(nmbf01) when log = TRUE")

expect_equal(
    pnmbf01(k = 3, n = 1000, usd = 1, null = 0, psd = 1,
            dpm = 0.5, dpsd = 0, lower.tail = FALSE),
    2.146670e-34, tolerance = 1e-6,
    info = "pnmbf01 should compute very small upper-tail probabilities directly"
)
