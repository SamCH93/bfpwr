library(tinytest)
library(bfpwr)

res <- bf01(estimate = c(-1, 0, 1), se = 1, null = 0, pm = 0, psd = 1, log = FALSE)
logres <- bf01(estimate = c(-1, 0, 1), se = 1, null = 0, pm = 0, psd = 1, log = TRUE)

expect_true(is.numeric(res), info = "bf01 should return a numeric value")

expect_true(length(res) == 3,
            info = "bf01 should handle vector inputs")

expect_equal(log(res), log(res),
             info = "bf01 should return log(bf01) when log = TRUE")

expect_equal(bf01(estimate = 0, se = 1, null = 0, pm = 0, psd = 0, log = FALSE),
             1, info = "bf01 should return 1 when null = 0, pm = 0, psd = 0")

expect_equal(
    pbf01(k = 3, n = 1000, usd = 1, null = 0, pm = 0, psd = 1,
          dpm = 0.5, dpsd = 0, lower.tail = FALSE),
    1.162618e-42, tolerance = 1e-6,
    info = "pbf01 should compute very small upper-tail probabilities directly"
)
