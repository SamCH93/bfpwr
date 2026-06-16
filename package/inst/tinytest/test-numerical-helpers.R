library(tinytest)
library(bfpwr)

logs <- log(c(0.2, 0.3))
expect_equal(
    bfpwr:::.bfpwr_logspace_sum(c(logs, -Inf)),
    log(0.5),
    tolerance = 1e-12,
    info = "log-space summation should ignore -Inf zero-mass terms"
)

expect_equal(
    bfpwr:::.bfpwr_logspace_sum(c(-Inf, -Inf)),
    -Inf,
    info = "log-space summation should return -Inf when all terms are zero mass"
)

expect_error(
    bfpwr:::.bfpwr_logspace_sum(c(logs, NA_real_)),
    "invalid non-finite",
    info = "log-space summation should reject NA values"
)

expect_error(
    bfpwr:::.bfpwr_logspace_sum(c(logs, NaN)),
    "invalid non-finite",
    info = "log-space summation should reject NaN values"
)

expect_error(
    bfpwr:::.bfpwr_logspace_sum(c(logs, Inf)),
    "invalid non-finite",
    info = "log-space summation should reject positive infinity"
)

expect_equal(
    bfpwr:::.bfpwr_lpnorm_interval(
        lower = c(-Inf, -1, 1),
        upper = c(0, 1, Inf),
        mean = c(0, 0.5, 0),
        sd = c(1, 2, 1)
    ),
    log(c(
        stats::pnorm(0),
        stats::pnorm(1, mean = 0.5, sd = 2) -
            stats::pnorm(-1, mean = 0.5, sd = 2),
        stats::pnorm(1, lower.tail = FALSE)
    )),
    tolerance = 1e-12,
    info = "normal interval helper should handle vector inputs"
)
