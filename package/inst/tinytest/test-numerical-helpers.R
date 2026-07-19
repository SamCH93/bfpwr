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

## Opposite-direction roots that bracket the proposed split certify a genuine
## two-sided pair. Certification should not evaluate the function again at the
## split, which can be expensive for t Bayes factors.
crossing_fun <- function(x) {
    if (length(x) == 1 && isTRUE(all.equal(x, 0))) {
        stop("the proposed split was evaluated")
    }
    1 - x^2
}
root_pair <- bfpwr:::.bfpwr_two_sided_root_pair(
    f = crossing_fun,
    lowerInterval = c(-2, -0.2),
    upperInterval = c(0.2, 2),
    split = 0
)
expect_true(root_pair$valid,
            info = "opposite roots should bracket the proposed split")
expect_equal(c(root_pair$lower, root_pair$upper), c(-1, 1),
             tolerance = 1e-6,
             info = "two-sided root helper should return both crossings")

outside_pair <- bfpwr:::.bfpwr_two_sided_root_pair(
    f = function(x) 1 - x^2,
    lowerInterval = c(-2, -0.2),
    upperInterval = c(0.2, 2),
    split = 2
)
expect_false(outside_pair$valid,
             info = "roots should be rejected when they do not bracket split")
