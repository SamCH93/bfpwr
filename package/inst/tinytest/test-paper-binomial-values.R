library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "paper binomial-value checks are extended"
    ))
}

## Checks for values printed in Kelter and Pawel (2025), "Sample Size
## Determination for Bayes Factor Analysis of Binomial Data",
## arXiv:2502.02914. Each example calls the package function for that value.

## Single-arm phase II proof-of-concept trial, flat directional priors.
expect_equal(
    nbinbf01(k = 1/10, power = 0.9, p0 = 0.2, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0.2, du = 1),
    110,
    info = "phase II strong-evidence H1 sample size matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 110, p0 = 0.2, type = "direction",
                         a = 1, b = 1, da = 1, db = 1, dl = 0.2, du = 1), 2),
    90.05,
    info = "phase II strong-evidence H1 power matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 110, p0 = 0.2, type = "direction",
                         a = 1, b = 1, da = 1, db = 1, dl = 0, du = 0.2), 2),
    0.16,
    info = "phase II strong-evidence H0 type-I error matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 110, p0 = 0.2, type = "direction",
                         a = 1, b = 1, dp = 0.4), 2),
    99.63,
    info = "phase II frequentist power at p1 = 0.4 matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 110, p0 = 0.2, type = "direction",
                         a = 1, b = 1, dp = 0.2), 2),
    2.47,
    info = "phase II frequentist type-I error at p0 = 0.2 matches the paper"
)

expect_equal(
    nbinbf01(k = 10, power = 0.9, p0 = 0.2, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0, du = 0.2,
             lower.tail = FALSE),
    245,
    info = "phase II strong-evidence H0 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 1/3, power = 0.9, p0 = 0.2, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0.2, du = 1),
    61,
    info = "phase II moderate-evidence H1 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 3, power = 0.9, p0 = 0.2, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0, du = 0.2,
             lower.tail = FALSE),
    60,
    info = "phase II moderate-evidence H0 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 1/3, power = 0.9, p0 = 0.2, type = "direction",
             a = 1, b = 1, dp = 0.4),
    36,
    info = "phase II point-design moderate-evidence H1 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 1/10, power = 0.9, p0 = 0.2, type = "direction",
             a = 1, b = 1, dp = 0.4),
    53,
    info = "phase II point-design strong-evidence H1 sample size matches the paper"
)

## Therapeutic touch experiment: 70 correct responses out of 150 trials.
expect_equal(
    round(binbf01(x = 70, n = 150, p0 = 0.5, type = "point",
                  a = 1, b = 1), 2),
    7.05,
    info = "therapeutic touch point-null BF01 matches the paper"
)

expect_equal(
    round(binbf01(x = 70, n = 150, p0 = 0.5, type = "direction",
                  a = 1, b = 1), 2),
    3.81,
    info = "therapeutic touch directional BF01 matches the paper"
)

expect_equal(
    nbinbf01(k = 1/10, power = 0.8, p0 = 0.5, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0.5, du = 1),
    50,
    info = "therapeutic touch directional H1 sample size matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 50, p0 = 0.5, type = "direction",
                         a = 1, b = 1, da = 1, db = 1, dl = 0.5, du = 1), 2),
    81.68,
    info = "therapeutic touch directional H1 power matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 50, p0 = 0.5, type = "direction",
                         a = 1, b = 1, da = 1, db = 1, dl = 0, du = 0.5), 3),
    0.674,
    info = "therapeutic touch directional H0 type-I error matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 50, p0 = 0.5, type = "direction",
                         a = 1, b = 1, dp = 0.5), 2),
    10.13,
    info = "therapeutic touch directional point-design type-I error matches the paper"
)

expect_equal(
    nbinbf01(k = 3.81, power = 0.8, p0 = 0.5, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0, du = 0.5,
             lower.tail = FALSE),
    27,
    info = "therapeutic touch directional retrospective H0 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 3, power = 0.8, p0 = 0.5, type = "direction",
             a = 1, b = 1, da = 1, db = 1, dl = 0, du = 0.5,
             lower.tail = FALSE),
    22,
    info = "therapeutic touch directional moderate-evidence H0 sample size matches the footnote"
)

expect_equal(
    nbinbf01(k = 1/10, power = 0.8, p0 = 0.5, type = "point",
             a = 1, b = 1, da = 1, db = 1, dl = 0, du = 1),
    245,
    info = "therapeutic touch two-sided strong-evidence H1 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 10, power = 0.8, p0 = 0.5, type = "point",
             a = 1, b = 1, dp = 0.5, lower.tail = FALSE,
             nrange = c(1, 1000)),
    853,
    info = "therapeutic touch two-sided strong-evidence H0 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 1/3, power = 0.8, p0 = 0.5, type = "point",
             a = 1, b = 1, da = 1, db = 1, dl = 0, du = 1),
    180,
    info = "therapeutic touch two-sided moderate-evidence H1 sample size matches the paper"
)

expect_equal(
    nbinbf01(k = 3, power = 0.8, p0 = 0.5, type = "point",
             a = 1, b = 1, dp = 0.5, lower.tail = FALSE,
             nrange = c(1, 1000)),
    90,
    info = "therapeutic touch two-sided moderate-evidence H0 sample size matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/10, n = 150, p0 = 0.5, type = "point",
                         a = 1, b = 1, da = 1, db = 1, dl = 0, du = 1), 2),
    75.50,
    info = "therapeutic touch two-sided observed-sample strong-evidence power matches the paper"
)

expect_equal(
    round(100 * pbinbf01(k = 1/3, n = 150, p0 = 0.5, type = "point",
                         a = 1, b = 1, da = 1, db = 1, dl = 0, du = 1), 2),
    79.47,
    info = "therapeutic touch two-sided observed-sample moderate-evidence power matches the paper"
)
