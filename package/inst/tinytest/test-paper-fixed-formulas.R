library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "paper fixed-design formula checks are extended"
    ))
}

## Checks for values printed in paper/bfssd.Rnw. Each example calls the
## package function and compares the rounded or integer result to the
## manuscript value.

## Mirtazapine example: "mirtazapine-example" and
## "mirtazapine-example-design".
mirt_est <- -1.74
mirt_ci <- c(-7.17, 3.69)
mirt_se <- (mirt_ci[2] - mirt_ci[1])/(2*stats::qnorm(0.975))
mirt_pm <- -6
mirt_usd <- sqrt(2)*15

expect_equal(
    round(bf01(estimate = mirt_est, se = mirt_se, null = 0,
               pm = mirt_pm, psd = 0), 1),
    2.7,
    info = "mirtazapine example BF01 rounds to the manuscript value"
)

expect_equal(
    nbf01(k = 1/10, power = 0.8, usd = mirt_usd, null = 0,
          pm = mirt_pm, psd = 0, dpm = mirt_pm, dpsd = 0),
    124,
    info = "mirtazapine point-design sample size matches the manuscript"
)

expect_equal(
    nbf01(k = 1/10, power = 0.8, usd = mirt_usd, null = 0,
          pm = mirt_pm, psd = 0, dpm = mirt_pm, dpsd = 2),
    195,
    info = "mirtazapine uncertain-design sample size matches the manuscript"
)

expect_equal(
    nbf01(k = 10, power = 0.8, usd = mirt_usd, null = 0,
          pm = mirt_pm, psd = 0, dpm = 0, dpsd = 0,
          lower.tail = FALSE),
    124,
    info = "mirtazapine true-null sample size matches the manuscript symmetry"
)

## Schoenbrodt and Wagenmakers standardized-mean-difference example:
## "BFDA-comparison".
normal_prior_n <- c(
    nbf01(k = 1/6, power = 0.95, usd = sqrt(2), null = 0,
          pm = 0, psd = 1/sqrt(2), dpm = 0.5, dpsd = 0),
    nbf01(k = 1/6, power = 0.95, usd = sqrt(2), null = 0,
          pm = 0, psd = 1/sqrt(2), dpm = 0.5, dpsd = 0.1),
    nbf01(k = 6, power = 0.95, usd = sqrt(2), null = 0,
          pm = 0, psd = 1/sqrt(2), dpm = 0, dpsd = 0,
          lower.tail = FALSE)
)

expect_equal(
    normal_prior_n,
    c(153, 211, 6691),
    info = "normal-prior sample sizes match the manuscript"
)

## One-sided JZS t-test example: "extensions-t-example".
expect_equal(
    as.numeric(ntbf01(k = 1/6, power = 0.95, null = 0,
                      plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                      alternative = "greater", type = "two.sample",
                      dpm = 0.5, dpsd = c(0, 0.1))),
    c(143, 195),
    info = "one-sided JZS t-test sample sizes match the manuscript"
)

expect_equal(
    as.numeric(ntbf01(k = 6, power = 0.95, null = 0,
                      plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                      alternative = "greater", type = "two.sample",
                      dpm = 0, dpsd = 0, lower.tail = FALSE,
                      nrange = c(2, 10000))),
    4920,
    info = "one-sided JZS true-null t-test sample size remains stable"
)

## Normal-moment example: "normal-moment-example".
moment_psd <- 0.5/sqrt(2)
moment_n <- c(
    nnmbf01(k = 1/6, power = 0.95, usd = 2, null = 0,
            psd = moment_psd, dpm = 0.5, dpsd = 0),
    nnmbf01(k = 1/6, power = 0.95, usd = 2, null = 0,
            psd = moment_psd, dpm = 0.5, dpsd = 0.1),
    nnmbf01(k = 6, power = 0.95, usd = 2, null = 0,
            psd = moment_psd, dpm = 0, dpsd = 0,
            lower.tail = FALSE),
    nbf01(k = 6, power = 0.95, usd = sqrt(2), null = 0,
          pm = 0, psd = 1/sqrt(2), dpm = 0, dpsd = 0,
          lower.tail = FALSE)
)

expect_equal(
    moment_n,
    c(302, 429, 997, 6691),
    info = "normal-moment sample sizes match the manuscript"
)
