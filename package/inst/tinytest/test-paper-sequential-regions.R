library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "paper sequential-region checks are extended"
    ))
}

## Checks for numbers printed in the BFGSD paper. Each example calls the
## package function and compares the result to the value reported in the paper.

expect_numeric_equal <- function(value, expected, tolerance, info) {
    expect_equal(as.numeric(value), expected, tolerance = tolerance, info = info)
}

## Low-PV interim Bayes factors.
lowpv_p0 <- 0.5
lowpv_p1 <- 0.75
lowpv_pm <- log((lowpv_p1/(1 - lowpv_p1))/(lowpv_p0/(1 - lowpv_p0)))
lowpv_psd <- 0

lowpv_logOR1 <- log(21*(26 - 15)/((24 - 21)*15))
lowpv_se1 <- sqrt(1/21 + 1/(24 - 21) + 1/15 + 1/(26 - 15))
lowpv_logOR2 <- log(42*(50 - 30)/((50 - 42)*30))
lowpv_se2 <- sqrt(1/42 + 1/(50 - 42) + 1/30 + 1/(50 - 30))

lowpv_bf <- c(
    bf01(estimate = lowpv_logOR1, se = lowpv_se1, null = 0,
         pm = lowpv_pm, psd = lowpv_psd),
    bf01(estimate = lowpv_logOR2, se = lowpv_se2, null = 0,
         pm = lowpv_pm, psd = lowpv_psd)
)

expect_numeric_equal(
    lowpv_bf,
    c(0.1090023, 0.03582541),
    tolerance = 5e-7,
    info = "Low-PV interim BF01 values match the paper"
)
expect_equal(
    round(1/lowpv_bf, 1),
    c(9.2, 27.9),
    info = "Low-PV reciprocal BF values round to the paper values"
)

## Low-PV three-look design.
lowpv_n <- c(25, 50, 75)
lowpv_k1 <- 1/10
lowpv_k0 <- 10
lowpv_se_h1 <- sqrt(1/(lowpv_p0*(1 - lowpv_p0)*lowpv_n) +
                    1/(lowpv_p1*(1 - lowpv_p1)*lowpv_n))
lowpv_se_h0 <- sqrt(1/(lowpv_p0*(1 - lowpv_p0)*lowpv_n) +
                    1/(lowpv_p0*(1 - lowpv_p0)*lowpv_n))

lowpv_h1 <- pbf01seq(k1 = lowpv_k1, k0 = lowpv_k0, se = lowpv_se_h1,
                     n = lowpv_n, pm = lowpv_pm, psd = lowpv_psd,
                     dpm = lowpv_pm, dpsd = 0, type = "normal")
lowpv_h0 <- pbf01seq(k1 = lowpv_k1, k0 = lowpv_k0, se = lowpv_se_h0,
                     n = lowpv_n, pm = lowpv_pm, psd = lowpv_psd,
                     dpm = 0, dpsd = 0, type = "normal")

expect_numeric_equal(
    lowpv_h1$cumpH1,
    c(0.3513772, 0.6701477, 0.8266490),
    tolerance = 5e-7,
    info = "Low-PV H1-design cumulative H1 probabilities match the paper"
)
expect_numeric_equal(
    lowpv_h1$cumpH0,
    c(0.01464240, 0.02500717, 0.03001821),
    tolerance = 5e-7,
    info = "Low-PV H1-design cumulative H0 probabilities match the paper"
)
expect_numeric_equal(
    lowpv_h1$EN,
    48.47064,
    tolerance = 5e-5,
    info = "Low-PV H1-design expected sample size matches the paper"
)
expect_numeric_equal(
    lowpv_h0$cumpH0,
    c(0.4150487, 0.7304677, 0.8678707),
    tolerance = 5e-7,
    info = "Low-PV H0-design cumulative H0 probabilities match the paper"
)
expect_numeric_equal(
    lowpv_h0$cumpH1,
    c(0.01551580, 0.02464182, 0.02853864),
    tolerance = 5e-7,
    info = "Low-PV H0-design cumulative H1 probabilities match the paper"
)
expect_numeric_equal(
    lowpv_h0$EN,
    45.35815,
    tolerance = 5e-5,
    info = "Low-PV H0-design expected sample size matches the paper"
)

## BFGSD appendix, one-sided JZS sequential t design.
## The paper code uses step <- 1 from n = 40 to n = 100, which is one new
## observation per group at each interim look. This is now fast enough to keep
## as an extended reference check.
jzs_n <- seq(40, 100, 1)
jzs_h1 <- ptbf01seq(k1 = 1/30, k0 = 6, n = jzs_n, plocation = 0,
                    pscale = 1/sqrt(2), pdf = 1, dpm = 0.5, dpsd = 0.1,
                    type = "two.sample", alternative = "greater")
jzs_h0 <- ptbf01seq(k1 = 1/30, k0 = 6, n = jzs_n, plocation = 0,
                    pscale = 1/sqrt(2), pdf = 1, dpm = 0, dpsd = 0,
                    type = "two.sample", alternative = "greater")

expect_equal(length(jzs_n), 61)
expect_numeric_equal(
    c(tail(jzs_h1$cumpH1, 1), tail(jzs_h1$cumpH0, 1), jzs_h1$EN1),
    c(0.7026386, 0.0178849, 69.40229),
    tolerance = 5e-5,
    info = "JZS H1-design final values match the paper schedule"
)
expect_numeric_equal(
    c(tail(jzs_h0$cumpH0, 1), tail(jzs_h0$cumpH1, 1), jzs_h0$EN1),
    c(0.7129341, 0.004836611, 65.74841),
    tolerance = 5e-5,
    info = "JZS H0-design final values match the paper schedule"
)
