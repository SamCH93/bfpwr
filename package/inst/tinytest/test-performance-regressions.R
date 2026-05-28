library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
bfpwr_exit_if_not_extended("performance regression checks are extended")

## These are deliberately loose wall-clock checks. They are meant to catch
## large algorithmic regressions, not small machine-to-machine timing noise.

pt_low_n <- bfpwr_expect_elapsed_under(
    "low-n one-sided ptbf01 impossible H0 boundary",
    seconds = 3,
    expr = suppressWarnings(
        ptbf01(k = 6, n = 2, dpm = 0, dpsd = 0,
               alternative = "greater", lower.tail = FALSE)
    )
)
expect_equal(
    pt_low_n,
    0,
    info = "low-n one-sided ptbf01 impossible H0 boundary returns zero probability"
)

nt_h0 <- bfpwr_expect_elapsed_under(
    "one-sided JZS true-null ntbf01 example",
    seconds = 5,
    expr = ntbf01(k = 6, power = 0.95, dpm = 0, dpsd = 0,
                  alternative = "greater", lower.tail = FALSE,
                  nrange = c(2, 10000))
)
expect_equal(
    as.numeric(nt_h0),
    4920,
    info = "one-sided JZS true-null ntbf01 example returns the stabilized sample size"
)

jzs_n <- seq(40, 100, 1)
jzs_h1 <- bfpwr_expect_elapsed_under(
    "one-sided JZS sequential t-test paper schedule",
    seconds = 15,
    expr = ptbf01seq(k1 = 1/30, k0 = 6, n = jzs_n, plocation = 0,
                     pscale = 1/sqrt(2), pdf = 1, dpm = 0.5, dpsd = 0.1,
                     type = "two.sample", alternative = "greater")
)
expect_equal(
    c(tail(jzs_h1$cumpH1, 1), tail(jzs_h1$cumpH0, 1), jzs_h1$EN1),
    c(0.7026386, 0.0178849, 69.40229),
    tolerance = 5e-5,
    info = "timed sequential t-test check still matches the paper schedule"
)
