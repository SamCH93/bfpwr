library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "performance regression checks are extended"
    ))
}

## These are deliberately loose wall-clock checks. They are meant to catch
## large algorithmic regressions, not small machine-to-machine timing noise.

expect_finite_numeric <- function(value, label) {
    expect_true(
        is.numeric(value) && length(value) >= 1 && all(is.finite(value)),
        info = paste(label, "returns finite numeric values")
    )
}

expect_probability <- function(value, label) {
    expect_true(
        is.numeric(value) && length(value) >= 1 &&
            all(is.finite(value)) &&
            all(value >= -1e-12) && all(value <= 1 + 1e-12),
        info = paste(label, "returns probabilities")
    )
}

expect_power_object <- function(value, label) {
    expect_true(inherits(value, "power.bftest"),
                info = paste(label, "returns a power.bftest object"))
    expect_probability(value$power, paste(label, "power"))
    expect_finite_numeric(value$n, paste(label, "n"))
}

expect_seq_object <- function(value, label) {
    expect_true(inherits(value, "bfseqdesign"),
                info = paste(label, "returns a bfseqdesign object"))
    expect_probability(value$cumpH1, paste(label, "cumpH1"))
    expect_probability(value$cumpH0, paste(label, "cumpH0"))
    expect_probability(value$cumpInc, paste(label, "cumpInc"))
}

## Bayes factor evaluators.
bf_z <- bfpwr_expect_elapsed_under(
    "bf01 point-null z BF",
    seconds = 1,
    expr = bf01(estimate = 0.2, se = 0.05, null = 0, pm = 0, psd = 2)
)
expect_finite_numeric(bf_z, "bf01")

bf_dir <- bfpwr_expect_elapsed_under(
    "dirbf01 directional z BF",
    seconds = 1,
    expr = dirbf01(estimate = 0.2, se = 0.2, null = 0, pm = 0, psd = 2)
)
expect_finite_numeric(bf_dir, "dirbf01")

bf_nm <- bfpwr_expect_elapsed_under(
    "nmbf01 normal-moment BF",
    seconds = 1,
    expr = nmbf01(estimate = 0.25, se = 0.05, null = 0,
                  psd = 0.5/sqrt(2))
)
expect_finite_numeric(bf_nm, "nmbf01")

bf_t <- bfpwr_expect_elapsed_under(
    "tbf01 two-value JZS BF",
    seconds = 3,
    expr = tbf01(t = c(0.69, 3.20), n = 100, pscale = 1,
                 type = "one.sample")
)
expect_finite_numeric(bf_t, "tbf01")

bf_bin <- bfpwr_expect_elapsed_under(
    "binbf01 point-null binomial BF",
    seconds = 1,
    expr = binbf01(x = 70, n = 150, p0 = 0.5, type = "point",
                   a = 1, b = 1)
)
expect_finite_numeric(bf_bin, "binbf01")

## Fixed-sample power/CDF functions.
p_z <- bfpwr_expect_elapsed_under(
    "pbf01 z-test power",
    seconds = 1,
    expr = pbf01(k = 1/10, n = 200, usd = 2, null = 0,
                 pm = 0.5, psd = 0)
)
expect_probability(p_z, "pbf01")

p_nm <- bfpwr_expect_elapsed_under(
    "pnmbf01 normal-moment power",
    seconds = 1,
    expr = pnmbf01(k = 1/10, n = 200, usd = 2, null = 0,
                   psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0)
)
expect_probability(p_nm, "pnmbf01")

p_t <- bfpwr_expect_elapsed_under(
    "ptbf01 one-sided t-test H1 power",
    seconds = 4,
    expr = ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0,
                  alternative = "greater")
)
expect_probability(p_t, "ptbf01")

p_t_impossible <- bfpwr_expect_elapsed_under(
    "ptbf01 low-n impossible H0 boundary",
    seconds = 3,
    expr = suppressWarnings(
        ptbf01(k = 6, n = 2, dpm = 0, dpsd = 0,
               alternative = "greater", lower.tail = FALSE)
    )
)
expect_equal(
    p_t_impossible,
    0,
    info = "low-n one-sided ptbf01 impossible H0 boundary returns zero probability"
)

p_bin <- bfpwr_expect_elapsed_under(
    "pbinbf01 directional binomial power",
    seconds = 3,
    expr = pbinbf01(k = 1/10, n = 50, p0 = 0.5, type = "direction",
                    a = 1, b = 1, da = 1, db = 1, dl = 0.5, du = 1)
)
expect_probability(p_bin, "pbinbf01")

## Sample-size searches.
n_z <- bfpwr_expect_elapsed_under(
    "nbf01 z-test sample-size search",
    seconds = 2,
    expr = nbf01(k = 1/10, power = 0.8, usd = 1, null = 0,
                 pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
                 nrange = c(1, 2000))
)
expect_finite_numeric(n_z, "nbf01")

n_nm <- bfpwr_expect_elapsed_under(
    "nnmbf01 normal-moment sample-size search",
    seconds = 3,
    expr = nnmbf01(k = 1/10, power = 0.9, usd = 1, null = 0,
                   psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0)
)
expect_finite_numeric(n_nm, "nnmbf01")

n_t <- bfpwr_expect_elapsed_under(
    "ntbf01 one-sided true-null t-test sample-size search",
    seconds = 8,
    expr = ntbf01(k = 6, power = 0.95, dpm = 0, dpsd = 0,
                  alternative = "greater", lower.tail = FALSE,
                  nrange = c(2, 10000))
)
expect_equal(
    as.numeric(n_t),
    4920,
    info = "one-sided JZS true-null ntbf01 example returns the stabilized sample size"
)

n_bin <- bfpwr_expect_elapsed_under(
    "nbinbf01 directional binomial sample-size search",
    seconds = 15,
    expr = nbinbf01(k = 1/10, power = 0.8, p0 = 0.2,
                    type = "direction", dl = 0.2,
                    nrange = c(1, 250))
)
expect_finite_numeric(n_bin, "nbinbf01")

## User-facing power wrappers, covering both power-at-n and sample-size modes.
power_z <- bfpwr_expect_elapsed_under(
    "powerbf01 power-at-n wrapper",
    seconds = 2,
    expr = powerbf01(n = 100, pm = 0, psd = 1, dpm = 0.5, dpsd = 0)
)
expect_power_object(power_z, "powerbf01 n mode")

size_z <- bfpwr_expect_elapsed_under(
    "powerbf01 sample-size wrapper",
    seconds = 4,
    expr = powerbf01(power = 0.8, pm = 0, psd = 1, dpm = 0.5,
                     dpsd = 0, nrange = c(1, 2000))
)
expect_power_object(size_z, "powerbf01 power mode")

power_nm <- bfpwr_expect_elapsed_under(
    "powernmbf01 power-at-n wrapper",
    seconds = 2,
    expr = powernmbf01(n = 100, psd = 1, dpm = 0.5, dpsd = 0)
)
expect_power_object(power_nm, "powernmbf01 n mode")

size_nm <- bfpwr_expect_elapsed_under(
    "powernmbf01 sample-size wrapper",
    seconds = 5,
    expr = powernmbf01(power = 0.8, psd = 1, dpm = 0.5, dpsd = 0,
                       nrange = c(1, 2000))
)
expect_power_object(size_nm, "powernmbf01 power mode")

power_t <- bfpwr_expect_elapsed_under(
    "powertbf01 power-at-n wrapper",
    seconds = 4,
    expr = powertbf01(n = 146, k = 1/6, dpm = 0.5, dpsd = 0,
                      alternative = "greater")
)
expect_power_object(power_t, "powertbf01 n mode")

size_t <- bfpwr_expect_elapsed_under(
    "powertbf01 sample-size wrapper",
    seconds = 8,
    expr = powertbf01(power = 0.8, k = 1/6, dpm = 0.5, dpsd = 0,
                      alternative = "greater", nrange = c(2, 1000))
)
expect_power_object(size_t, "powertbf01 power mode")

power_bin <- bfpwr_expect_elapsed_under(
    "powerbinbf01 power-at-n wrapper",
    seconds = 3,
    expr = powerbinbf01(n = 100, type = "point")
)
expect_power_object(power_bin, "powerbinbf01 n mode")

size_bin <- bfpwr_expect_elapsed_under(
    "powerbinbf01 sample-size wrapper",
    seconds = 15,
    expr = powerbinbf01(power = 0.8, p0 = 0.2, type = "direction",
                        dl = 0.2, nrange = c(1, 250))
)
expect_power_object(size_bin, "powerbinbf01 power mode")

## Sequential designs.
z_n <- seq(50, 200, 50)
seq_z <- bfpwr_expect_elapsed_under(
    "pbf01seq z-test sequential design",
    seconds = 5,
    expr = pbf01seq(k1 = 1/10, k0 = 3, se = sqrt(2/z_n), n = z_n,
                    pm = 0, psd = 1, dpm = 0.5, dpsd = 0.05,
                    type = "normal")
)
expect_seq_object(seq_z, "pbf01seq")

jzs_n <- seq(40, 100, 1)
seq_t <- bfpwr_expect_elapsed_under(
    "ptbf01seq one-sided JZS sequential design",
    seconds = 15,
    expr = ptbf01seq(k1 = 1/30, k0 = 6, n = jzs_n,
                     plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                     dpm = 0.5, dpsd = 0.1, type = "two.sample",
                     alternative = "greater")
)
expect_seq_object(seq_t, "ptbf01seq")
expect_equal(
    c(tail(seq_t$cumpH1, 1), tail(seq_t$cumpH0, 1), seq_t$EN1),
    c(0.7026386, 0.0178849, 69.40229),
    tolerance = 5e-5,
    info = "timed sequential t-test check still matches the paper schedule"
)

## S3 print and plot methods. Use plot = FALSE to exercise the method while
## avoiding graphics-device noise in timing measurements.
printed_power <- bfpwr_expect_elapsed_under(
    "print.power.bftest method",
    seconds = 1,
    expr = capture.output(print(size_z))
)
expect_true(length(printed_power) > 0,
            info = "print.power.bftest produces console output")

plotted_power <- bfpwr_expect_elapsed_under(
    "plot.power.bftest method",
    seconds = 5,
    expr = plot(size_z, nlim = c(2, 50), ngrid = 10, plot = FALSE,
                nullplot = FALSE)
)
expect_true(is.list(plotted_power) && length(plotted_power) > 0,
            info = "plot.power.bftest returns plotting data")

printed_seq <- bfpwr_expect_elapsed_under(
    "print.bfseqdesign method",
    seconds = 1,
    expr = capture.output(print(seq_z))
)
expect_true(length(printed_seq) > 0,
            info = "print.bfseqdesign produces console output")

plotted_seq <- bfpwr_expect_elapsed_under(
    "plot.bfseqdesign method",
    seconds = 2,
    expr = plot(seq_z, plot = FALSE, nullplot = FALSE)
)
expect_true(is.list(plotted_seq) && length(plotted_seq) > 0,
            info = "plot.bfseqdesign returns plotting data")
