library(tinytest)
library(bfpwr)

## powertbf01() should honor the caller's sample-size search range.
wide_t <- suppressWarnings(
    powertbf01(k = 1/10, power = 0.7, plocation = 0, pscale = 0.707,
               pdf = 1, type = "two.sample", alternative = "greater",
               dpm = 0.05, dpsd = 0, nrange = c(2, 20000))
)
expect_true(is.finite(wide_t$n) && wide_t$n > 10000 && wide_t$n < 20000,
            info = "powertbf01 should pass custom nrange to ntbf01")

## H0 plot marker calculations use a point null design prior, even for nonzero
## nulls. A negative null used to be passed as dpsd and failed validation.
z_design <- suppressWarnings(
    powerbf01(k = 1/10, power = 0.5, null = -0.1, pm = 0.2, psd = 0.5,
              dpm = 0.2, dpsd = 0, nrange = c(1, 10000))
)
z_plot <- try(plot(z_design, plot = FALSE, nlim = c(10, 100), ngrid = 5),
              silent = TRUE)
expect_false(inherits(z_plot, "try-error"),
             info = "plot.power.bftest should handle negative null H0 marker")
