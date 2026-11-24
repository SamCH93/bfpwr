library(tinytest)
library(bfpwr)

## The closed-form sample size for a point analysis prior does not use nrange;
## a result outside nrange is still returned but flagged.
pointDesign <- function(effect, nrange, ...) {
    nbf01(k = 1/6, power = 0.8, usd = 1, null = 0, pm = effect, psd = 0,
          dpm = effect, dpsd = 0, nrange = nrange, ...)
}

## required sample size above the upper bound
expect_warning(pointDesign(0.3, nrange = c(10, 50)),
               pattern = "outside the sample size search range")
expect_equal(suppressWarnings(pointDesign(0.3, nrange = c(10, 50))),
             pointDesign(0.3, nrange = c(1, 10^5)))
expect_true(is.nan(suppressWarnings(pointDesign(0.3, nrange = c(10, 50),
                                                analytical = FALSE))))

## required sample size below the lower bound
expect_warning(pointDesign(1, nrange = c(10, 50)),
               pattern = "outside the sample size search range")

## no warning inside the range, where both methods agree
expect_silent(pointDesign(0.3, nrange = c(10, 200)))
expect_equal(pointDesign(0.3, nrange = c(10, 200), integer = FALSE),
             pointDesign(0.3, nrange = c(10, 200), integer = FALSE,
                         analytical = FALSE),
             tolerance = 1e-6)
