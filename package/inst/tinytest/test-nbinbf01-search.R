library(tinytest)
library(bfpwr)

## Binomial power oscillates with n. nbinbf01() intentionally returns a sample
## size after which power stays at or above the target for the next 10 sample
## sizes, not the first sample size that reaches the target.
binArgs <- list(k = 1/6, power = 0.8, p0 = 0.5, type = "point", a = 1, b = 1,
                dp = 0.7)
binPower <- function(n) {
    pbinbf01(k = 1/6, n = n, p0 = 0.5, type = "point", a = 1, b = 1, dp = 0.7)
}
binN <- function(nrange) do.call(nbinbf01, c(binArgs, list(nrange = nrange)))

nStable <- binN(c(10, 10^4))
expect_equal(nStable, 79)
expect_true(all(binPower(nStable:(nStable + 10)) >= 0.8))
expect_true(binPower(nStable - 1) < 0.8)
expect_true(binPower(74) >= 0.8,
            info = "an earlier sample size reaches the target before power settles")

## searches that cannot return a sample size give NaN with a warning instead
## of an error, as in the continuous-data sample-size functions
expect_warning(binN(c(74, 10^4)), pattern = "lower bound")
expect_true(is.nan(suppressWarnings(binN(c(74, 10^4)))))
expect_warning(binN(c(10, 50)), pattern = "upper bound")
expect_true(is.nan(suppressWarnings(binN(c(10, 50)))))

## power reaches the target at nrange[2] = 77 but falls below it at n = 78
expect_warning(binN(c(10, 77)), pattern = "above the upper bound")
expect_true(is.nan(suppressWarnings(binN(c(10, 77)))))

## the wrapper passes the NaN through
expect_warning(powerbinbf01(power = 0.8, k = 1/6, p0 = 0.5, type = "point",
                            dp = 0.7, nrange = c(74, 10^4)),
               pattern = "lower bound")
wrapped <- suppressWarnings(
    powerbinbf01(power = 0.8, k = 1/6, p0 = 0.5, type = "point", dp = 0.7,
                 nrange = c(74, 10^4))
)
expect_true(is.nan(wrapped$n))
