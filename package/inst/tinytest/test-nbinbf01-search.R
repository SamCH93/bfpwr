library(tinytest)
library(bfpwr)

## Binomial power can oscillate with n. nbinbf01() intentionally checks that
## power does not drop below the target for the next sample sizes after a
## crossing instead of returning the first sample size that reaches the target.
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

## the window check can extend beyond nrange[2]: power reaches the target at
## nrange[2] = 77 but drops at 78, and the search returns 79 without a warning
expect_silent(nAbove <- binN(c(10, 77)))
expect_equal(nAbove, 79)

## when only the crossing itself is below the target, the next sample size is
## returned and only 9 further sample sizes are checked; power can drop after
dirPower <- function(n) {
    pbinbf01(k = 1/6, n = n, p0 = 0.5, type = "direction", dp = 0.55)
}
nDir <- nbinbf01(k = 1/6, power = 0.75, p0 = 0.5, type = "direction",
                 dp = 0.55, nrange = c(1, 400))
expect_equal(nDir, 307)
expect_true(all(dirPower(nDir:(nDir + 9)) >= 0.75))
expect_true(dirPower(nDir + 10) < 0.75)

## searches that fail give NaN with the search warning instead of an error, as
## in the continuous-data sample-size functions
expect_warning(binN(c(74, 10^4)), pattern = "lower bound")
expect_true(is.nan(suppressWarnings(binN(c(74, 10^4)))))
expect_warning(binN(c(10, 50)), pattern = "upper bound")
expect_true(is.nan(suppressWarnings(binN(c(10, 50)))))

## vectorised input keeps the successful element
expect_equal(suppressWarnings(
    nbinbf01(k = 1/6, power = c(0.8, 0.95), p0 = 0.5, type = "point",
             dp = 0.7, nrange = c(10, 100))
), c(79, NaN))

## restarts exhausted: warn and return the last crossing
dipEvery7 <- function(n) ifelse(ceiling(n) >= 10 & ceiling(n) %% 7 != 0, 1, -1)
expect_warning(bfpwr:::searchNoscil(dipEvery7, nrange = c(1, 1000)),
               pattern = "extend root-finding cycles")
expect_true(is.finite(suppressWarnings(
    bfpwr:::searchNoscil(dipEvery7, nrange = c(1, 1000))
)))

## downstream: powerbinbf01() passes the NaN through, and plotting no longer
## fails when the H0 search (lower.tail = FALSE) fails
expect_warning(powerbinbf01(power = 0.8, k = 1/6, p0 = 0.5, type = "point",
                            dp = 0.7, nrange = c(74, 10^4)),
               pattern = "lower bound")
wrapped <- suppressWarnings(
    powerbinbf01(power = 0.8, k = 1/6, p0 = 0.5, type = "point", dp = 0.7,
                 nrange = c(74, 10^4))
)
expect_true(is.nan(wrapped$n))
designH1 <- powerbinbf01(power = 0.8, k = 1/6, p0 = 0.5, type = "point",
                         dp = 0.7, nrange = c(10, 100))
expect_equal(designH1$n, 79)
expect_warning(curves <- plot(designH1, plot = FALSE), pattern = "upper bound")
expect_equal(names(curves), c("powDFH1", "powDFH0"))
