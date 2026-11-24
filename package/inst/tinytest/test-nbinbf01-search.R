library(tinytest)
library(bfpwr)

## Binomial power can oscillate with n. nbinbf01() intentionally checks that
## power does not drop below the target for the next sample sizes after a
## crossing instead of returning the first sample size that reaches the target.
## In this example, an earlier size (74) reaches the target, power dips at 78,
## and the search returns 79 with power at or above the target through 89.
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
expect_true(binPower(77) >= 0.8)
expect_silent(nAbove <- binN(c(10, 77)))
expect_equal(nAbove, 79)

## when only the crossing itself is below the target, the next sample size is
## returned and only 9 further sample sizes are checked; power can drop at
## the next unchecked size
dirPower <- function(n) {
    pbinbf01(k = 1/6, n = n, p0 = 0.5, type = "direction", dp = 0.55)
}
expect_silent(nDir <- nbinbf01(k = 1/6, power = 0.75, p0 = 0.5,
                               type = "direction", dp = 0.55,
                               nrange = c(1, 400)))
expect_equal(nDir, 307)
expect_true(dirPower(nDir - 1) < 0.75)
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

## restarts exhausted: after the initial search and 5 restarts, warn and
## return the last candidate (searchN() stubbed to make this deterministic)
searchCalls <- 0
searchStub <- function(...) {
    searchCalls <<- searchCalls + 1
    c(10, 20, 30, 40, 50, 60)[searchCalls]
}
searchNoscilStub <- bfpwr:::searchNoscil
environment(searchNoscilStub) <- list2env(
    list(searchN = searchStub), parent = environment(bfpwr:::searchNoscil)
)
dipAt9 <- function(n) ifelse(n %% 10 == 9, -1, 0)
expect_warning(nExhausted <- searchNoscilStub(dipAt9, nrange = c(1, 1000)),
               pattern = "extend root-finding cycles")
expect_equal(nExhausted, 60)
expect_equal(searchCalls, 6)

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
grDevices::pdf(NULL)
expect_warning(curves <- plot(designH1), pattern = "upper bound")
invisible(grDevices::dev.off())
expect_equal(names(curves), c("powDFH1", "powDFH0"))
expect_true(nrow(curves$powDFH0) > 0 && all(is.finite(curves$powDFH0$power)))
