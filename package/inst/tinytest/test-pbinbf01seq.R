library(tinytest)
library(bfpwr)

## Tests the two-stage binomial power functions pbinbf01seq() and
## powerbinbf01seq(). The study continues after the interim stage if
## BF01(x1) <= kf and has efficacy at the final stage if BF01(x1 + x2) <= k,
## where x1 and x2 are the successes among the first n1 and the remaining
## n2 - n1 observations. No bfssd manuscript formula covers binomial BF power;
## the reference is therefore the exhaustive enumeration of all (x1, x2) pairs
## with the same floating-point predicate on the log Bayes factors as
## pbinbf01(), weighted with the predictive masses under the design prior.

ref_pbinbf01seq <- function(n1, n2, k, kf, p0 = 0.5, type = "direction",
                            a = 1, b = 1, dp = NA, da = a, db = b, dl = 0,
                            du = 1) {
    stopifnot(type == "direction")
    g <- expand.grid(x1 = 0:n1, x2 = 0:(n2 - n1))
    s <- g$x1 + g$x2
    logbf1 <- binbf01(x = g$x1, n = n1, p0 = p0, type = "direction", a = a,
                      b = b, log = TRUE)
    logbf2 <- binbf01(x = s, n = n2, p0 = p0, type = "direction", a = a,
                      b = b, log = TRUE)
    event <- logbf1 <= log(kf) & logbf2 <= log(k)
    if (!is.na(dp)) {
        logp <- dbinom(g$x1, size = n1, prob = dp, log = TRUE) +
            dbinom(g$x2, size = n2 - n1, prob = dp, log = TRUE)
    } else {
        logint <- vapply(s, function(si) {
            bfpwr:::.bfpwr_lpbeta_interval(lower = dl, upper = du,
                                           shape1 = da + si,
                                           shape2 = db + n2 - si)
        }, numeric(1))
        lognorm <- bfpwr:::.bfpwr_lpbeta_interval(lower = dl, upper = du,
                                                  shape1 = da, shape2 = db)
        logp <- lchoose(n1, g$x1) + lchoose(n2 - n1, g$x2) +
            lbeta(da + s, db + n2 - s) - lbeta(da, db) + logint - lognorm
    }
    sum(exp(logp[event]))
}

expect_matches_reference <- function(label, ...) {
    args <- list(...)
    expect_equal(
        do.call(pbinbf01seq, args), do.call(ref_pbinbf01seq, args),
        tolerance = 1e-12,
        info = paste(label, "should match the enumeration of all (x1, x2)")
    )
}

## point design priors, under the alternative and under the null
expect_matches_reference("point H1", n1 = 10, n2 = 25, k = 1/10, kf = 3,
                         p0 = 0.3, dp = 0.5)
expect_matches_reference("point H0", n1 = 10, n2 = 25, k = 1/10, kf = 3,
                         p0 = 0.3, dp = 0.3)

## truncated beta design priors: alternative, null, and an interior interval
expect_matches_reference("truncated beta H1", n1 = 8, n2 = 20, k = 1/6,
                         kf = 2, p0 = 0.4, a = 1, b = 1, da = 2, db = 3,
                         dl = 0.4, du = 1)
expect_matches_reference("truncated beta H0", n1 = 8, n2 = 20, k = 1/6,
                         kf = 2, p0 = 0.4, a = 1, b = 1, da = 2, db = 3,
                         dl = 0, du = 0.4)
expect_matches_reference("interior interval", n1 = 8, n2 = 20, k = 1/6,
                         kf = 2, p0 = 0.4, a = 1, b = 1, da = 2, db = 3,
                         dl = 0.45, du = 0.8)

## non-default analysis prior and a futility threshold below the efficacy one
expect_matches_reference("informative analysis prior, kf < k", n1 = 9, n2 = 15,
                         k = 1/2, kf = 1/3, p0 = 0.5, a = 2, b = 1.5,
                         dp = 0.65)

## calls of the examples in help("pbinbf01seq")
expect_matches_reference("example, truncated beta", n1 = 10, n2 = 20,
                         k = 1/10, kf = 3, p0 = 0.2, a = 1, b = 1, dl = 0.2,
                         du = 1, type = "direction")
expect_matches_reference("example, point", n1 = 10, n2 = 20, k = 1/10, kf = 3,
                         p0 = 0.2, a = 1, b = 1, dp = 0.2, type = "direction")

## Critical counts follow the integer grid of pbinbf01(): a threshold that is
## the BF01 of x = 6 successes in 10 observations includes x = 6 exactly
## (reference by exhaustive enumeration: 0.3408203125).
k <- binbf01(x = 6, n = 10, type = "direction")
expect_equal(
    pbinbf01seq(n1 = 2, n2 = 10, k = k, kf = 3, dp = 0.5), 0.3408203125,
    tolerance = 1e-12,
    info = "pbinbf01seq should use the integer grid of critical counts"
)

## Final efficacy that is unattainable, and interim continuation that is
## impossible, have probability zero.
expect_equal(
    pbinbf01seq(n1 = 2, n2 = 3, k = 0.01, kf = 3, dp = 0.5), 0,
    info = "unattainable final efficacy should have probability 0"
)
expect_equal(
    pbinbf01seq(n1 = 2, n2 = 10, k = 0.05, kf = 0.1, p0 = 0.5, dp = 0.5), 0,
    info = "impossible interim continuation should have probability 0"
)

## Tiny predictive masses of the truncated beta design prior stay on the log
## scale. Reference values from the enumeration (log-scale masses).
expect_equal(
    pbinbf01seq(n1 = 10, n2 = 20, p0 = 0.2, k = 1/10, kf = 10, da = 1,
                db = 100, dl = 0.5, du = 1),
    0.87741631847, tolerance = 1e-10,
    info = "pbinbf01seq should handle a tiny upper-tail design mass"
)
expect_equal(
    pbinbf01seq(n1 = 10, n2 = 40, k = 1/10, kf = 3, p0 = 0.5, da = 2000,
                db = 1, dl = 0, du = 0.5),
    0.0754139872159, tolerance = 1e-10,
    info = "pbinbf01seq should handle a tiny lower-tail design mass"
)

## A futility threshold so large that the study never stops at the interim
## stage reduces to the fixed design with n2 observations.
expect_equal(
    pbinbf01seq(n1 = 10, n2 = 25, k = 1/10, kf = 1e6, p0 = 0.3, dp = 0.5),
    pbinbf01(k = 1/10, n = 25, p0 = 0.3, type = "direction", dp = 0.5),
    tolerance = 1e-12,
    info = "pbinbf01seq with a huge kf should equal pbinbf01"
)
