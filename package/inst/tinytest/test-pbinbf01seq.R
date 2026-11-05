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

## Fractional numbers of observations are rounded up once and the rounded
## sizes are used everywhere, including the expected sample sizes.
expect_equal(
    pbinbf01seq(n1 = 10.1, n2 = 20.1, p0 = 0.2, dp = 0.4),
    pbinbf01seq(n1 = 11, n2 = 21, p0 = 0.2, dp = 0.4),
    tolerance = 1e-14,
    info = "pbinbf01seq should round n1 and n2 up"
)
res <- powerbinbf01seq(n1 = 10.1, n2 = 20.1, p0 = 0.2, dp = 0.4)
expect_identical(c(res$n1, res$n2), c(11, 21),
                 info = "powerbinbf01seq should store the rounded sizes")
expect_equal(res$nexp, res$pintfut*11 + (1 - res$pintfut)*21,
             tolerance = 1e-14,
             info = "nexp should use the rounded sizes")
expect_equal(res$nexp0, res$pintfut0*11 + (1 - res$pintfut0)*21,
             tolerance = 1e-14,
             info = "nexp0 should use the rounded sizes")

## The interim stage has to come before the final stage after rounding.
expect_error(
    pbinbf01seq(n1 = 10.2, n2 = 10.8, dp = 0.5), pattern = "n1 < n2",
    info = "pbinbf01seq should require n1 < n2 after rounding"
)
expect_error(
    powerbinbf01seq(n1 = 10.2, n2 = 10.8, dp = 0.5), pattern = "n1 < n2",
    info = "powerbinbf01seq should require n1 < n2 after rounding"
)

## Argument checks give informative errors.
expect_error(pbinbf01seq(n1 = 5, n2 = 10, kf = -1, dp = 0.5), pattern = "kf",
             info = "negative kf should be rejected")
expect_error(pbinbf01seq(n1 = 5, n2 = 10, kf = NA, dp = 0.5), pattern = "kf",
             info = "missing kf should be rejected")
expect_error(pbinbf01seq(n1 = 5, n2 = 10, k = -1, dp = 0.5), pattern = "0 < k",
             info = "negative k should be rejected")
expect_error(pbinbf01seq(n1 = 5, n2 = 10, p0 = 1, dp = 0.5), pattern = "p0",
             info = "p0 = 1 should be rejected")
expect_error(pbinbf01seq(n1 = 5, n2 = 10, a = 0, dp = 0.5), pattern = "0 < a",
             info = "a = 0 should be rejected")
expect_error(pbinbf01seq(n1 = 5, n2 = 10, b = NA, dp = 0.5), pattern = "b",
             info = "missing b should be rejected")
expect_error(powerbinbf01seq(n1 = 5, n2 = 10, kf = NA, dp = 0.5),
             pattern = "kf",
             info = "powerbinbf01seq should check kf")
expect_error(powerbinbf01seq(n1 = 5, n2 = 10, type = "point", dp = 0.5),
             pattern = "direction",
             info = "powerbinbf01seq should only accept type = \"direction\"")

## powerbinbf01seq: the expected sample sizes weight the interim and final
## sample sizes with the probabilities of stopping for futility, the names of
## the returned object are unique, and the probability to continue and reach
## efficacy cannot exceed the probability not to stop for futility.
for (dp in c(0.4, NA)) {
    res <- powerbinbf01seq(n1 = 12, n2 = 30, k = 1/10, kf = 3, p0 = 0.2,
                           dp = dp, dl = 0.2, du = 1)
    label <- if (is.na(dp)) "truncated beta design prior" else "point design prior"
    expect_true(inherits(res, "power.bftestseq"),
                info = paste(label, "should return a power.bftestseq object"))
    expect_true(!anyDuplicated(names(res)),
                info = paste(label, "should return uniquely named elements"))
    expect_equal(res$nexp, res$pintfut*12 + (1 - res$pintfut)*30,
                 tolerance = 1e-14, info = paste(label, "nexp identity"))
    expect_equal(res$nexp0, res$pintfut0*12 + (1 - res$pintfut0)*30,
                 tolerance = 1e-14, info = paste(label, "nexp0 identity"))
    expect_equal(res$pintfut,
                 pbinbf01(k = 3, n = 12, p0 = 0.2, type = "direction", dp = dp,
                          dl = 0.2, du = 1, lower.tail = FALSE),
                 tolerance = 1e-14,
                 info = paste(label, "pintfut should be the interim futility probability"))
    expect_equal(res$pfineff,
                 pbinbf01seq(n1 = 12, n2 = 30, k = 1/10, kf = 3, p0 = 0.2,
                             dp = dp, dl = 0.2, du = 1),
                 tolerance = 1e-14,
                 info = paste(label, "pfineff should be the two-stage probability"))
    expect_true(res$pfineff <= 1 - res$pintfut + 1e-12 &&
                res$pfineff0 <= 1 - res$pintfut0 + 1e-12,
                info = paste(label, "pfineff cannot exceed 1 - pintfut"))
}

## At kf = 1 and p0 = 0.5 with a flat analysis prior, BF01 = 1 exactly for
## x = n/2, which the floating-point predicate classifies as BF01 > 1 in
## pbinbf01() and pbinbf01seq() alike, so that the two stay coherent.
res <- powerbinbf01seq(n1 = 12, n2 = 24, k = 1/10, kf = 1, p0 = 0.5, dp = 0.7)
expect_true(res$pfineff <= 1 - res$pintfut + 1e-12,
            info = "interim ties should be classified alike in both functions")
expect_equal(res$pfineff,
             ref_pbinbf01seq(n1 = 12, n2 = 24, k = 1/10, kf = 1, p0 = 0.5,
                             dp = 0.7),
             tolerance = 1e-12,
             info = "pfineff at the interim tie should match the enumeration")

## print method: the object is returned unchanged, and the output uses the
## thresholds in the table and in the labels.
for (dp in c(0.4, NA)) {
    res <- powerbinbf01seq(n1 = 12, n2 = 30, k = 1/10, kf = 3, p0 = 0.2,
                           dp = dp, dl = 0.2, du = 1)
    out <- capture.output(printed <- print(res))
    expect_identical(printed, res,
                     info = "print should return the object unchanged")
    expect_true(is.numeric(printed$k) && is.numeric(printed$kf),
                info = "print should not turn the thresholds into strings")
    condition <- if (is.na(dp)) "design prior" else "p = 0.4"
    expect_true(
        any(grepl("BF efficacy threshold k = 1/10", out, fixed = TRUE)) &&
        any(grepl("BF futility threshold kf = 3", out, fixed = TRUE)),
        info = "print should show the thresholds"
    )
    expect_true(
        any(grepl(paste0("Pr(BF01 > 3 at interim | ", condition, ")"), out,
                  fixed = TRUE)) &&
        any(grepl("Pr(BF01 <= 1/10 at final | p = 0.2)", out, fixed = TRUE)),
        info = "print labels should state the events of the design"
    )
    expect_true(
        any(grepl("NOTE: BF oriented in favor of H0 (BF01 < 1 indicates evidence for H1 over H0)",
                  out, fixed = TRUE)),
        info = "print should show the note on a single line"
    )
    expect_false(any(grepl("Plot power", out, fixed = TRUE)),
                 info = "print should not suggest a power plot")
}
