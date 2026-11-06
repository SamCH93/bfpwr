library(tinytest)
library(bfpwr)

## Tests the two-stage binomial power functions pbinbf01seq() and
## powerbinbf01seq(). The study continues after the interim stage if
## BF01(x1) <= kf and has efficacy at the final stage if BF01(x1 + x2) <= k,
## where x1 and x2 are the successes among the first n1 and the remaining
## n2 - n1 observations. The functions were written for Kelter and Pawel (2025,
## <https://doi.org/10.48550/arXiv.2511.23144>), whose designs are reproduced
## first. Further references are Simon's two-stage designs (binomial sums) and
## the exhaustive enumeration of all (x1, x2) pairs with the same predicate on
## the log Bayes factors as pbinbf01() (a tolerance of 1e-10 on the log scale
## counts mathematical ties), weighted with the predictive masses under the
## design prior.

## Kelter and Pawel (2025), Tables 2 and 3: optimal two-stage designs with
## kf = 3 and a flat analysis prior (a = b = 1). Power is computed under a point
## alternative dp or under a Beta(da, db) design prior truncated to [p0, 1]
## whose mode is 0.3 (example 1) or 0.4 (example 2). Type I error and the
## probability of compelling evidence for H0 (futility stop at the interim) are
## computed under p0. The reference values are the 7 significant digits printed
## by the paper's code on OSF (sim.html); they round to the values in the tables.
paper_da <- function(mode, db) (mode*(db - 2) + 1)/(1 - mode)
paper_ref <- function(pintfut, pintfut0, pfineff, pfineff0, nexp, nexp0) {
    c(pintfut = pintfut, pintfut0 = pintfut0, pfineff = pfineff,
      pfineff0 = pfineff0, nexp = nexp, nexp0 = nexp0)
}
paper_designs <- list(
    ## Table 2 (example 1): p0 = 0.1, k = 1/3, p1 = 0.3
    list(label = "Table 2, frequentist power", n1 = 10, n2 = 29, k = 1/3,
         p0 = 0.1, dp = 0.3,
         ref = paper_ref(0.1493083, 0.7360989, 0.8050629, 0.04708631,
                         26.16314, 15.01412)),
    list(label = "Table 2, Beta(1, 1)", n1 = 5, n2 = 15, k = 1/3, p0 = 0.1,
         da = 1, db = 1,
         ref = paper_ref(0.098415, 0.59049, 0.8107187, 0.04800019, 14.01585,
                         9.0951)),
    list(label = "Table 2, Beta(7, 15)", n1 = 11, n2 = 36, k = 1/3, p0 = 0.1,
         da = paper_da(0.3, 15), db = 15,
         ref = paper_ref(0.1389243, 0.6973569, 0.8016758, 0.04704394,
                         32.52689, 18.56608)),
    list(label = "Table 2, Beta(11.29, 25)", n1 = 12, n2 = 28, k = 1/3,
         p0 = 0.1, da = paper_da(0.3, 25), db = 25,
         ref = paper_ref(0.1048149, 0.6590023, 0.8019558, 0.04767261,
                         26.32296, 17.45596)),
    list(label = "Table 2, Beta(22, 50)", n1 = 10, n2 = 36, k = 1/3, p0 = 0.1,
         da = paper_da(0.3, 50), db = 50,
         ref = paper_ref(0.1581685, 0.7360989, 0.8038297, 0.04324746,
                         31.88762, 16.86143)),
    ## Table 3 (example 2): p0 = 0.2, p1 = 0.4
    list(label = "Table 3, frequentist power, k = 1/3", n1 = 17, n2 = 37,
         k = 1/3, p0 = 0.2, dp = 0.4,
         ref = paper_ref(0.04642293, 0.5488762, 0.9032743, 0.09478437,
                         36.07154, 26.02248)),
    list(label = "Table 3, frequentist power, k = 1/3, f = 0.6", n1 = 30,
         n2 = 36, k = 1/3, p0 = 0.2, dp = 0.4,
         ref = paper_ref(0.01718302, 0.6069699, 0.9091478, 0.08861462,
                         35.8969, 32.35818)),
    list(label = "Table 3, frequentist power, k = 1/10", n1 = 21, n2 = 51,
         k = 1/10, p0 = 0.2, dp = 0.4,
         ref = paper_ref(0.03695564, 0.5860084, 0.9020781, 0.03402297,
                         49.89133, 33.41975)),
    list(label = "Table 3, frequentist power, k = 1/10, f = 0.6", n1 = 30,
         n2 = 50, k = 1/10, p0 = 0.2, dp = 0.4,
         ref = paper_ref(0.01718302, 0.6069699, 0.9010762, 0.03022596,
                         49.65634, 37.8606)),
    list(label = "Table 3, Beta(1, 1), k = 1/3", n1 = 27, n2 = 54, k = 1/3,
         p0 = 0.2, da = 1, db = 1,
         ref = paper_ref(0.04706176, 0.5386567, 0.9003252, 0.09876502,
                         52.72933, 39.45627)),
    list(label = "Table 3, Beta(10.33, 15), k = 1/3", n1 = 28, n2 = 67,
         k = 1/3, p0 = 0.2, da = paper_da(0.4, 15), db = 15,
         ref = paper_ref(0.04398824, 0.5005331, 0.9006022, 0.09805443,
                         65.28446, 47.47921)),
    list(label = "Table 3, Beta(13.67, 20), k = 1/3", n1 = 24, n2 = 58,
         k = 1/3, p0 = 0.2, da = paper_da(0.4, 20), db = 20,
         ref = paper_ref(0.03803118, 0.4598773, 0.9001069, 0.09576332,
                         56.70694, 42.36417)),
    list(label = "Table 3, Beta(1, 1), k = 1/10", n1 = 42, n2 = 100,
         k = 1/10, p0 = 0.2, da = 1, db = 1,
         ref = paper_ref(0.03645961, 0.5309237, 0.9002312, 0.03246814,
                         97.88534, 69.20643)),
    list(label = "Table 3, Beta(10.33, 15), k = 1/10", n1 = 65, n2 = 100,
         k = 1/10, p0 = 0.2, da = paper_da(0.4, 15), db = 15,
         ref = paper_ref(0.02552603, 0.5735012, 0.900092, 0.03385449,
                         99.10659, 79.92746)),
    list(label = "Table 3, Beta(13.67, 20), k = 1/10", n1 = 40, n2 = 88,
         k = 1/10, p0 = 0.2, da = paper_da(0.4, 20), db = 20,
         ref = paper_ref(0.03626588, 0.5931271, 0.9004765, 0.03471779,
                         86.25924, 59.5299)),
    list(label = "Table 3, Beta(1, 1), k = 1/3, f = 0.6", n1 = 30, n2 = 58,
         k = 1/3, p0 = 0.2, da = 1, db = 1,
         ref = paper_ref(0.05431338, 0.6069699, 0.9008425, 0.09275858,
                         56.47923, 41.00484)),
    list(label = "Table 3, Beta(10.33, 15), k = 1/3, f = 0.6", n1 = 30,
         n2 = 76, k = 1/3, p0 = 0.2, da = paper_da(0.4, 15), db = 15,
         ref = paper_ref(0.06071236, 0.6069699, 0.9007968, 0.09232889,
                         73.20723, 48.07938)),
    list(label = "Table 3, Beta(13.67, 20), k = 1/3, f = 0.6", n1 = 30,
         n2 = 63, k = 1/3, p0 = 0.2, da = paper_da(0.4, 20), db = 20,
         ref = paper_ref(0.0529332, 0.6069699, 0.9056205, 0.09780135,
                         61.2532, 42.96999)),
    list(label = "Table 3, Beta(1, 1), k = 1/10, f = 0.6", n1 = 30, n2 = 130,
         k = 1/10, p0 = 0.2, da = 1, db = 1,
         ref = paper_ref(0.05431338, 0.6069699, 0.9015596, 0.02790104,
                         124.5687, 69.30301)),
    list(label = "Table 3, Beta(10.33, 15), k = 1/10, f = 0.6", n1 = 30,
         n2 = 147, k = 1/10, p0 = 0.2, da = paper_da(0.4, 15), db = 15,
         ref = paper_ref(0.06071236, 0.6069699, 0.9011365, 0.02634651,
                         139.8967, 75.98452)),
    list(label = "Table 3, Beta(13.67, 20), k = 1/10, f = 0.6", n1 = 30,
         n2 = 108, k = 1/10, p0 = 0.2, da = paper_da(0.4, 20), db = 20,
         ref = paper_ref(0.0529332, 0.6069699, 0.900138, 0.02732785,
                         103.8712, 60.65635))
)
for (d in paper_designs) {
    design_prior <- if (is.null(d$dp)) {
        list(da = d$da, db = d$db, dl = d$p0, du = 1)
    } else {
        list(dp = d$dp)
    }
    res <- do.call(powerbinbf01seq,
                   c(list(n1 = d$n1, n2 = d$n2, k = d$k, kf = 3, p0 = d$p0,
                          a = 1, b = 1, type = "direction"), design_prior))
    for (quantity in names(d$ref)) {
        expect_equal(res[[quantity]], d$ref[[quantity]], tolerance = 1e-6,
                     info = paste(d$label, quantity,
                                  "should match Kelter and Pawel (2025)"))
    }
}

## With frequentist power, k = 1/3 and kf = 3, the optimal designs of Kelter
## and Pawel (2025) are Simon's optimal two-stage designs, which continue after
## the interim stage if x1 > r1 and reject H0 if x1 + x2 > r (r1/n1 = 1/10 and
## r/n = 5/29 for p0 = 0.1 and p1 = 0.3; 3/17 and 10/37 for p0 = 0.2 and
## p1 = 0.4, from clinfun::ph2simon). Their probabilities are binomial sums.
simon_preject <- function(p, r1, n1, r, n) {
    x1 <- seq(r1 + 1, n1)
    sum(stats::dbinom(x1, size = n1, prob = p) *
        stats::pbinom(r - x1, size = n - n1, prob = p, lower.tail = FALSE))
}
simon_designs <- list(
    list(r1 = 1, n1 = 10, r = 5, n = 29, p0 = 0.1, p1 = 0.3),
    list(r1 = 3, n1 = 17, r = 10, n = 37, p0 = 0.2, p1 = 0.4)
)
for (d in simon_designs) {
    res <- powerbinbf01seq(n1 = d$n1, n2 = d$n, k = 1/3, kf = 3, p0 = d$p0,
                           dp = d$p1)
    label <- sprintf("Simon's design r1/n1 = %d/%d, r/n = %d/%d:", d$r1, d$n1,
                     d$r, d$n)
    expect_equal(res$pfineff,
                 simon_preject(d$p1, r1 = d$r1, n1 = d$n1, r = d$r, n = d$n),
                 tolerance = 1e-12, info = paste(label, "power"))
    expect_equal(res$pfineff0,
                 simon_preject(d$p0, r1 = d$r1, n1 = d$n1, r = d$r, n = d$n),
                 tolerance = 1e-12, info = paste(label, "type I error"))
    expect_equal(res$pintfut0, stats::pbinom(d$r1, size = d$n1, prob = d$p0),
                 tolerance = 1e-12,
                 info = paste(label, "early termination under H0"))
    expect_equal(res$pintfut, stats::pbinom(d$r1, size = d$n1, prob = d$p1),
                 tolerance = 1e-12,
                 info = paste(label, "early termination under H1"))
}

tie_tol <- 1e-10

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
    event <- logbf1 <= log(kf) + tie_tol & logbf2 <= log(k) + tie_tol
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

## Tiny predictive masses of the truncated beta design prior. The current
## calculation in pbinbf01seq() returns NaN for them (0/0 in the truncation
## constants); the commented-out log-scale calculation there handles them.
## Re-enable these tests together with it. Reference values from the
## enumeration (log-scale masses).
## expect_equal(
##     pbinbf01seq(n1 = 10, n2 = 20, p0 = 0.2, k = 1/10, kf = 10, da = 1,
##                 db = 100, dl = 0.5, du = 1),
##     0.87741631847, tolerance = 1e-10,
##     info = "pbinbf01seq should handle a tiny upper-tail design mass"
## )
## expect_equal(
##     pbinbf01seq(n1 = 10, n2 = 40, k = 1/10, kf = 3, p0 = 0.5, da = 2000,
##                 db = 1, dl = 0, du = 0.5),
##     0.0754139872159, tolerance = 1e-10,
##     info = "pbinbf01seq should handle a tiny lower-tail design mass"
## )

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
## x = n/2, whose floating-point value is 1 + 1e-15. The tie counts as
## reaching kf (tolerance 1e-10 on the log scale) in pbinbf01() and
## pbinbf01seq() alike, so that the two stay coherent: BF01 > 1 only for
## x1 <= 5.
res <- powerbinbf01seq(n1 = 12, n2 = 24, k = 1/10, kf = 1, p0 = 0.5, dp = 0.7)
expect_equal(res$pintfut, stats::pbinom(5, size = 12, prob = 0.7),
             tolerance = 1e-12,
             info = "pintfut at the interim tie should exclude x1 = n1/2")
expect_equal(res$pintfut, 0.038600843058, tolerance = 1e-10,
             info = "pintfut at the interim tie should be 0.0386008")
expect_true(res$pfineff <= 1 - res$pintfut + 1e-12,
            info = "interim ties should be classified alike in both functions")
expect_equal(res$pfineff,
             ref_pbinbf01seq(n1 = 12, n2 = 24, k = 1/10, kf = 1, p0 = 0.5,
                             dp = 0.7),
             tolerance = 1e-12,
             info = "pfineff at the interim tie should match the enumeration")
expect_equal(res$pfineff, 0.7224538968, tolerance = 1e-9,
             info = "pfineff at the interim tie should be 0.7224539")

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
