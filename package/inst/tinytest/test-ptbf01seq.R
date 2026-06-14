library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "sequential t BF boundary checks are extended"
    ))
}

## Tests ptbf01seq one-stage equivalence to ptbf01 and impossible H0 boundaries.
## Related manuscript source: t BF section in paper/bfssd.Rnw 1481-1530 and the
## BFGSD appendix JZS sequence; these specific fixtures are package regressions.

## One-stage sequential designs should agree with the non-sequential t-test
## power calculation. This also exercises the internal tcrit() root search.
regression_n <- 9596.363636363636
common_args <- list(n = regression_n, plocation = 0, pscale = 0.707,
                    pdf = 1, type = "two.sample", dpm = 0, dpsd = 0)

for (alt in c("greater", "less", "two.sided")) {
    seqres <- suppressWarnings(
        do.call(ptbf01seq, c(list(k1 = 1/10, k0 = 10, alternative = alt),
                             common_args))
    )
    pH1 <- suppressWarnings(
        do.call(ptbf01, c(list(k = 1/10, alternative = alt), common_args))
    )
    pH0 <- suppressWarnings(
        do.call(ptbf01, c(list(k = 10, alternative = alt, lower.tail = FALSE),
                          common_args))
    )
    expect_true(is.finite(seqres$cumpH1) && seqres$cumpH1 >= 0 &&
                    seqres$cumpH1 < 0.01,
                info = paste(alt, "ptbf01seq H1 probability should be finite"))
    if (alt == "two.sided") {
        expect_true(seqres$cumpH1 > 0,
                    info = "two-sided ptbf01seq H1 probability should be positive")
    }
    expect_true(is.finite(seqres$cumpH0) && seqres$cumpH0 > 0.9 &&
                    seqres$cumpH0 < 1,
                info = paste(alt, "ptbf01seq H0 probability should be finite"))
    expect_true(abs(seqres$cumpH1 - pH1) < 1e-6,
                info = paste(alt, "one-stage H1 probability should match ptbf01"))
    expect_true(abs(seqres$cumpH0 - pH0) < 5e-4,
                info = paste(alt, "one-stage H0 probability should match ptbf01"))
}

limit_warning <- NULL
missing_h0 <- try(
    withCallingHandlers(
        ptbf01seq(k1 = 1/10, k0 = 10, n = c(5, 10), plocation = 0,
                  pscale = 0.707, pdf = 1, type = "two.sample",
                  alternative = "greater", dpm = 0, dpsd = 0),
        warning = function(w) {
            limit_warning <<- conditionMessage(w)
            invokeRestart("muffleWarning")
        }
    ),
    silent = TRUE
)
expect_false(inherits(missing_h0, "try-error"),
             info = "ptbf01seq should handle stages where H0 boundary is impossible")
expect_true(all(missing_h0$cumpH0 < 1e-12),
            info = "impossible one-sided H0 boundaries should have zero stop probability")
expect_true(grepl("Adaptive t critical-value search reached", limit_warning,
                  fixed = TRUE),
            info = "ptbf01seq should warn when adaptive tcrit search reaches its limit")
expect_true(grepl("marginal tail probability <= 0.001", limit_warning,
                  fixed = TRUE),
            info = "ptbf01seq adaptive-limit warning should report the per-boundary tail-eps cutoff")
expect_equal(missing_h0$tail.eps, 1e-3,
             info = "ptbf01seq should store the sequential tail.eps control")

missing_h1_continuation <- bfpwr:::genregions1(
    zcrit0 = c(NaN, 0),
    zcrit1 = c(NaN, 2),
    direction = "positive"
)
expect_equal(missing_h1_continuation$H1[[2]][[1]][, 1], c(-Inf, Inf),
             info = paste("missing earlier H1 boundaries should become",
                          "unbounded positive-direction continuation"))
expect_equal(missing_h1_continuation$H0[[2]][[1]][, 1], c(-Inf, Inf),
             info = paste("missing earlier H1 boundaries should not zero",
                          "later H0 continuation"))

missing_h1_continuation_less <- bfpwr:::genregions1(
    zcrit0 = c(NaN, 0),
    zcrit1 = c(NaN, -2),
    direction = "negative"
)
expect_equal(missing_h1_continuation_less$H1[[2]][[1]][, 1], c(-Inf, Inf),
             info = paste("missing earlier H1 boundaries should become",
                           "unbounded negative-direction continuation"))

expect_error(
    bfpwr:::genregions1(zcrit0 = c(NA_real_, 0), zcrit1 = c(2, 2)),
    "Critical values cannot contain NA",
    info = "one-sided region generation should reject plain NA boundaries"
)
expect_error(
    bfpwr:::genregions2(
        zcrit0 = matrix(c(NA_real_, 0, 1, 1), nrow = 2),
        zcrit1 = matrix(c(-2, -2, 2, 2), nrow = 2)
    ),
    "Critical values cannot contain NA",
    info = "two-sided region generation should reject plain NA boundaries"
)

empty_region_probability <- bfpwr:::.bfseq_intstage_sum(
    stageregions = list(matrix(c(NaN, NaN), nrow = 2)),
    mean = 0,
    sigma = matrix(1)
)
expect_equal(empty_region_probability, 0,
             info = "generated NaN regions should remain empty stopping events")
na_bound_probability <- try(
    bfpwr:::.bfseq_intstage_sum(
        stageregions = list(matrix(c(NA_real_, 1), nrow = 2)),
        mean = 0,
        sigma = matrix(1)
    ),
    silent = TRUE
)
expect_true(
    inherits(na_bound_probability, "try-error") &&
        grepl("NA bounds",
              conditionMessage(attr(na_bound_probability, "condition")),
              fixed = TRUE),
    info = "non-empty sequential regions with NA bounds should error"
)

custom_tail <- suppressWarnings(
    ptbf01seq(k1 = 1/10, k0 = 10, n = 20, plocation = 0,
              pscale = 0.707, pdf = 1, type = "two.sample",
              alternative = "greater", dpm = 0.5, dpsd = 0.1,
              tail.eps = 1e-2)
)
expect_equal(custom_tail$tail.eps, 1e-2,
             info = "ptbf01seq should preserve custom tail.eps values")

custom_limit_warning <- NULL
suppressWarnings(
    withCallingHandlers(
        ptbf01seq(k1 = 1/10, k0 = 10, n = c(5, 10), plocation = 0,
                  pscale = 0.707, pdf = 1, type = "two.sample",
                  alternative = "greater", dpm = 0, dpsd = 0,
                  tail.eps = 1e-2),
        warning = function(w) {
            custom_limit_warning <<- conditionMessage(w)
        }
    )
)
expect_true(grepl("marginal tail probability <= 0.01",
                  custom_limit_warning, fixed = TRUE),
            info = "ptbf01seq adaptive-limit warning should use custom tail.eps")

finite_zcrit0 <- matrix(rep(c(-1, 1), 10), nrow = 2)
region_count <- bfpwr:::.count_strict_two_sided_regions(finite_zcrit0)
expect_equal(region_count$total, 3069,
             info = "strict two-sided region count should match exact branching")
expect_equal(region_count$firstH0, 1,
             info = "region count should identify first finite H0 boundary")

slow_warning <- NULL
slow_exact <- try(
    withCallingHandlers(
        ptbf01seq(k1 = 1/10, k0 = 3, n = seq(40, 130, 10),
                  plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                  type = "two.sample", alternative = "two.sided",
                  dpm = 0.5, dpsd = 0.1, strict = TRUE),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("strict = TRUE with two-sided sequential t testing",
                      msg, fixed = TRUE)) {
                slow_warning <<- msg
                stop("caught expected strict two-sided warning")
            }
        }
    ),
    silent = TRUE
)
expect_true(inherits(slow_exact, "try-error"),
            info = "test should abort as soon as the slow exact warning appears")
expect_true(grepl("integrate 3,069 regions", slow_warning, fixed = TRUE),
            info = "ptbf01seq should warn immediately before slow exact integration")

explicit_trange <- ptbf01seq(k1 = 1/10, k0 = 10, n = 100, plocation = 0,
                             pscale = 0.707, pdf = 1, type = "two.sample",
                             alternative = "greater", dpm = 0.5, dpsd = 0.1,
                             trange = c(-2, 6))
expect_true(is.finite(explicit_trange$cumpH1) &&
                is.finite(explicit_trange$cumpH0),
            info = "ptbf01seq should accept explicit t-statistic trange")
expect_equal(explicit_trange$trange, c(-2, 6),
             info = "ptbf01seq should store the explicit t-statistic trange")
expect_equal(explicit_trange$tail.eps, 1e-3,
             info = "ptbf01seq should store tail.eps even when explicit trange is used")

narrow_trange <- try(
    ptbf01seq(k1 = 1/10, k0 = 10, n = 100, plocation = 0,
              pscale = 0.707, pdf = 1, type = "two.sample",
              alternative = "greater", dpm = 0.5, dpsd = 0.1,
              trange = c(-1, 1)),
    silent = TRUE
)
expect_true(
    inherits(narrow_trange, "try-error") &&
        grepl("Failed to compute H1 sequential t stopping boundary",
              conditionMessage(attr(narrow_trange, "condition")),
              fixed = TRUE),
    info = "ptbf01seq should not convert failed numeric H1 boundary searches to empty regions"
)

narrow_h0_trange <- try(
    ptbf01seq(k1 = 1/10, k0 = 2, n = 100, plocation = 0,
              pscale = 0.707, pdf = 1, type = "two.sample",
              alternative = "two.sided", dpm = 0.5, dpsd = 0.1,
              trange = c(2.5, 3)),
    silent = TRUE
)
expect_true(
    inherits(narrow_h0_trange, "try-error") &&
        grepl("Failed to compute H0 sequential t stopping boundary",
              conditionMessage(attr(narrow_h0_trange, "condition")),
              fixed = TRUE),
    info = "ptbf01seq should not treat numeric-range H0 misses as impossible boundaries"
)

impossible_h0_numeric_trange <- suppressWarnings(
    ptbf01seq(k1 = 1/30, k0 = 1e9, n = c(10, 20), plocation = 0,
              pscale = 0.707, pdf = 1, type = "two.sample",
              alternative = "two.sided", dpm = 0.5, dpsd = 0.1,
              strict = FALSE, trange = c(-10, 10))
)
expect_true(all(impossible_h0_numeric_trange$cumpH0 == 0),
            info = "ptbf01seq should allow proven impossible H0 boundaries with numeric trange")
expect_true(all(is.finite(impossible_h0_numeric_trange$cumpH1)),
            info = "proven impossible numeric-trange H0 boundaries should not block H1 probabilities")

bad_tail <- try(
    ptbf01seq(k1 = 1/10, k0 = 10, n = 100, plocation = 0,
              pscale = 0.707, pdf = 1, type = "two.sample",
              alternative = "greater", dpm = 0.5, dpsd = 0.1,
              tail.eps = 0.5),
    silent = TRUE
)
expect_true(inherits(bad_tail, "try-error"),
            info = "ptbf01seq should validate tail.eps")

old_drange <- try(
    ptbf01seq(k1 = 1/10, k0 = 10, n = 100, plocation = 0,
              pscale = 0.707, pdf = 1, type = "two.sample",
              alternative = "greater", dpm = 0.5, dpsd = 0.1,
              drange = c(-2, 6)),
    silent = TRUE
)
expect_true(
    inherits(old_drange, "try-error") &&
        grepl("renamed to 'trange'",
              conditionMessage(attr(old_drange, "condition")), fixed = TRUE),
    info = "ptbf01seq should no longer accept drange"
)

## ## do not run these tests for the moment, because they are there to verify
## ## the power with simulation which takes a long time to run

## library(tinytest)
## library(bfpwr)

## ## verify that computed power equal to simulated power
## grid <- expand.grid(plocation = 0,
##                     pscale = 0.71,
##                     pdf = 1,
##                     dpm = c(0.5, 0),
##                     dpsd = c(0, 0.1),
##                     alternative = c("two.sided", "greater", "less"),
##                     type = c("one.sample", "two.sample"),
##                     stringsAsFactors = FALSE)
## n <- seq(40, 100, 10)
## k0 <- 6
## k1 <- 1/30
## set.seed(42)
## nsim <- 10000
## pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
## pH1sim <- numeric(nrow(grid))
## pH0sim <- numeric(nrow(grid))
## pIncsim <- numeric(nrow(grid))
## pH1 <- numeric(nrow(grid))
## pH0 <- numeric(nrow(grid))
## pInc <- numeric(nrow(grid))
## for (i in seq(1, nrow(grid))) {
##     setTxtProgressBar(pb, i)
##     test <- paste(i, paste(colnames(grid), grid[i,], sep = "=", collapse = ", "))
##     ## simulate
##     results <- replicate(n = nsim, expr = {
##         smd <- rnorm(n = 1, mean = grid$dpm[i], sd = grid$dpsd[i])
##         if (grid$type[i] == "one.sample") {
##             y <- rnorm(n = max(n), mean = smd, sd = 1)
##             t <- sapply(seq_along(n), FUN = function(j) {
##                 ttest <- t.test(y[1:n[j]], alternative = "two.sided")$statistic
##             })
##         } else {
##             y1 <- rnorm(n = max(n), mean = 0, sd = 1)
##             y2 <- rnorm(n = max(n), mean = smd, sd = 1)
##             t <- sapply(seq_along(n), FUN = function(j) {
##                 ttest <- t.test(y2[1:n[j]], y1[1:n[j]], var.equal = TRUE,
##                                 alternative = "two.sided")$statistic
##             })
##         }
##         bf <- sapply(seq_along(n), FUN = function(j) {
##             tbf01(t = t[j], n = n[j], plocation = grid$plocation[i],
##                   pscale = grid$pscale[i], pdf = grid$pdf[i],
##                   type = grid$type[i], alternative = grid$alternative[i])
##         })
##         result <- "inconclusive"
##         for (j in seq_along(bf)) {
##             if (bf[j] >= k0) {
##                 result <- "H0"
##                 break
##             }
##             if (bf[j] <= k1) {
##                 result <- "H1"
##                 break
##             }
##         }
##         result
##     })
##     pH1sim[i] <- mean(results == "H1")
##     pH0sim[i] <- mean(results == "H0")
##     pIncsim[i] <- mean(results == "inconclusive")
##     res <- ptbf01seq(k1 = k1, k0 = k0, n1 = n, n2 = n,
##                      plocation = grid$plocation[i], pscale = grid$pscale[i],
##                      pdf = grid$pdf[i], type = grid$type[i], dpm = grid$dpm[i],
##                      dpsd = grid$dpsd[i], alternative = grid$alternative[i],
##                      strict = FALSE)
##     pH1[i] <- tail(res$cumpH1, n = 1)
##     pH0[i] <- tail(res$cumpH0, n = 1)
##     pInc[i] <- tail(res$cumpInc, n = 1)
##     cat("\n")
##     print(expect_equal(pH1[i], pH1sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pH0[i], pH0sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pInc[i], pIncsim[i], info = test, tolerance = 0.02))
## }


## for (i in seq(1, nrow(grid))) {
##     cat("\n")
##     print(expect_equal(pH1[i], pH1sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pH0[i], pH0sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pInc[i], pIncsim[i], info = test, tolerance = 0.02))
## }

## ## percentage scale
## round(cbind(pH1, pH1sim, pH0, pH0sim, pInc, pIncsim)*100, 2)

## ## percentage differences
## round(cbind(pH1 - pH1sim, pH0 - pH0sim, pInc - pIncsim)*100, 2)
