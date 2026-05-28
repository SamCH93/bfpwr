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
    expect_true(is.finite(seqres$cumpH1) && seqres$cumpH1 > 0 &&
                    seqres$cumpH1 < 0.01,
                info = paste(alt, "ptbf01seq H1 probability should be finite"))
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
