library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "t BF power-boundary checks are extended"
    ))
}

## Tests ptbf01 adaptive boundary selection, impossible-boundary handling,
## shifted-null recentering, and small tail probabilities. Manuscript source:
## tBF/power example in paper/bfssd.Rnw 1481-1530; edge cases are regressions.

## Regression test for adaptive root selection in one-sided t-test designs.
## This fractional n occurs in the default plot grid for nlim = c(10, 10000).
regression_n <- 9596.363636363636
common_args <- list(n = regression_n, null = 0, plocation = 0, pscale = 0.707,
                    pdf = 1, type = "two.sample", dpm = 0, dpsd = 0)

greater_adaptive <- suppressWarnings(
    do.call(ptbf01, c(list(k = 1/10, alternative = "greater"), common_args))
)
greater_positive_range <- suppressWarnings(
    do.call(ptbf01, c(list(k = 1/10, alternative = "greater", drange = c(0, 2)),
                      common_args))
)
expect_true(is.finite(greater_adaptive) && greater_adaptive < 0.01,
            info = "greater one-sided adaptive search should not select the lower root")
expect_true(abs(greater_adaptive - greater_positive_range) < 5e-6,
            info = "greater one-sided adaptive search should match positive-side search")

less_adaptive <- suppressWarnings(
    do.call(ptbf01, c(list(k = 1/10, alternative = "less"), common_args))
)
less_negative_range <- suppressWarnings(
    do.call(ptbf01, c(list(k = 1/10, alternative = "less", drange = c(-2, 0)),
                      common_args))
)
expect_true(is.finite(less_adaptive) && less_adaptive < 0.01,
            info = "less one-sided adaptive search should not select the upper root")
expect_true(abs(less_adaptive - less_negative_range) < 5e-6,
            info = "less one-sided adaptive search should match negative-side search")

h0_adaptive <- suppressWarnings(
    do.call(ptbf01, c(list(k = 10, alternative = "greater", lower.tail = FALSE),
                      common_args))
)
h0_positive_range <- suppressWarnings(
    do.call(ptbf01, c(list(k = 10, alternative = "greater", lower.tail = FALSE,
                           drange = c(0, 2)), common_args))
)
expect_true(is.finite(h0_adaptive) && h0_adaptive > 0.9 && h0_adaptive < 1,
            info = "greater one-sided H0 power should remain on the positive root")
expect_true(abs(h0_adaptive - h0_positive_range) < 5e-4,
            info = "greater one-sided H0 adaptive search should match positive-side search")

limit_warning <- NULL
impossible_h0 <- withCallingHandlers(
    ptbf01(k = 10, n = 5, plocation = 0, pscale = 0.707, pdf = 1,
           type = "two.sample", alternative = "greater",
           dpm = 0, dpsd = 0, lower.tail = FALSE),
    warning = function(w) {
        limit_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    }
)
expect_equal(impossible_h0, 0,
             info = "one-sided ptbf01 should return zero H0 probability when threshold is not reachable")
expect_true(grepl("Adaptive t power-boundary search reached", limit_warning,
                  fixed = TRUE),
            info = "one-sided ptbf01 should warn when adaptive boundary search reaches its limit")

## Regression cases for one-sided H0 roots that previously trusted a fast
## wrong-tail scout root. Reference values are certified by tbf01() below.
wrong_tail_tcrit_cases <- data.frame(
    n = c(28, 41),
    expected_tcrit = c(26.089665, 7.588086)
)
for (i in seq_len(nrow(wrong_tail_tcrit_cases))) {
    n <- wrong_tail_tcrit_cases$n[i]
    crit <- suppressWarnings(
        bfpwr:::tcrit(k = 30, n1 = n, n2 = n,
                      plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                      type = "two.sample", alternative = "less",
                      trange = "adaptive")
    )
    residual <- suppressWarnings(
        tbf01(t = crit, n1 = n, n2 = n,
              plocation = 0, pscale = 1/sqrt(2), pdf = 1,
              type = "two.sample", alternative = "less",
              log = TRUE) - log(30)
    )

    expect_equal(
        as.numeric(crit), wrong_tail_tcrit_cases$expected_tcrit[i],
        tolerance = 0.05,
        info = paste("less one-sided H0 tcrit should stay near the exact wrong-tail root at n =", n)
    )
    expect_true(
        is.finite(residual) && abs(residual) < 1e-4,
        info = paste("less one-sided H0 tcrit should satisfy BF01 = 30 at n =", n)
    )
}

## For nonzero nulls, one-sided H0 evidence must be computed after recentering
## the analysis prior around the tested null.
shifted_args <- list(n = 100, n1 = 100, n2 = 110, null = 0.2,
                     plocation = 0.5, pscale = 0.35, pdf = 30,
                     type = "two.sample", dpm = 0, dpsd = 0,
                     alternative = "greater")
shifted_h0 <- suppressWarnings(
    do.call(ptbf01, c(list(k = 30, lower.tail = FALSE), shifted_args))
)
neff <- 1 / (1 / shifted_args$n1 + 1 / shifted_args$n2)
se <- 1 / sqrt(neff)
root_fun <- function(est) {
    tbf01(t = (est - shifted_args$null) / se,
          n1 = shifted_args$n1,
          n2 = shifted_args$n2,
          plocation = shifted_args$plocation - shifted_args$null,
          pscale = shifted_args$pscale,
          pdf = shifted_args$pdf,
          type = shifted_args$type,
          alternative = shifted_args$alternative,
          log = TRUE) - log(30)
}
upper_root <- suppressWarnings(stats::uniroot(root_fun, c(-0.5, -0.2))$root)
shifted_h0_expected <- stats::pnorm(upper_root, mean = shifted_args$dpm,
                                    sd = se)
expect_true(is.finite(shifted_h0) && shifted_h0 > 0.001,
            info = "greater one-sided shifted-null H0 probability should not collapse to zero")
expect_true(abs(shifted_h0 - shifted_h0_expected) < 1e-5,
            info = "greater one-sided shifted-null H0 probability should use recentered critical value")
shifted_h1 <- suppressWarnings(
    do.call(ptbf01, c(list(k = 30, lower.tail = TRUE), shifted_args))
)
expect_equal(shifted_h0 + shifted_h1, 1, tolerance = 1e-12,
             info = "greater one-sided shifted-null H0/H1 probabilities should complement")

twosided_adaptive <- suppressWarnings(
    do.call(ptbf01, c(list(k = 1/10, alternative = "two.sided"), common_args))
)
twosided_wide_range <- suppressWarnings(
    do.call(ptbf01, c(list(k = 1/10, alternative = "two.sided", drange = c(-2, 2)),
                      common_args))
)
expect_true(is.finite(twosided_adaptive) && twosided_adaptive > 0 &&
                twosided_adaptive < 0.01,
            info = "two-sided adaptive search should use finite lower and upper roots")
expect_true(abs(twosided_adaptive - twosided_wide_range) < 5e-6,
            info = "two-sided adaptive search should match a bracketing two-root search")

tiny_upper_tail <- suppressWarnings(
    ptbf01(k = 3, n = 1000, plocation = 0, pscale = 1/sqrt(2), pdf = 1,
           dpm = 0.5, dpsd = 0, type = "two.sample",
           alternative = "two.sided", lower.tail = FALSE, drange = c(-2, 2))
)
expect_equal(tiny_upper_tail, 1.384599e-20, tolerance = 1e-6,
             info = "ptbf01 should compute very small upper-tail probabilities directly")

## ## do not run these tests for the moment, because they are there to verify
## ## the power with simulation which takes a long time to run

## library(tinytest)
## library(bfpwr)

## ## verify that computed power equal to simulated power
## grid <- expand.grid(n = c(50, 500, 5000), null = 0, plocation = 0,
##                     pscale = 0.71, pdf = 1, dpm = c(0.3, 0), dpsd = c(0, 0.1),
##                     alternative = c("two.sided", "greater", "less"),
##                     type = c("one.sample", "two.sample"),
##                     stringsAsFactors = FALSE)

## k <- c(1/4, 3)
## set.seed(42)
## nsim <- 1000
## pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
## pow1 <- numeric(nrow(grid))
## pow2 <- numeric(nrow(grid))
## powSim1 <- numeric(nrow(grid))
## powSim2 <- numeric(nrow(grid))
## for (i in seq(1, nrow(grid))) {
##     setTxtProgressBar(pb, i)
##     test <- paste(i, paste(colnames(grid), grid[i,], sep = "=", collapse = ", "))
##     ## simulate
##     if (grid$type[i] == "one.sample") {
##         se <- sqrt(1/grid$n[i])
##     } else {
##         se <- sqrt(2/grid$n[i])
##     }
##     estimatesim <- rnorm(n = nsim, mean = grid$dpm[i], sd = sqrt(grid$dpsd[i]^2 + se^2))
##     tsim <- (estimatesim - grid$null[i])/se
##     tbfsim <- tbf01(t = tsim, n = grid$n[i], plocation = grid$plocation[i],
##                     pscale = grid$pscale[i], pdf = grid$pdf[i],
##                     alternative = grid$alternative[i], type = grid$type[i])
##     if (grid$alternative[i] == "two.sided") {
##         powSim1[i] <- mean(tbfsim <= k[1])
##         powSim2[i] <- mean(tbfsim <= k[2])
##     } else if (grid$alternative[i] == "less") {
##         powSim1[i] <- mean(tbfsim <= k[1] & tsim <= 0)
##         powSim2[i] <- mean(tbfsim <= k[2] & tsim <= 0)
##     } else {
##         powSim1[i] <- mean(tbfsim <= k[1] & tsim >= 0)
##         powSim2[i] <- mean(tbfsim <= k[2] & tsim >= 0)
##     }

##     ## compute numerically
##     pows <- ptbf01(k = k, n = grid$n[i], null = grid$null[i],
##                    plocation = grid$plocation[i], pscale = grid$pscale[i],
##                    pdf = grid$pdf[i], dpm = grid$dpm[i], dpsd = grid$dpsd[i],
##                    alternative = grid$alternative[i], type = grid$type[i],
##                    lower.tail = TRUE, drange = "adaptive")
##     pow1[i] <- pows[1]
##     pow2[i] <- pows[2]
##     cat("\n")
##     print(expect_equal(pow1[i], powSim1[i], info = test, tolerance = 0.04))
##     print(expect_equal(pow2[i], powSim2[i], info = test, tolerance = 0.04))
## }
