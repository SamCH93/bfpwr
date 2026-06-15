library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)

## Tests tbf01 API behavior and stable one-sided/two-sided tail calculations.
## Manuscript source: informed/JZS t BF section in paper/bfssd.Rnw 1481-1530 and
## the one-sided example at 1609-1627; extreme-tail numbers are package regressions.

res <- tbf01(t = c(-1, 0, 1), n = 100, plocation = 0, pscale = 1, pdf = 1,
             type = "one.sample", alternative = "two.sided", log = FALSE)
logres <- tbf01(t = c(-1, 0, 1), n = 100, plocation = 0, pscale = 1, pdf = 1,
                type = "one.sample", alternative = "two.sided", log = TRUE)

expect_true(is.numeric(res), info = "tbf01 should return a numeric value")

expect_true(length(res) == 3, info = "tbf01 should handle vector inputs")

expect_equal(log(res), logres,
             info = "tbf01 should return log(tbf01) when log = TRUE")

tcrit_search_limit <- bfpwr:::.bfpwr_one_sided_tail_limits(
    origin = 0, step_scale = 1, mean = 0, sd = 20, tail.eps = 1e-6
)

less_crit <- suppressWarnings(
    bfpwr:::tcrit(k = 30, n1 = 41, n2 = 41, plocation = 0,
                  pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
                  alternative = "less", trange = "adaptive",
                  search_limit = tcrit_search_limit)
)
greater_crit <- suppressWarnings(
    bfpwr:::tcrit(k = 30, n1 = 41, n2 = 41, plocation = 0,
                  pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
                  alternative = "greater", trange = "adaptive",
                  search_limit = tcrit_search_limit)
)
less_residual <- suppressWarnings(
    tbf01(t = less_crit, n1 = 41, n2 = 41, plocation = 0,
          pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
          alternative = "less", log = TRUE) - log(30)
)
greater_residual <- suppressWarnings(
    tbf01(t = greater_crit, n1 = 41, n2 = 41, plocation = 0,
          pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
          alternative = "greater", log = TRUE) - log(30)
)
expect_true(less_crit > 4 && greater_crit < -4,
            info = "one-sided adaptive tcrit should search the wrong tail for H0 evidence")
expect_equal(less_crit, -greater_crit, tolerance = 1e-5,
             info = "mirrored one-sided adaptive tcrit roots should agree")
expect_true(max(abs(c(less_residual, greater_residual))) < 1e-4,
            info = "one-sided adaptive tcrit roots should satisfy BF01 threshold")

tcrit_missing_limit_warning <- NULL
tcrit_missing_limit <- withCallingHandlers(
    bfpwr:::tcrit(k = 30, n1 = 41, n2 = 41, plocation = 0,
                  pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
                  alternative = "greater", trange = "adaptive"),
    warning = function(w) {
        tcrit_missing_limit_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    }
)
expect_true(is.nan(tcrit_missing_limit),
            info = "one-sided adaptive tcrit should require an explicit search limit")
expect_true(grepl("requires 'search_limit'", tcrit_missing_limit_warning,
                  fixed = TRUE),
            info = "one-sided adaptive tcrit should warn clearly when search_limit is missing")

seen_x <- numeric(0)
always_positive <- function(x) {
    seen_x <<- c(seen_x, x)
    1
}
no_root <- bfpwr:::.bfpwr_one_sided_adaptive_root(
    certify_fun = always_positive,
    scout_fun = always_positive,
    alternative = "greater",
    origin = 0,
    step_scale = 1,
    search_limit = 1
)
expect_true(inherits(no_root$root, "try-error"),
            info = "one-sided helper should report no root when expected side has no crossing")
expect_true(all(seen_x >= 0),
            info = "one-sided helper should not retry the opposite side by default")

nonfinite_limit_fun <- function(x) {
    if (x >= 10) {
        return(NaN)
    }
    1
}
nonfinite_limit <- bfpwr:::.bfpwr_one_sided_adaptive_root(
    certify_fun = nonfinite_limit_fun,
    scout_fun = nonfinite_limit_fun,
    alternative = "greater",
    origin = 0,
    step_scale = 1,
    search_limit = 10,
    steps = 1,
    scout_tail_steps = 1,
    tail_steps = 10
)
expect_true(inherits(nonfinite_limit$root, "try-error"),
            info = "one-sided helper should report no root with non-finite limit values")
expect_false(nonfinite_limit$search_limit_reached,
             info = "one-sided helper should not convert non-finite limit values to tail cutoffs")

search_root <- function(x) 1 - x
final_root <- function(x) 2 - x
full_certified_root <- bfpwr:::.bfpwr_one_sided_adaptive_root(
    certify_fun = final_root,
    search_fun = search_root,
    scout_fun = search_root,
    alternative = "greater",
    origin = 0,
    step_scale = 1,
    search_limit = 3
)
expect_equal(
    full_certified_root$root,
    2,
    tolerance = 1e-8,
    info = "one-sided helper should return the root certified by the final function"
)
expect_false(
    full_certified_root$search_limit_reached,
    info = "one-sided helper should not treat a relaxed-search root as final certification"
)

search_finite_limit <- function(x) 1
final_nonfinite_limit <- function(x) {
    if (x >= 10) {
        return(NaN)
    }
    1
}
nonfinite_final_limit <- bfpwr:::.bfpwr_one_sided_adaptive_root(
    certify_fun = final_nonfinite_limit,
    search_fun = search_finite_limit,
    scout_fun = search_finite_limit,
    alternative = "greater",
    origin = 0,
    step_scale = 1,
    search_limit = 10,
    steps = 1,
    scout_tail_steps = 1,
    tail_steps = 10
)
expect_true(
    inherits(nonfinite_final_limit$root, "try-error"),
    info = "one-sided helper should fail when the final function is non-finite at the limit"
)
expect_false(
    nonfinite_final_limit$search_limit_reached,
    info = "one-sided helper should not report a tail cutoff from only the relaxed search function"
)

flat_missing_root <- bfpwr:::.bfpwr_one_sided_adaptive_root(
    certify_fun = function(x) -1 - 1/(abs(x) + 1),
    scout_fun = function(x) -1 - 1/(abs(x) + 1),
    alternative = "greater",
    origin = 0,
    step_scale = 1,
    search_limit = 64,
    steps = 1,
    scout_tail_steps = 1,
    tail_steps = 64
)
expect_equal(
    flat_missing_root$status,
    "impossible",
    info = "one-sided helper should distinguish flat unattainable roots from finite cutoffs"
)

near_limit_missing_root <- bfpwr:::.bfpwr_one_sided_adaptive_root(
    certify_fun = function(x) -1/(abs(x) + 1),
    scout_fun = function(x) -1/(abs(x) + 1),
    alternative = "greater",
    origin = 0,
    step_scale = 1,
    search_limit = 64,
    steps = 1,
    scout_tail_steps = 1,
    tail_steps = 64
)
expect_equal(
    near_limit_missing_root$status,
    "tail_cutoff",
    info = "one-sided helper should not call near-threshold finite cutoffs impossible"
)

tcrit_impossible_warning <- NULL
tcrit_impossible <- withCallingHandlers(
    bfpwr:::tcrit(k = 20, n1 = 15, n2 = 15, plocation = 0,
                  pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
                  alternative = "greater", trange = "adaptive",
                  search_limit = 64),
    warning = function(w) {
        tcrit_impossible_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    }
)
expect_true(is.nan(tcrit_impossible),
            info = "one-sided tcrit should return NaN when BF01 = k is unattainable")
expect_true(grepl("appears unattainable", tcrit_impossible_warning,
                  fixed = TRUE),
            info = "one-sided tcrit should not recommend widening trange for unattainable roots")

tcrit_impossible_result <- suppressWarnings(
    bfpwr:::.bfpwr_tcrit_result(
        k = 20, n1 = 15, n2 = 15, plocation = 0,
        pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
        alternative = "greater", trange = "adaptive",
        search_limit = 64
    )
)
expect_equal(
    tcrit_impossible_result$status,
    "impossible",
    info = "tcrit result should classify one-sided unattainable roots structurally"
)
expect_equal(
    tcrit_impossible_result$reason,
    "one_sided_unattainable",
    info = "tcrit result should expose the status reason code"
)

seq_impossible_warning <- NULL
withCallingHandlers(
    bfpwr:::.bfseq_warn_t_boundary_statuses(
        results0 = list(list(value = NaN, status = "impossible",
                             warnings = "BF01 = k appears unattainable")),
        results1 = list(list(value = 1, status = "ok",
                             warnings = character())),
        tail.eps = 1e-3
    ),
    warning = function(w) {
        seq_impossible_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    }
)
expect_true(grepl("No H0 sequential t stopping boundary exists",
                  seq_impossible_warning, fixed = TRUE),
            info = "sequential t diagnostics should aggregate impossible H0 boundaries separately")
expect_false(grepl("Pass a wider", seq_impossible_warning, fixed = TRUE),
             info = "impossible H0 boundary diagnostics should not recommend wider trange")

expect_equal(
    bfpwr:::.bfpwr_tcrit_status_from_result(
        value = 1,
        issues = list(bfpwr:::.bfpwr_tcrit_issue(
            code = "critical_value_failed",
            status = "search_failed",
            message = "Numerical problems finding critical value"
        ))
    ),
    "search_failed",
    info = "tcrit status should not ignore structured failure issues for finite values"
)
expect_equal(
    bfpwr:::.bfpwr_tcrit_status_from_result(value = numeric(0),
                                            issues = list()),
    "search_failed",
    info = "tcrit status should not treat empty results as valid"
)

if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "remaining tbf01 numerical-stability checks are extended"
    ))
}

expect_equal(
    tbf01(t = -20, n1 = 7880, n2 = 7880, alternative = "greater",
          type = "two.sample", log = TRUE),
    7.2302204, tolerance = 1e-5,
    info = "tbf01 should use a stable wrong-tail one-sided calculation"
)

expect_equal(
    tbf01(t = -4, n1 = 53, n2 = 57, plocation = 0.350, pscale = 0.102,
          pdf = 3, alternative = "greater", type = "two.sample", log = TRUE),
    4.3357318, tolerance = 1e-5,
    info = "tbf01 should handle shifted informed priors in the wrong tail"
)

expect_equal(
    tbf01(t = 1, n1 = 1e6, n2 = 1e6, pscale = 1, type = "two.sample",
          log = TRUE),
    6.2869769, tolerance = 1e-4,
    info = "tbf01 should return finite large-n two-sided log Bayes factors"
)

opposite_less <- tbf01(t = 7.792904, n1 = 500, n2 = 500,
                       plocation = 0, pscale = 1 / sqrt(2), pdf = 1,
                       type = "two.sample", alternative = "less",
                       log = TRUE)
opposite_greater <- tbf01(t = -7.792904, n1 = 500, n2 = 500,
                          plocation = 0, pscale = 1 / sqrt(2), pdf = 1,
                          type = "two.sample", alternative = "greater",
                          log = TRUE)
expect_true(is.finite(opposite_less) && opposite_less > 0,
            info = "one-sided opposite-direction tbf01 should be finite on the log scale")
expect_equal(opposite_less, opposite_greater, tolerance = 1e-8,
             info = "mirrored one-sided opposite-direction tbf01 values should agree")
