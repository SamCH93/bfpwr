# Directed checks of arXiv:2601.02851v2, Table 1 and Sections 2/3.1/3.2.
# Run with an installed bfpwr on the library path; first argument is output dir.
library(bfpwr)
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else file.path(tempdir(), "moment-paper-check")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(digits = 17, warn = 1)
cat(R.version.string, "\nPackage:", find.package("bfpwr"), "\nLibraries:\n")
print(.libPaths())

# Literal Table 1 expression; no package boundary or BF helpers are used.
table_logbf <- function(z, r) {
    q <- z/sqrt(1 + 1/r^2)
    denominator <- 2*exp(q^2/2)*(1 + q^2)*pnorm(q) + q*sqrt(2/pi)
    if (any(!is.finite(denominator) | denominator <= 0)) return(NA_real_)
    1.5*log1p(r^2) - log(denominator)
}
quadrature_logbf <- function(z, r) {
    q <- z/sqrt(1 + 1/r^2)
    marginal <- integrate(function(t) t^2*exp(q*t - t^2/2), 0, Inf,
                          rel.tol = 1e-11, subdivisions = 1000L)$value
    1.5*log1p(r^2) + .5*log(2*pi) - log(2) - log(marginal)
}
write_result <- function(x, name) {
    write.csv(x, file.path(output, paste0(name, ".csv")), row.names = FALSE)
    invisible(x)
}
checks <- list()
record <- function(name, error, tolerance, count = length(error)) {
    checks[[name]] <<- data.frame(check = name, count = count,
        max_abs_error = max(abs(error)), tolerance = tolerance,
        passed = all(is.finite(error) & abs(error) <= tolerance))
}

# 1. Formula, independent marginal integral, normalization, and identities.
bf <- expand.grid(r = c(.01, .1, .5, 1, 3, 10),
    z = c(-5.001, -5, -4.999, -2, -1, 0, 1, 2, 5, 8),
    alternative = c("greater", "less"), null = c(0, 2),
    stringsAsFactors = FALSE)
bf$table <- mapply(table_logbf, bf$z, bf$r)
bf$quadrature <- mapply(quadrature_logbf, bf$z, bf$r)
bf$actual <- mapply(function(z, r, alternative, null) {
    direction <- if (alternative == "greater") 1 else -1
    nmbf01(estimate = null + direction*z, se = 1, null = null, psd = r,
           log = TRUE, alternative = alternative)
}, bf$z, bf$r, bf$alternative, bf$null)
bf$table_error <- bf$table - bf$quadrature
bf$actual_error <- bf$actual - bf$quadrature
write_result(bf, "formula")
record("Table 1 versus quadrature (log BF)", bf$table_error, 1e-8)
record("nmbf01 versus quadrature (log BF)", bf$actual_error, 1e-8)
normalization <- data.frame(r = c(.01, .1, .5, 1, 3, 10))
normalization$integral <- vapply(normalization$r, function(r) {
    integrate(function(delta) 2*delta^2*dnorm(delta, sd = r)/r^2,
              0, Inf, rel.tol = 1e-11)$value
}, numeric(1))
write_result(normalization, "normalization")
record("Normalized positive half-moment prior", normalization$integral - 1, 1e-10)
zero <- bf[bf$z == 0, ]
record("z=0 closed form (log BF)", zero$actual - 1.5*log1p(zero$r^2), 1e-10)
harmonic <- unique(bf[c("r", "z", "null")])
harmonic$error <- mapply(function(r, z, null) {
    logs <- nmbf01(estimate = null + z, se = 1, null = null, psd = r,
                  log = TRUE, alternative = c("greater", "less"))
    m <- max(-logs)
    log_harmonic <- log(2) - m - log(sum(exp(-logs - m)))
    log_harmonic - nmbf01(null + z, se = 1, null = null, psd = r, log = TRUE)
}, harmonic$r, harmonic$z, harmonic$null)
write_result(harmonic, "harmonic")
record("Harmonic one-sided/two-sided identity (log BF)", harmonic$error, 1e-8)
cat("Formula checks complete.\n")

# Validate endpoints independently; use the integral only if literal tail
# cancellation is observed. Roots never use the package's zcrit helper.
roots <- expand.grid(n = c(20, 50, 60, 150), k = c(.1, 3))
roots$r <- .5/sqrt(2/roots$n)
roots$left_table <- mapply(table_logbf, -20, roots$r)
roots$right_table <- mapply(table_logbf, 20, roots$r)
roots$left_quad <- mapply(quadrature_logbf, -20, roots$r)
roots$right_quad <- mapply(quadrature_logbf, 20, roots$r)
roots$use_quadrature <- !is.finite(roots$left_table) |
    !is.finite(roots$right_table) |
    abs(roots$left_table - roots$left_quad) > 1e-8 |
    abs(roots$right_table - roots$right_quad) > 1e-8
stopifnot(all(roots$left_quad > log(roots$k)),
          all(roots$right_quad < log(roots$k)))
roots$z <- mapply(function(r, k, use_quadrature) {
    bf_fun <- if (use_quadrature) quadrature_logbf else table_logbf
    uniroot(function(z) bf_fun(z, r) - log(k), c(-20, 20), tol = 1e-12)$root
}, roots$r, roots$k, roots$use_quadrature)
roots$residual <- mapply(quadrature_logbf, roots$z, roots$r) - log(roots$k)
write_result(roots, "roots")
record("Independent critical-value quadrature residual", roots$residual, 1e-8)
critical <- function(n, k) roots$z[roots$n == n & roots$k == k]
fixed_probability <- function(n, k, null, dpsd, alternative) {
    se <- sqrt(2/n)
    direction <- if (alternative == "greater") 1 else -1
    Zcrit <- (null + direction*se*critical(n, k) - (null + .3))/
        sqrt(se^2 + dpsd^2)
    pnorm(Zcrit, lower.tail = alternative == "less")
}

# 2. Fixed sample probabilities, both BF tails, and one-look agreement.
fixed <- expand.grid(k = c(.1, 3), n = c(20, 60, 150), null = c(0, 2),
    dpsd = c(0, .2), alternative = c("greater", "less"),
    lower.tail = c(TRUE, FALSE), stringsAsFactors = FALSE)
fixed$expected <- mapply(function(n, k, null, dpsd, alternative, lower.tail) {
    p <- fixed_probability(n, k, null, dpsd, alternative)
    if (lower.tail) p else 1 - p
}, fixed$n, fixed$k, fixed$null, fixed$dpsd, fixed$alternative, fixed$lower.tail)
fixed$actual <- mapply(function(n, k, null, dpsd, alternative, lower.tail) {
    pnmbf01(k, n, usd = sqrt(2), null = null, psd = .5, dpm = null + .3,
            dpsd = dpsd, alternative = alternative, lower.tail = lower.tail)
}, fixed$n, fixed$k, fixed$null, fixed$dpsd, fixed$alternative, fixed$lower.tail)
fixed$error <- fixed$actual - fixed$expected
write_result(fixed, "fixed")
record("pnmbf01 independent normal-tail probabilities", fixed$error, 1e-8)
one_look <- unique(fixed[c("n", "null", "dpsd", "alternative")])
one_look$H1_expected <- one_look$H0_expected <- NA_real_
one_look$H1_actual <- one_look$H0_actual <- NA_real_
for (j in seq_len(nrow(one_look))) {
    x <- one_look[j, ]
    actual <- pbf01seq(k1 = .1, k0 = 3, se = sqrt(2/x$n), n = x$n,
        null = x$null, psd = .5, dpm = x$null + .3, dpsd = x$dpsd,
        type = "moment", alternative = x$alternative)
    one_look$H1_expected[j] <- fixed_probability(x$n, .1, x$null, x$dpsd, x$alternative)
    one_look$H0_expected[j] <- 1 - fixed_probability(x$n, 3, x$null, x$dpsd, x$alternative)
    one_look$H1_actual[j] <- actual$cumpH1
    one_look$H0_actual[j] <- actual$cumpH0
}
write_result(one_look, "one-look")
record("One-look pbf01seq H1/H0 probabilities", c(one_look$H1_actual - one_look$H1_expected,
    one_look$H0_actual - one_look$H0_expected), 1e-8)
cat("Fixed and one-look checks complete.\n")

# 3. Eq. 3.1 covariance and conditional-normal two-look integration.
sequential <- expand.grid(null = c(0, 2), dpsd = c(0, .2),
    alternative = c("greater", "less"), stringsAsFactors = FALSE)
sequential_rows <- moment_rows <- list()
for (j in seq_len(nrow(sequential))) {
    x <- sequential[j, ]
    n <- c(20, 50)
    se <- sqrt(2/n)
    direction <- if (x$alternative == "greater") 1 else -1
    mu <- direction*.3/se
    V <- outer(n, n, function(a, b) sqrt(pmin(a, b)/pmax(a, b))) +
        x$dpsd^2/outer(se, se)
    c0 <- vapply(n, critical, numeric(1), k = 3)
    c1 <- vapply(n, critical, numeric(1), k = .1)
    conditional_sd <- sqrt(V[2, 2] - V[1, 2]^2/V[1, 1])
    second <- function(boundary, lower.tail) {
        integrate(function(z1) {
            conditional_mean <- mu[2] + V[1, 2]/V[1, 1]*(z1 - mu[1])
            dnorm(z1, mu[1], sqrt(V[1, 1]))*
                pnorm(boundary, conditional_mean, conditional_sd,
                      lower.tail = lower.tail)
        }, c0[1], c1[1], rel.tol = 1e-10)$value
    }
    h1 <- c(pnorm(c1[1], mu[1], sqrt(V[1, 1]), lower.tail = FALSE),
            second(c1[2], FALSE))
    h0 <- c(pnorm(c0[1], mu[1], sqrt(V[1, 1])), second(c0[2], TRUE))
    actual <- pbf01seq(k1 = .1, k0 = 3, se = se, n = n, null = x$null,
        psd = .5, dpm = x$null + .3, dpsd = x$dpsd,
        type = "moment", alternative = x$alternative, ngrid = 50000)
    expected <- c(h1, h0, cumsum(h1), cumsum(h0), 1 - cumsum(h1) - cumsum(h0))
    observed <- c(diff(c(0, actual$cumpH1)), diff(c(0, actual$cumpH0)),
        actual$cumpH1, actual$cumpH0, actual$cumpInc)
    sequential_rows[[j]] <- data.frame(null = x$null, dpsd = x$dpsd,
        alternative = x$alternative,
        quantity = rep(c("stageH1", "stageH0", "cumpH1", "cumpH0", "cumpInc"), each = 2),
        look = rep(1:2, 5), expected = expected, actual = observed,
        error = observed - expected)
    stop1 <- h1[1] + h0[1]
    expected_moments <- c(n[1]*stop1 + n[2]*(1 - stop1),
        (n[2] - n[1])^2*stop1*(1 - stop1))
    moment_rows[[j]] <- data.frame(null = x$null, dpsd = x$dpsd,
        alternative = x$alternative, quantity = c("EN", "VarN"),
        expected = expected_moments, actual = c(actual$EN, actual$VarN),
        error = c(actual$EN, actual$VarN) - expected_moments,
        tolerance = 1e-5*c(max(n), max(n)^2))
    cat("Two-look case", j, "of", nrow(sequential), "complete.\n")
}
sequential_rows <- write_result(do.call(rbind, sequential_rows), "two-look")
moment_rows <- write_result(do.call(rbind, moment_rows), "moments")
record("Two-look stage/cumulative/inconclusive probabilities", sequential_rows$error, 1e-5)
record("Two-look expected sample size", subset(moment_rows, quantity == "EN")$error, 5e-4)
record("Two-look sample-size variance", subset(moment_rows, quantity == "VarN")$error, .025)
translation <- merge(subset(sequential_rows, null == 0), subset(sequential_rows, null == 2),
    by = c("dpsd", "alternative", "quantity", "look"), suffixes = c("_zero", "_two"))
translation$error <- translation$actual_two - translation$actual_zero
write_result(translation, "translation")
record("Translated-null sequential identity", translation$error, 1e-10)

# 4. Check the exhaustive first crossing and record the fast solver result.
search <- data.frame(method = c("nnmbf01", "nbf01seq exhaustive one look"),
    expected = c(NA_real_, 18), actual = c(
    nnmbf01(.1, .016, 1, psd = .3, dpm = .03, dpsd = .03, alternative = "greater"),
    nbf01seq(k1 = .1, k0 = 10, power = .016, usd = 1, psd = .3, dpm = .03,
        dpsd = .03, type = "moment", alternative = "greater", looks = 1,
        search = "exhaustive", nrange = c(2, 1000))))
write_result(search, "held-search")
record("Exhaustive one-look search returns the first crossing",
    search$actual[2L] - search$expected[2L], 0)
summary <- write_result(do.call(rbind, checks), "summary")
cat("\nSummary:\n")
print(summary, row.names = FALSE)
cat("\nRoot endpoint diagnostics:\n")
print(roots)
cat("\nFixed and exhaustive search results:\n")
print(search, row.names = FALSE)
saveRDS(list(formula = bf, normalization = normalization, harmonic = harmonic,
    roots = roots, fixed = fixed, one_look = one_look, sequential = sequential_rows,
    moments = moment_rows, translation = translation, search = search, summary = summary),
    file.path(output, "results.rds"))
if (!all(summary$passed)) stop("Directed manuscript check failed; inspect saved results.")
cat("\nAll directed checks passed.\n")
