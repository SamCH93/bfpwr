library(tinytest)
library(bfpwr)

## Strict fixed-sample manuscript checks. The helper functions below rederive the
## paper formulas outside bfpwr, then the test blocks compare package results to
## those references. Source line references are to paper/bfssd.Rnw unless noted.

## Independent numeric helpers for log-space arithmetic, normal/beta interval
## probabilities, direct BF formulas, and root searches. These are test oracles,
## not manuscript examples.
logspace_sub <- function(logx, logy) {
    if (logy > logx) return(NaN)
    if (logy == -Inf) return(logx)
    if (logx == logy) return(-Inf)
    logx + log1p(-exp(logy - logx))
}

lpnorm_interval <- function(lower, upper, mean = 0, sd = 1) {
    if (is.infinite(lower) && lower < 0 && is.infinite(upper) && upper > 0) {
        return(0)
    }
    if (is.infinite(lower) && lower < 0) {
        return(stats::pnorm(upper, mean = mean, sd = sd, log.p = TRUE))
    }
    if (is.infinite(upper) && upper > 0) {
        return(stats::pnorm(lower, mean = mean, sd = sd, lower.tail = FALSE,
                            log.p = TRUE))
    }
    left <- logspace_sub(
        stats::pnorm(upper, mean = mean, sd = sd, log.p = TRUE),
        stats::pnorm(lower, mean = mean, sd = sd, log.p = TRUE)
    )
    right <- logspace_sub(
        stats::pnorm(lower, mean = mean, sd = sd, lower.tail = FALSE,
                     log.p = TRUE),
        stats::pnorm(upper, mean = mean, sd = sd, lower.tail = FALSE,
                     log.p = TRUE)
    )
    max(left, right, na.rm = TRUE)
}

lpbeta_interval <- function(lower, upper, shape1, shape2) {
    if (lower <= 0 && upper >= 1) return(0)
    if (lower <= 0) {
        return(stats::pbeta(upper, shape1 = shape1, shape2 = shape2,
                            log.p = TRUE))
    }
    if (upper >= 1) {
        return(stats::pbeta(lower, shape1 = shape1, shape2 = shape2,
                            lower.tail = FALSE, log.p = TRUE))
    }
    left <- logspace_sub(
        stats::pbeta(upper, shape1 = shape1, shape2 = shape2, log.p = TRUE),
        stats::pbeta(lower, shape1 = shape1, shape2 = shape2, log.p = TRUE)
    )
    right <- logspace_sub(
        stats::pbeta(lower, shape1 = shape1, shape2 = shape2,
                     lower.tail = FALSE, log.p = TRUE),
        stats::pbeta(upper, shape1 = shape1, shape2 = shape2,
                     lower.tail = FALSE, log.p = TRUE)
    )
    max(left, right, na.rm = TRUE)
}

ref_log_bf01 <- function(estimate, se, null, pm, psd) {
    stats::dnorm(estimate, mean = null, sd = se, log = TRUE) -
        stats::dnorm(estimate, mean = pm, sd = sqrt(se^2 + psd^2),
                     log = TRUE)
}

ref_log_nmbf01 <- function(estimate, se, null, psd) {
    marginal_sd <- sqrt(se^2 + psd^2)
    post_var <- 1/(1/se^2 + 1/psd^2)
    post_mean_diff <- psd^2/(se^2 + psd^2) * (estimate - null)
    moment_factor <- (post_var + post_mean_diff^2)/psd^2
    stats::dnorm(estimate, mean = null, sd = se, log = TRUE) -
        (stats::dnorm(estimate, mean = null, sd = marginal_sd, log = TRUE) +
         log(moment_factor))
}

normal_probability_from_logbf <- function(logbf, k, mean, sd,
                                          lower.tail = TRUE) {
    f <- function(x) {
        vapply(x, function(xi) logbf(xi) - log(k), numeric(1))
    }
    grid <- seq(mean - 12*sd, mean + 12*sd, length.out = 2001)
    values <- f(grid)
    roots <- numeric()
    for (i in seq_len(length(grid) - 1)) {
        f0 <- values[i]
        f1 <- values[i + 1]
        if (!is.finite(f0) || !is.finite(f1)) next
        if (f0 == 0) {
            roots <- c(roots, grid[i])
        } else if (f0*f1 < 0) {
            roots <- c(roots, stats::uniroot(function(x) f(x), grid[c(i, i + 1)],
                                             tol = 1e-12)$root)
        }
    }
    roots <- sort(unique(round(roots, 12)))
    cuts <- c(-Inf, roots, Inf)
    prob <- 0
    for (i in seq_len(length(cuts) - 1)) {
        lower <- cuts[i]
        upper <- cuts[i + 1]
        midpoint <- if (is.infinite(lower) && lower < 0) {
            upper - max(1, sd)
        } else if (is.infinite(upper) && upper > 0) {
            lower + max(1, sd)
        } else {
            (lower + upper)/2
        }
        if (f(midpoint) <= 0) {
            prob <- prob + exp(lpnorm_interval(lower, upper, mean = mean,
                                               sd = sd))
        }
    }
    prob <- min(max(prob, 0), 1)
    if (lower.tail) prob else 1 - prob
}

ref_pbf01 <- function(k, n, usd, null, pm, psd, dpm, dpsd,
                      lower.tail = TRUE) {
    se <- usd/sqrt(n)
    design_sd <- sqrt(se^2 + dpsd^2)
    normal_probability_from_logbf(
        logbf = function(estimate) ref_log_bf01(estimate, se, null, pm, psd),
        k = k, mean = dpm, sd = design_sd, lower.tail = lower.tail
    )
}

ref_pnmbf01 <- function(k, n, usd, null, psd, dpm, dpsd,
                        lower.tail = TRUE) {
    se <- usd/sqrt(n)
    design_sd <- sqrt(se^2 + dpsd^2)
    normal_probability_from_logbf(
        logbf = function(estimate) ref_log_nmbf01(estimate, se, null, psd),
        k = k, mean = dpm, sd = design_sd, lower.tail = lower.tail
    )
}

ref_log_tbf01 <- function(t, n1, n2, plocation, pscale, pdf,
                          type = "two.sample",
                          alternative = "two.sided") {
    if (type == "two.sample") {
        df <- n1 + n2 - 2
        neff <- 1/(1/n1 + 1/n2)
    } else {
        df <- n1 - 1
        neff <- n1
    }
    lower <- -Inf
    upper <- Inf
    norm_const <- 1
    if (alternative == "greater") {
        lower <- 0
        norm_const <- stats::pt((0 - plocation)/pscale, df = pdf,
                                lower.tail = FALSE)
    }
    if (alternative == "less") {
        upper <- 0
        norm_const <- stats::pt((0 - plocation)/pscale, df = pdf)
    }
    prior <- function(delta) {
        stats::dt((delta - plocation)/pscale, df = pdf)/
            (pscale*norm_const)
    }
    integrand <- function(delta) {
        suppressWarnings(stats::dt(t, df = df, ncp = sqrt(neff)*delta)) *
            prior(delta)
    }
    f1 <- stats::integrate(integrand, lower = lower, upper = upper,
                           rel.tol = 1e-10, subdivisions = 2000)$value
    stats::dt(t, df = df, log = TRUE) - log(f1)
}

ref_ptbf01_greater <- function(k, n, null, plocation, pscale, pdf,
                               dpm, dpsd, type = "two.sample",
                               lower.tail = TRUE) {
    n1 <- n2 <- n
    neff <- if (type == "two.sample") 1/(1/n1 + 1/n2) else n1
    se <- 1/sqrt(neff)
    root <- stats::uniroot(
        function(t) {
            ref_log_tbf01(t, n1 = n1, n2 = n2, plocation = plocation,
                          pscale = pscale, pdf = pdf, type = type,
                          alternative = "greater") - log(k)
        },
        interval = c(0, 10),
        tol = 1e-8
    )$root
    crit_est <- null + root*se
    logpow <- stats::pnorm(crit_est, mean = dpm,
                           sd = sqrt(se^2 + dpsd^2),
                           lower.tail = FALSE, log.p = TRUE)
    if (lower.tail) exp(logpow) else 1 - exp(logpow)
}

ref_log_binbf01 <- function(x, n, p0, type, a, b) {
    if (type == "point") {
        x*log(p0) + (n - x)*log1p(-p0) + lbeta(a, b) -
            lbeta(a + x, b + n - x)
    } else {
        lpost0 <- stats::pbeta(p0, shape1 = a + x, shape2 = b + n - x,
                               log.p = TRUE)
        lpost1 <- stats::pbeta(p0, shape1 = a + x, shape2 = b + n - x,
                               lower.tail = FALSE, log.p = TRUE)
        lprior0 <- stats::pbeta(p0, shape1 = a, shape2 = b, log.p = TRUE)
        lprior1 <- stats::pbeta(p0, shape1 = a, shape2 = b,
                                lower.tail = FALSE, log.p = TRUE)
        lpost0 - lpost1 + lprior1 - lprior0
    }
}

ref_pbinbf01 <- function(k, n, p0, type, a, b, dp = NA, da = a, db = b,
                         dl = 0, du = 1, lower.tail = TRUE) {
    n <- ceiling(n)
    x <- 0:n
    success <- ref_log_binbf01(x, n, p0, type, a, b) <= log(k)
    if (!is.na(dp)) {
        logpmf <- stats::dbinom(x, size = n, prob = dp, log = TRUE)
    } else {
        log_norm <- lpbeta_interval(dl, du, shape1 = da, shape2 = db)
        logpmf <- lchoose(n, x) + lbeta(da + x, db + n - x) -
            lbeta(da, db) +
            vapply(x, function(xi) {
                lpbeta_interval(dl, du, shape1 = da + xi,
                                shape2 = db + n - xi)
            }, numeric(1)) -
            log_norm
    }
    prob <- sum(exp(logpmf[success]))
    prob <- min(max(prob, 0), 1)
    if (lower.tail) prob else 1 - prob
}

ref_stable_n <- function(pfun, power, nrange = c(1, 500), nextend = 10) {
    vals <- vapply(seq(nrange[1], nrange[2]), function(n) pfun(n),
                   numeric(1))
    ns <- seq(nrange[1], nrange[2])
    candidates <- which(vals >= power)
    candidates <- candidates[candidates + nextend <= length(vals)]
    ns[candidates[vapply(candidates, function(i) {
        all(vals[i:(i + nextend)] >= power)
    }, logical(1))][1]]
}

## Tests BF01 value and log value from eq. BF01 and the package reference at
## bfssd.Rnw 354-369. Inputs are synthetic formula checks, not rounded text.
expect_equal(
    bf01(estimate = 0.23, se = 0.14, null = 0.1, pm = 0.4, psd = 0.25),
    exp(ref_log_bf01(estimate = 0.23, se = 0.14, null = 0.1, pm = 0.4,
                     psd = 0.25)),
    tolerance = 1e-12,
    info = "bf01 agrees with the normal density-ratio formula"
)
expect_equal(
    bf01(estimate = 0.23, se = 0.14, null = 0.1, pm = 0.4, psd = 0.25,
         log = TRUE),
    ref_log_bf01(estimate = 0.23, se = 0.14, null = 0.1, pm = 0.4,
                 psd = 0.25),
    tolerance = 1e-12,
    info = "bf01 log output agrees with the log density-ratio formula"
)

## Tests fixed-n pbf01 lower/upper tails from the power functions prLR/prBF
## (bfssd.Rnw 449-466 and 590-606), as used in the plot-power chunk
## (bfssd.Rnw 501-558).
pbf_expected <- ref_pbf01(k = 1/6, n = 80, usd = sqrt(2), null = 0,
                          pm = 0.3, psd = 0.2, dpm = 0.35, dpsd = 0.05)
expect_equal(
    pbf01(k = 1/6, n = 80, usd = sqrt(2), null = 0, pm = 0.3,
          psd = 0.2, dpm = 0.35, dpsd = 0.05),
    pbf_expected,
    tolerance = 1e-10,
    info = "pbf01 agrees with root-partitioned normal predictive integration"
)
expect_equal(
    pbf01(k = 1/6, n = 80, usd = sqrt(2), null = 0, pm = 0.3,
          psd = 0.2, dpm = 0.35, dpsd = 0.05, lower.tail = FALSE),
    1 - pbf_expected,
    tolerance = 1e-10,
    info = "pbf01 upper tail complements the reference lower tail"
)

## Tests nbf01 point-prior sample-size root against the sample-size section and
## verify-equations chunk (bfssd.Rnw 640-696 and 711-719).
nbf_expected <- stats::uniroot(
    function(n) {
        ref_pbf01(k = 1/10, n = n, usd = sqrt(2), null = 0, pm = 0.5,
                  psd = 0, dpm = 0.5, dpsd = 0) - 0.8
    },
    interval = c(1, 200)
)$root
expect_equal(
    nbf01(k = 1/10, power = 0.8, usd = sqrt(2), null = 0, pm = 0.5,
          psd = 0, dpm = 0.5, dpsd = 0, integer = FALSE,
          analytical = TRUE),
    nbf_expected,
    tolerance = 1e-8,
    info = "nbf01 analytical point-prior sample size agrees with reference root"
)

## Tests powerbf01 wrappers around the same pbf01/nbf01 references; the wrapper
## is introduced in the package-illustration section (bfssd.Rnw 1932-1966).
power_z <- powerbf01(n = 80, k = 1/6, sd = 1, null = 0, pm = 0.3,
                     psd = 0.2, dpm = 0.35, dpsd = 0.05,
                     type = "two.sample")
expect_equal(power_z$power, pbf_expected, tolerance = 1e-10,
             info = "powerbf01 fixed-n wrapper agrees with pbf01 reference")
ssd_z <- powerbf01(power = 0.8, k = 1/10, sd = 1, null = 0, pm = 0.5,
                   psd = 0, dpm = 0.5, dpsd = 0, type = "two.sample",
                   nrange = c(1, 200))
expect_equal(ssd_z$n, nbf_expected, tolerance = 1e-8,
             info = "powerbf01 sample-size wrapper agrees with nbf01 reference")

## Tests nmbf01 value/log output from the normal-moment BF formula nlBF
## (bfssd.Rnw 1672-1688). The helper computes the marginal likelihood directly.
expect_equal(
    nmbf01(estimate = 0.25, se = 0.12, null = 0, psd = 0.5/sqrt(2)),
    exp(ref_log_nmbf01(estimate = 0.25, se = 0.12, null = 0,
                       psd = 0.5/sqrt(2))),
    tolerance = 1e-12,
    info = "nmbf01 agrees with the convolved normal-moment marginal likelihood"
)
expect_equal(
    nmbf01(estimate = 0.25, se = 0.12, null = 0, psd = 0.5/sqrt(2),
           log = TRUE),
    ref_log_nmbf01(estimate = 0.25, se = 0.12, null = 0,
                   psd = 0.5/sqrt(2)),
    tolerance = 1e-12,
    info = "nmbf01 log output agrees with the log marginal-likelihood ratio"
)

## Tests pnmbf01, nnmbf01, and powernmbf01 from the normal-moment power formula
## pnlBF (bfssd.Rnw 1693-1709; appendix.Rnw 237-259) and the normal-moment
## example chunk (bfssd.Rnw 1713-1778).
pnm_expected <- ref_pnmbf01(k = 1/6, n = 90, usd = sqrt(2), null = 0,
                            psd = 0.5/sqrt(2), dpm = 0.4, dpsd = 0.1)
expect_equal(
    pnmbf01(k = 1/6, n = 90, usd = sqrt(2), null = 0,
            psd = 0.5/sqrt(2), dpm = 0.4, dpsd = 0.1),
    pnm_expected,
    tolerance = 1e-10,
    info = "pnmbf01 agrees with root-partitioned normal predictive integration"
)
nnm_expected <- stats::uniroot(
    function(n) {
        ref_pnmbf01(k = 1/6, n = n, usd = sqrt(2), null = 0,
                    psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0) - 0.7
    },
    interval = c(1, 500)
)$root
expect_equal(
    nnmbf01(k = 1/6, power = 0.7, usd = sqrt(2), null = 0,
            psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0, integer = FALSE),
    nnm_expected,
    tolerance = 1e-6,
    info = "nnmbf01 sample size agrees with reference root"
)
power_nm <- powernmbf01(n = 90, k = 1/6, sd = 1, null = 0,
                        psd = 0.5/sqrt(2), dpm = 0.4, dpsd = 0.1,
                        type = "two.sample")
expect_equal(power_nm$power, pnm_expected, tolerance = 1e-10,
             info = "powernmbf01 fixed-n wrapper agrees with pnmbf01 reference")
ssd_nm <- powernmbf01(power = 0.7, k = 1/6, sd = 1, null = 0,
                      psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0,
                      type = "two.sample", nrange = c(1, 500))
expect_equal(ssd_nm$n, nnm_expected, tolerance = 1e-6,
             info = "powernmbf01 sample-size wrapper agrees with nnmbf01 reference")

## Tests informed/JZS t BF value, power, sample-size, and wrapper behavior from
## the tBF section and one-sided example (bfssd.Rnw 1481-1530 and 1609-1627).
tbf_expected <- ref_log_tbf01(t = 1.4, n1 = 35, n2 = 40, plocation = 0,
                              pscale = 1/sqrt(2), pdf = 1,
                              type = "two.sample",
                              alternative = "two.sided")
expect_equal(
    tbf01(t = 1.4, n1 = 35, n2 = 40, plocation = 0, pscale = 1/sqrt(2),
          pdf = 1, type = "two.sample", alternative = "two.sided",
          log = TRUE),
    tbf_expected,
    tolerance = 1e-7,
    info = "tbf01 agrees with direct noncentral-t prior integration"
)
t_power_expected <- ref_ptbf01_greater(k = 1/6, n = 40, null = 0,
                                       plocation = 0, pscale = 1/sqrt(2),
                                       pdf = 1, dpm = 0.5, dpsd = 0.05)
expect_equal(
    ptbf01(k = 1/6, n = 40, null = 0, plocation = 0, pscale = 1/sqrt(2),
           pdf = 1, dpm = 0.5, dpsd = 0.05, alternative = "greater",
           type = "two.sample"),
    t_power_expected,
    tolerance = 1e-6,
    info = "ptbf01 one-sided power agrees with direct t-boundary integration"
)
t_n_expected <- stats::uniroot(
    function(n) {
        ref_ptbf01_greater(k = 1/6, n = n, null = 0, plocation = 0,
                           pscale = 1/sqrt(2), pdf = 1, dpm = 0.5,
                           dpsd = 0.05) - 0.5
    },
    interval = c(10, 100)
)$root
expect_equal(
    ntbf01(k = 1/6, power = 0.5, null = 0, plocation = 0,
           pscale = 1/sqrt(2), pdf = 1, dpm = 0.5, dpsd = 0.05,
           alternative = "greater", type = "two.sample",
           integer = FALSE, nrange = c(10, 100)),
    t_n_expected,
    tolerance = 1e-4,
    info = "ntbf01 sample size agrees with reference ptbf01 root"
)
power_t <- powertbf01(n = 40, k = 1/6, null = 0, plocation = 0,
                      pscale = 1/sqrt(2), pdf = 1, dpm = 0.5,
                      dpsd = 0.05, alternative = "greater",
                      type = "two.sample")
expect_equal(power_t$power, t_power_expected, tolerance = 1e-6,
             info = "powertbf01 fixed-n wrapper agrees with ptbf01 reference")
ssd_t <- powertbf01(power = 0.5, k = 1/6, null = 0, plocation = 0,
                    pscale = 1/sqrt(2), pdf = 1, dpm = 0.5,
                    dpsd = 0.05, alternative = "greater",
                    type = "two.sample", nrange = c(10, 100))
expect_equal(ssd_t$n, t_n_expected, tolerance = 1e-4,
             info = "powertbf01 sample-size wrapper agrees with ntbf01 reference")

## Package-only coverage: binomial BF APIs are not derived in bfssd.Rnw; binary
## outcomes are explicitly listed as future work (bfssd.Rnw 1854-1856).
expect_equal(
    binbf01(x = 17, n = 25, p0 = 0.5, type = "point", a = 2, b = 3),
    exp(ref_log_binbf01(x = 17, n = 25, p0 = 0.5, type = "point",
                        a = 2, b = 3)),
    tolerance = 1e-12,
    info = "binbf01 point-null value agrees with the beta-binomial formula"
)
expect_equal(
    binbf01(x = 17, n = 25, p0 = 0.5, type = "direction", a = 2, b = 3,
            log = TRUE),
    ref_log_binbf01(x = 17, n = 25, p0 = 0.5, type = "direction",
                    a = 2, b = 3),
    tolerance = 1e-12,
    info = "binbf01 directional log value agrees with posterior/prior odds"
)
bin_p_expected <- ref_pbinbf01(k = 1/10, n = 35, p0 = 0.5, type = "point",
                               a = 2, b = 3, dp = 0.65)
expect_equal(
    pbinbf01(k = 1/10, n = 35, p0 = 0.5, type = "point", a = 2, b = 3,
             dp = 0.65),
    bin_p_expected,
    tolerance = 1e-12,
    info = "pbinbf01 point-design probability agrees with exact enumeration"
)
bin_beta_expected <- ref_pbinbf01(k = 1/10, n = 35, p0 = 0.5,
                                  type = "direction", a = 1, b = 1,
                                  da = 1, db = 1, dl = 0.5, du = 1)
expect_equal(
    pbinbf01(k = 1/10, n = 35, p0 = 0.5, type = "direction", a = 1,
             b = 1, da = 1, db = 1, dl = 0.5, du = 1),
    bin_beta_expected,
    tolerance = 1e-12,
    info = "pbinbf01 truncated-beta probability agrees with exact enumeration"
)
bin_n_expected <- ref_stable_n(
    function(n) ref_pbinbf01(k = 1/10, n = n, p0 = 0.5,
                             type = "direction", a = 1, b = 1, dp = 0.65),
    power = 0.7,
    nrange = c(1, 200)
)
expect_equal(
    nbinbf01(k = 1/10, power = 0.7, p0 = 0.5, type = "direction",
             a = 1, b = 1, dp = 0.65, nrange = c(1, 200)),
    bin_n_expected,
    info = "nbinbf01 sample size agrees with brute-force exact enumeration"
)
power_bin <- powerbinbf01(n = 35, k = 1/10, p0 = 0.5, type = "point",
                          a = 2, b = 3, dp = 0.65)
expect_equal(power_bin$power, bin_p_expected, tolerance = 1e-12,
             info = "powerbinbf01 fixed-n wrapper agrees with pbinbf01 reference")
ssd_bin <- powerbinbf01(power = 0.7, k = 1/10, p0 = 0.5,
                        type = "direction", a = 1, b = 1, dp = 0.65,
                        nrange = c(1, 200))
expect_equal(ssd_bin$n, bin_n_expected,
             info = "powerbinbf01 sample-size wrapper agrees with nbinbf01 reference")

## Tests literal mirtazapine case-study calculations from the manuscript chunks:
## setup/BF at bfssd.Rnw 1060-1077 and 1080-1116, design/sample-size at
## bfssd.Rnw 1121-1196 and 1209-1224.
mirt_est <- -1.74
mirt_ci <- c(-7.17, 3.69)
mirt_se <- (mirt_ci[2] - mirt_ci[1])/(2*stats::qnorm(0.975))
mirt_pm <- -6
mirt_usd <- sqrt(2)*15
mirt_n_expected <- ceiling(
    (stats::qnorm(0.8) +
     sqrt(stats::qnorm(0.8)^2 - log((1/10)^2)*
          (mirt_pm + 0 - 2*mirt_pm)/(0 - mirt_pm)))^2/
        (mirt_pm + 0 - 2*mirt_pm)^2*mirt_usd^2
)
expect_equal(
    bf01(estimate = mirt_est, se = mirt_se, null = 0, pm = mirt_pm,
         psd = 0),
    exp(ref_log_bf01(estimate = mirt_est, se = mirt_se, null = 0,
                     pm = mirt_pm, psd = 0)),
    tolerance = 1e-12,
    info = "mirtazapine manuscript BF01 agrees with density-ratio reference"
)
expect_equal(
    ceiling(nbf01(k = 1/10, power = 0.8, usd = mirt_usd, null = 0,
                  pm = mirt_pm, psd = 0, dpm = mirt_pm, dpsd = 0)),
    mirt_n_expected,
    info = "mirtazapine manuscript sample size agrees with closed-form reference"
)
