library(tinytest)
library(bfpwr)

## Strict sequential checks. The helper oracles construct critical values and
## multivariate-normal stopping regions independently of bfpwr. Only the Low-PV
## block below is a direct external BFGSD paper example; the other blocks are
## package-level adversarial checks for region algebra, signs, and recentering.

## Independent MVN/root/integration helpers. These support the tests below but
## have no direct manuscript source.
pmv <- function(lower, upper, mean, sigma, ...) {
    as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean,
                                sigma = sigma, keepAttr = FALSE, ...))
}

manual_predpars <- function(se, dpm, dpsd) {
    information <- 1/se^2
    sigma <- outer(information, information,
                   function(a, b) sqrt(pmin(a, b)/pmax(a, b)))
    sigma <- sigma + dpsd^2 * sqrt(information) %*% t(sqrt(information))
    list(mean = dpm/se, sigma = sigma)
}

manual_zcrit_point <- function(k, se, mu) {
    (mu^2/se^2 - 2*log(k))/(2*mu/se)
}

manual_log_bf01_z <- function(z, se, pm, psd) {
    estimate <- z*se
    stats::dnorm(estimate, mean = 0, sd = se, log = TRUE) -
        stats::dnorm(estimate, mean = pm, sd = sqrt(se^2 + psd^2),
                     log = TRUE)
}

manual_zcrit_normal_roots <- function(k, se, pm, psd) {
    f <- function(z) manual_log_bf01_z(z, se = se, pm = pm, psd = psd) -
        log(k)
    width <- 4 + abs(pm/se) + if (psd > 0) 6*psd/se else 0
    for (multiplier in c(1, 2, 4, 8)) {
        grid <- seq(-width*multiplier, width*multiplier, length.out = 4001)
        values <- f(grid)
        roots <- numeric()
        for (i in seq_len(length(grid) - 1)) {
            f0 <- values[i]
            f1 <- values[i + 1]
            if (!is.finite(f0) || !is.finite(f1)) next
            if (f0 == 0) {
                roots <- c(roots, grid[i])
            } else if (f0*f1 < 0) {
                roots <- c(roots, stats::uniroot(f, grid[c(i, i + 1)],
                                                 tol = 1e-12)$root)
            }
        }
        roots <- sort(unique(round(roots, 12)))
        if (length(roots) > 0) return(roots)
    }
    numeric()
}

manual_zcrit_normal_one <- function(k, se, pm) {
    roots <- manual_zcrit_normal_roots(k = k, se = se, pm = pm, psd = 0)
    stopifnot(length(roots) == 1)
    roots
}

manual_zcrit_normal_two <- function(k, se, pm, psd) {
    roots <- manual_zcrit_normal_roots(k = k, se = se, pm = pm, psd = psd)
    if (length(roots) == 0) return(c(NaN, NaN))
    stopifnot(length(roots) == 2)
    roots
}

manual_zcrit_directional <- function(k, se, mu, tau) {
    logpriorodds <- stats::pnorm(mu/tau, lower.tail = FALSE, log.p = TRUE) -
        stats::pnorm(mu/tau, log.p = TRUE)
    postq <- stats::qnorm(1/(1 + exp(log(k) + logpriorodds)))
    (postq*sqrt(1/se^2 + 1/tau^2) - mu/tau^2)*se
}

manual_zcrit_moment <- function(k, se, tau) {
    y <- (2*lamW::lambertW0(exp(0.5)*(1 + tau^2/se^2)^1.5/(2*k)) - 1)*
        (1 + se^2/tau^2)
    c(-1, 1)*sqrt(y)
}

manual_one_boundary_seq <- function(zcrit0, zcrit1, se, n, dpm, dpsd,
                                    positive = TRUE) {
    pars <- manual_predpars(se, dpm, dpsd)
    mean <- pars$mean
    sigma <- pars$sigma
    if (positive) {
        pH1 <- c(
            pmv(zcrit1[1], Inf, mean[1], matrix(sigma[1, 1], 1)),
            pmv(c(zcrit0[1], zcrit1[2]), c(zcrit1[1], Inf),
                mean, sigma)
        )
        pH0 <- c(
            pmv(-Inf, zcrit0[1], mean[1], matrix(sigma[1, 1], 1)),
            pmv(c(zcrit0[1], -Inf), c(zcrit1[1], zcrit0[2]),
                mean, sigma)
        )
    } else {
        pH1 <- c(
            pmv(-Inf, zcrit1[1], mean[1], matrix(sigma[1, 1], 1)),
            pmv(c(zcrit1[1], -Inf), c(zcrit0[1], zcrit1[2]),
                mean, sigma)
        )
        pH0 <- c(
            pmv(zcrit0[1], Inf, mean[1], matrix(sigma[1, 1], 1)),
            pmv(c(zcrit1[1], zcrit0[2]), c(zcrit0[1], Inf),
                mean, sigma)
        )
    }
    list(cumpH1 = cumsum(pH1),
         cumpH0 = cumsum(pH0),
         cumpInc = 1 - cumsum(pH1) - cumsum(pH0),
         EN = sum((pH1 + pH0)*n) + (1 - sum(pH1 + pH0))*max(n))
}

manual_one_boundary_seq_any <- function(zcrit0, zcrit1, se, n, dpm, dpsd,
                                        positive = TRUE, ...) {
    pars <- manual_predpars(se, dpm, dpsd)
    mean <- pars$mean
    sigma <- pars$sigma
    m <- length(n)
    stage_prob <- function(i, stop_for) {
        lower <- numeric(i)
        upper <- numeric(i)
        for (j in seq_len(i)) {
            if (j < i) {
                if (positive) {
                    lower[j] <- if (is.nan(zcrit0[j])) -Inf else zcrit0[j]
                    upper[j] <- zcrit1[j]
                } else {
                    lower[j] <- zcrit1[j]
                    upper[j] <- if (is.nan(zcrit0[j])) Inf else zcrit0[j]
                }
            } else if (identical(stop_for, "H1")) {
                if (positive) {
                    lower[j] <- zcrit1[j]
                    upper[j] <- Inf
                } else {
                    lower[j] <- -Inf
                    upper[j] <- zcrit1[j]
                }
            } else {
                if (positive) {
                    lower[j] <- -Inf
                    upper[j] <- zcrit0[j]
                } else {
                    lower[j] <- zcrit0[j]
                    upper[j] <- Inf
                }
            }
        }
        pmv(lower, upper, mean[seq_len(i)], sigma[seq_len(i), seq_len(i),
                                                  drop = FALSE], ...)
    }
    pH1 <- vapply(seq_len(m), stage_prob, numeric(1), stop_for = "H1")
    pH0 <- vapply(seq_len(m), stage_prob, numeric(1), stop_for = "H0")
    list(cumpH1 = cumsum(pH1),
         cumpH0 = cumsum(pH0),
         cumpInc = 1 - cumsum(pH1) - cumsum(pH0),
         EN = sum((pH1 + pH0)*n) + (1 - sum(pH1 + pH0))*max(n))
}

manual_two_boundary_seq <- function(zcrit0, zcrit1, se, n, dpm, dpsd) {
    pars <- manual_predpars(se, dpm, dpsd)
    mean <- pars$mean
    sigma <- pars$sigma
    lower_h1 <- zcrit1[1, ]
    upper_h1 <- zcrit1[2, ]
    lower_h0 <- zcrit0[1, ]
    upper_h0 <- zcrit0[2, ]
    pH1 <- c(
        pmv(-Inf, lower_h1[1], mean[1], matrix(sigma[1, 1], 1)) +
            pmv(upper_h1[1], Inf, mean[1], matrix(sigma[1, 1], 1)),
        pmv(c(lower_h1[1], -Inf), c(lower_h0[1], lower_h1[2]),
            mean, sigma) +
            pmv(c(lower_h1[1], upper_h1[2]), c(lower_h0[1], Inf),
                mean, sigma) +
            pmv(c(upper_h0[1], -Inf), c(upper_h1[1], lower_h1[2]),
                mean, sigma) +
            pmv(c(upper_h0[1], upper_h1[2]), c(upper_h1[1], Inf),
                mean, sigma)
    )
    pH0 <- c(
        pmv(lower_h0[1], upper_h0[1], mean[1], matrix(sigma[1, 1], 1)),
        pmv(c(lower_h1[1], lower_h0[2]), c(lower_h0[1], upper_h0[2]),
            mean, sigma) +
            pmv(c(upper_h0[1], lower_h0[2]), c(upper_h1[1], upper_h0[2]),
                mean, sigma)
    )
    list(cumpH1 = cumsum(pH1),
         cumpH0 = cumsum(pH0),
         cumpInc = 1 - cumsum(pH1) - cumsum(pH0),
         EN = sum((pH1 + pH0)*n) + (1 - sum(pH1 + pH0))*max(n))
}

manual_two_boundary_seq_any <- function(zcrit0, zcrit1, se, n, dpm, dpsd,
                                        ...) {
    pars <- manual_predpars(se, dpm, dpsd)
    mean <- pars$mean
    sigma <- pars$sigma
    m <- length(n)
    stage_prob <- function(i, stop_for) {
        intervals <- vector("list", i)
        for (j in seq_len(i)) {
            h0_missing <- any(is.nan(zcrit0[, j]))
            if (j < i) {
                intervals[[j]] <- if (h0_missing) {
                    list(c(zcrit1[1, j], zcrit1[2, j]))
                } else {
                    list(c(zcrit1[1, j], zcrit0[1, j]),
                         c(zcrit0[2, j], zcrit1[2, j]))
                }
            } else if (identical(stop_for, "H1")) {
                intervals[[j]] <- list(c(-Inf, zcrit1[1, j]),
                                       c(zcrit1[2, j], Inf))
            } else {
                intervals[[j]] <- list(c(zcrit0[1, j], zcrit0[2, j]))
            }
        }
        combos <- do.call(expand.grid, c(lapply(intervals, seq_along),
                                         KEEP.OUT.ATTRS = FALSE))
        sum(apply(combos, 1, function(row) {
            row <- as.integer(row)
            bounds <- vapply(seq_along(row), function(j) {
                intervals[[j]][[row[j]]]
            }, numeric(2))
            if (any(is.nan(bounds))) return(0)
            pmv(bounds[1, ], bounds[2, ], mean[seq_len(i)],
                sigma[seq_len(i), seq_len(i), drop = FALSE], ...)
        }))
    }
    pH1 <- vapply(seq_len(m), stage_prob, numeric(1), stop_for = "H1")
    pH0 <- vapply(seq_len(m), stage_prob, numeric(1), stop_for = "H0")
    pstop <- pH1 + pH0
    EN <- sum(pstop*n) + (1 - sum(pstop))*max(n)
    EN2 <- sum(pstop*n^2) + (1 - sum(pstop))*max(n)^2
    list(cumpH1 = cumsum(pH1),
         cumpH0 = cumsum(pH0),
         cumpInc = 1 - cumsum(pH1) - cumsum(pH0),
         EN = EN,
         VarN = EN2 - EN^2)
}

manual_log_tbf01 <- function(t, n1, n2, plocation, pscale, pdf,
                             alternative = "greater") {
    df <- n1 + n2 - 2
    neff <- 1/(1/n1 + 1/n2)
    lower <- if (alternative == "greater") 0 else -Inf
    upper <- if (alternative == "greater") Inf else 0
    norm_const <- if (alternative == "greater") {
        stats::pt((0 - plocation)/pscale, df = pdf, lower.tail = FALSE)
    } else {
        stats::pt((0 - plocation)/pscale, df = pdf)
    }
    prior <- function(delta) {
        stats::dt((delta - plocation)/pscale, df = pdf)/
            (pscale*norm_const)
    }
    marginal <- stats::integrate(
        function(delta) {
            suppressWarnings(stats::dt(t, df = df, ncp = sqrt(neff)*delta))*
                prior(delta)
        },
        lower = lower,
        upper = upper,
        rel.tol = 1e-10,
        subdivisions = 2000
    )$value
    stats::dt(t, df = df, log = TRUE) - log(marginal)
}

manual_tcrit_greater <- function(k, n) {
    stats::uniroot(
        function(t) {
            manual_log_tbf01(t = t, n1 = n, n2 = n, plocation = 0,
                             pscale = 1/sqrt(2), pdf = 1,
                             alternative = "greater") - log(k)
        },
        interval = c(0, 10),
        tol = 1e-8
    )$root
}

normal_n <- c(40, 80)
normal_se <- 1/sqrt(normal_n)
normal_z1 <- manual_zcrit_point(k = 1/5, se = normal_se, mu = 0.25)
normal_z0 <- manual_zcrit_point(k = 4, se = normal_se, mu = 0.25)
normal_ref <- manual_one_boundary_seq(zcrit0 = normal_z0, zcrit1 = normal_z1,
                                      se = normal_se, n = normal_n,
                                      dpm = 0.2, dpsd = 0.05,
                                      positive = TRUE)
normal_pkg <- pbf01seq(k1 = 1/5, k0 = 4, se = normal_se, n = normal_n,
                       pm = 0.25, psd = 0, dpm = 0.2, dpsd = 0.05,
                       type = "normal", strict = TRUE, method = "pmvnorm")
## Package-only adversarial test: two-look point-normal pbf01seq should match
## manually assembled MVN stopping regions. There is no exact manuscript row.
expect_equal(normal_pkg$cumpH1, normal_ref$cumpH1, tolerance = 1e-10,
             info = "two-look point-normal sequential H1 matches manual MVN regions")
expect_equal(normal_pkg$cumpH0, normal_ref$cumpH0, tolerance = 1e-10,
             info = "two-look point-normal sequential H0 matches manual MVN regions")
expect_equal(normal_pkg$EN, normal_ref$EN, tolerance = 1e-10,
             info = "two-look point-normal expected sample size matches manual MVN regions")

directional_n <- c(50, 100)
directional_se <- 1/sqrt(directional_n)
directional_z1 <- manual_zcrit_directional(k = 1/10, se = directional_se,
                                           mu = 0, tau = 1/sqrt(2))
directional_z0 <- manual_zcrit_directional(k = 3, se = directional_se,
                                           mu = 0, tau = 1/sqrt(2))
directional_ref <- manual_one_boundary_seq(zcrit0 = directional_z0,
                                           zcrit1 = directional_z1,
                                           se = directional_se,
                                           n = directional_n,
                                           dpm = 0.25, dpsd = 0.1,
                                           positive = TRUE)
directional_pkg <- pbf01seq(k1 = 1/10, k0 = 3, se = directional_se,
                            n = directional_n, pm = 0, psd = 1/sqrt(2),
                            dpm = 0.25, dpsd = 0.1, type = "directional",
                            strict = TRUE, method = "pmvnorm")
## Package-only adversarial test: directional pbf01seq boundaries and stopping
## regions. Directional normal BF tests have no bfssd/BFGSD manuscript row.
expect_equal(directional_pkg$cumpH1, directional_ref$cumpH1,
             tolerance = 1e-10,
             info = "two-look directional sequential H1 matches manual MVN regions")
expect_equal(directional_pkg$cumpH0, directional_ref$cumpH0,
             tolerance = 1e-10,
             info = "two-look directional sequential H0 matches manual MVN regions")
expect_equal(directional_pkg$EN, directional_ref$EN, tolerance = 1e-10,
             info = "two-look directional expected sample size matches manual MVN regions")

moment_n <- c(45, 90)
moment_se <- 1/sqrt(moment_n)
moment_z1 <- sapply(moment_se, function(se) {
    manual_zcrit_moment(k = 1/8, se = se, tau = 0.5)
})
moment_z0 <- sapply(moment_se, function(se) {
    manual_zcrit_moment(k = 5, se = se, tau = 0.5)
})
moment_ref <- manual_two_boundary_seq(zcrit0 = moment_z0, zcrit1 = moment_z1,
                                      se = moment_se, n = moment_n,
                                      dpm = 0.35, dpsd = 0.05)
moment_pkg <- pbf01seq(k1 = 1/8, k0 = 5, se = moment_se, n = moment_n,
                       psd = 0.5, dpm = 0.35, dpsd = 0.05,
                       type = "moment", strict = TRUE, method = "pmvnorm")
## Package-only sequential extension of the fixed-sample normal-moment formulas
## in bfssd.Rnw 1672-1709 and appendix.Rnw 237-259.
expect_equal(moment_pkg$cumpH1, moment_ref$cumpH1, tolerance = 1e-10,
             info = "two-look moment sequential H1 matches manual MVN regions")
expect_equal(moment_pkg$cumpH0, moment_ref$cumpH0, tolerance = 1e-10,
             info = "two-look moment sequential H0 matches manual MVN regions")
expect_equal(moment_pkg$EN, moment_ref$EN, tolerance = 1e-10,
             info = "two-look moment expected sample size matches manual MVN regions")

normal2_n <- c(50, 100, 150)
normal2_se <- sqrt(2/normal2_n)
normal2_z1 <- sapply(normal2_se, function(se) {
    manual_zcrit_normal_two(k = 1/6, se = se, pm = 0, psd = 0.5)
})
normal2_z0 <- sapply(normal2_se, function(se) {
    manual_zcrit_normal_two(k = 3, se = se, pm = 0, psd = 0.5)
})
normal2_ref <- manual_two_boundary_seq_any(
    zcrit0 = normal2_z0, zcrit1 = normal2_z1, se = normal2_se,
    n = normal2_n, dpm = 0.35, dpsd = 0.05, algorithm = mvtnorm::Miwa())
normal2_pkg <- pbf01seq(k1 = 1/6, k0 = 3, se = normal2_se, n = normal2_n,
                        pm = 0, psd = 0.5, dpm = 0.35, dpsd = 0.05,
                        type = "normal", strict = TRUE,
                        method = "pmvnorm", algorithm = mvtnorm::Miwa())
## Package-only adversarial test: three-look normal-prior schedules, two-sided
## BF regions, VarN, and the default lpmvnorm path. Related fixed-sample formulas
## are bfssd.Rnw 590-606 and 640-653.
expect_equal(normal2_pkg$cumpH1, normal2_ref$cumpH1, tolerance = 1e-10,
             info = "three-look normal-prior H1 matches root-found MVN regions")
expect_equal(normal2_pkg$cumpH0, normal2_ref$cumpH0, tolerance = 1e-10,
             info = "three-look normal-prior H0 matches root-found MVN regions")
expect_equal(normal2_pkg$EN, normal2_ref$EN, tolerance = 1e-10,
             info = "three-look normal-prior expected n matches root-found MVN regions")
expect_equal(normal2_pkg$VarN, normal2_ref$VarN, tolerance = 1e-10,
             info = "three-look normal-prior sample-size variance matches root-found MVN regions")
normal2_default <- pbf01seq(k1 = 1/6, k0 = 3, se = normal2_se, n = normal2_n,
                            pm = 0, psd = 0.5, dpm = 0.35, dpsd = 0.05,
                            type = "normal", strict = TRUE)
expect_true(max(abs(normal2_default$cumpH1 - normal2_pkg$cumpH1),
                abs(normal2_default$cumpH0 - normal2_pkg$cumpH0),
                abs(normal2_default$EN - normal2_pkg$EN)/max(normal2_n)) <
                5e-4,
            info = "default lpmvnorm path stays close to exact pmvnorm reference")

shifted_null <- 0.2
shifted_pm <- 0.45 - shifted_null
shifted_dpm <- 0.35 - shifted_null
shifted_n <- c(60, 120)
shifted_se <- 1/sqrt(shifted_n)
shifted_z1 <- manual_zcrit_point(k = 1/5, se = shifted_se, mu = shifted_pm)
shifted_z0 <- manual_zcrit_point(k = 4, se = shifted_se, mu = shifted_pm)
shifted_ref <- manual_one_boundary_seq(zcrit0 = shifted_z0, zcrit1 = shifted_z1,
                                       se = shifted_se, n = shifted_n,
                                       dpm = shifted_dpm, dpsd = 0.05,
                                       positive = TRUE)
shifted_pkg <- pbf01seq(k1 = 1/5, k0 = 4, se = shifted_se, n = shifted_n,
                        pm = shifted_pm, psd = 0, dpm = shifted_dpm,
                        dpsd = 0.05, type = "normal", strict = TRUE,
                        method = "pmvnorm")
## Package-only adversarial test: nonzero-null/recentered normal pbf01seq. This
## guards implementation behavior and is not a literal paper example.
expect_equal(shifted_pkg$cumpH1, shifted_ref$cumpH1, tolerance = 1e-10,
             info = "two-look recentered nonzero-null H1 matches manual MVN regions")
expect_equal(shifted_pkg$cumpH0, shifted_ref$cumpH0, tolerance = 1e-10,
             info = "two-look recentered nonzero-null H0 matches manual MVN regions")

negative_n <- c(40, 80)
negative_se <- 1/sqrt(negative_n)
negative_pm <- -0.25
negative_z1 <- vapply(negative_se, function(se) {
    manual_zcrit_normal_one(k = 1/5, se = se, pm = negative_pm)
}, numeric(1))
negative_z0 <- vapply(negative_se, function(se) {
    manual_zcrit_normal_one(k = 4, se = se, pm = negative_pm)
}, numeric(1))
negative_ref <- manual_one_boundary_seq(zcrit0 = negative_z0,
                                        zcrit1 = negative_z1,
                                        se = negative_se, n = negative_n,
                                        dpm = -0.2, dpsd = 0.05,
                                        positive = FALSE)
negative_pkg <- pbf01seq(k1 = 1/5, k0 = 4, se = negative_se, n = negative_n,
                         pm = negative_pm, psd = 0, dpm = -0.2,
                         dpsd = 0.05, type = "normal", strict = TRUE,
                         method = "pmvnorm")
## Package-only adversarial test: negative-effect orientation and sign handling
## for point-normal sequential regions.
expect_equal(negative_pkg$cumpH1, negative_ref$cumpH1, tolerance = 1e-10,
             info = "two-look negative point-normal H1 matches manual MVN regions")
expect_equal(negative_pkg$cumpH0, negative_ref$cumpH0, tolerance = 1e-10,
             info = "two-look negative point-normal H0 matches manual MVN regions")
expect_equal(negative_pkg$EN, negative_ref$EN, tolerance = 1e-10,
             info = "two-look negative point-normal expected sample size matches manual MVN regions")

lowpv_p0 <- 0.5
lowpv_p1 <- 0.75
lowpv_pm <- log((lowpv_p1/(1 - lowpv_p1))/(lowpv_p0/(1 - lowpv_p0)))
lowpv_n <- c(25, 50, 75)
lowpv_se_h1 <- sqrt(1/(lowpv_p0*(1 - lowpv_p0)*lowpv_n) +
                    1/(lowpv_p1*(1 - lowpv_p1)*lowpv_n))
lowpv_se_h0 <- sqrt(1/(lowpv_p0*(1 - lowpv_p0)*lowpv_n) +
                    1/(lowpv_p0*(1 - lowpv_p0)*lowpv_n))
lowpv_z1_h1 <- manual_zcrit_point(k = 1/10, se = lowpv_se_h1,
                                  mu = lowpv_pm)
lowpv_z0_h1 <- manual_zcrit_point(k = 10, se = lowpv_se_h1,
                                  mu = lowpv_pm)
lowpv_z1_h0 <- manual_zcrit_point(k = 1/10, se = lowpv_se_h0,
                                  mu = lowpv_pm)
lowpv_z0_h0 <- manual_zcrit_point(k = 10, se = lowpv_se_h0,
                                  mu = lowpv_pm)
lowpv_ref_h1 <- manual_one_boundary_seq_any(
    zcrit0 = lowpv_z0_h1, zcrit1 = lowpv_z1_h1, se = lowpv_se_h1,
    n = lowpv_n, dpm = lowpv_pm, dpsd = 0, positive = TRUE,
    algorithm = mvtnorm::Miwa())
lowpv_ref_h0 <- manual_one_boundary_seq_any(
    zcrit0 = lowpv_z0_h0, zcrit1 = lowpv_z1_h0, se = lowpv_se_h0,
    n = lowpv_n, dpm = 0, dpsd = 0, positive = TRUE,
    algorithm = mvtnorm::Miwa())
lowpv_pkg_h1 <- pbf01seq(k1 = 1/10, k0 = 10, se = lowpv_se_h1,
                         n = lowpv_n, pm = lowpv_pm, psd = 0,
                         dpm = lowpv_pm, dpsd = 0, type = "normal",
                         method = "pmvnorm", algorithm = mvtnorm::Miwa())
lowpv_pkg_h0 <- pbf01seq(k1 = 1/10, k0 = 10, se = lowpv_se_h0,
                         n = lowpv_n, pm = lowpv_pm, psd = 0,
                         dpm = 0, dpsd = 0, type = "normal",
                         method = "pmvnorm", algorithm = mvtnorm::Miwa())
## Direct external-paper block: Low-PV three-look design from SamCH93/bfgsd
## paper/BFGSD.R 386-413 and paper/BFGSD.Rnw 1005-1025. The test strengthens
## the manuscript calculation by reconstructing MVN regions manually.
expect_equal(lowpv_pkg_h1$cumpH1, lowpv_ref_h1$cumpH1, tolerance = 1e-10,
             info = "Low-PV three-look H1 design matches manual MVN regions")
expect_equal(lowpv_pkg_h1$cumpH0, lowpv_ref_h1$cumpH0, tolerance = 1e-10,
             info = "Low-PV three-look H1 false-stop curve matches manual MVN regions")
expect_equal(lowpv_pkg_h1$EN, lowpv_ref_h1$EN, tolerance = 1e-10,
             info = "Low-PV three-look H1 expected sample size matches manual MVN regions")
expect_equal(lowpv_pkg_h0$cumpH0, lowpv_ref_h0$cumpH0, tolerance = 1e-10,
             info = "Low-PV three-look H0 design matches manual MVN regions")
expect_equal(lowpv_pkg_h0$cumpH1, lowpv_ref_h0$cumpH1, tolerance = 1e-10,
             info = "Low-PV three-look H0 false-stop curve matches manual MVN regions")
expect_equal(lowpv_pkg_h0$EN, lowpv_ref_h0$EN, tolerance = 1e-10,
             info = "Low-PV three-look H0 expected sample size matches manual MVN regions")

t_n <- c(30, 60)
t_se <- sqrt(2/t_n)
t_z1 <- vapply(t_n, function(n) manual_tcrit_greater(k = 1/6, n = n),
               numeric(1))
t_z0 <- vapply(t_n, function(n) manual_tcrit_greater(k = 3, n = n),
               numeric(1))
t_ref <- manual_one_boundary_seq(zcrit0 = t_z0, zcrit1 = t_z1,
                                 se = t_se, n = t_n, dpm = 0.5,
                                 dpsd = 0.05, positive = TRUE)
t_pkg <- ptbf01seq(k1 = 1/6, k0 = 3, n = t_n, plocation = 0,
                   pscale = 1/sqrt(2), pdf = 1, dpm = 0.5, dpsd = 0.05,
                   type = "two.sample", alternative = "greater",
                   strict = TRUE, method = "pmvnorm")
## Package-level reduced JZS sequential check. It is related to the BFGSD
## appendix one-sided JZS design (BFGSD.R 823-827; BFGSD.Rnw 1701-1704), but uses
## smaller thresholds and looks; the exact appendix reproduction is in
## simulations/scripts/verify_paper_values.R.
expect_equal(t_pkg$cumpH1, t_ref$cumpH1, tolerance = 5e-6,
             info = "two-look one-sided t sequential H1 matches manual MVN regions")
expect_equal(t_pkg$cumpH0, t_ref$cumpH0, tolerance = 5e-6,
             info = "two-look one-sided t sequential H0 matches manual MVN regions")
expect_equal(t_pkg$EN1, t_ref$EN, tolerance = 5e-6,
             info = "two-look one-sided t expected sample size matches manual MVN regions")
