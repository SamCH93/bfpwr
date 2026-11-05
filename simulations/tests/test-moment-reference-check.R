## Run from the repository root with ANALYSIS_RSCRIPT.
source("simulations/scripts/simulation_support.R")
source_simulation_library(getwd())
source_package_checkout(getwd())
source("simulations/R/moment-extension.R")
source("simulations/R/moment-reference-check.R")

check <- function(alternative, null = .2, usd = 1, dpm = .35, dpsd = .1, psd = .3) {
    n <- c(10, 30)
    boundary <- function(k) vapply(n, function(ni) zcrit(k, se = usd/sqrt(ni),
        null = null, tau = psd, type = "moment", alternative = alternative), numeric(1))
    b1 <- boundary(1/10)
    b0 <- boundary(10)
    lo <- pmin(b1, b0)
    hi <- pmax(b1, b0)
    r <- bfpwr_sim_moment_markov(n, lo, hi, dpm, dpsd, null, usd,
                                bfpwr_sim_moment_legendre(128L))
    mu <- (dpm - null)*sqrt(n)/usd
    sd1 <- sqrt(1 + n[1]*dpsd^2/usd^2)
    a <- (usd^2 + n[2]*dpsd^2)/(usd^2 + n[1]*dpsd^2)
    cm <- function(z) mu[2] + a*sqrt(n[1]/n[2])*(z - mu[1])
    cs <- sqrt(a*(n[2] - n[1])/n[2])
    conditional <- function(event) integrate(function(z) {
        probability <- switch(event,
            lower = pnorm(lo[2], cm(z), cs),
            upper = pnorm(hi[2], cm(z), cs, lower.tail = FALSE),
            continuation = pnorm(hi[2], cm(z), cs) - pnorm(lo[2], cm(z), cs))
        dnorm(z, mu[1], sd1)*probability
    }, lo[1], hi[1], abs.tol = 1e-12, rel.tol = 1e-12)$value
    stopifnot(abs(r$lower_stop[1] - pnorm(lo[1], mu[1], sd1)) < 1e-12,
        abs(r$upper_stop[1] - pnorm(hi[1], mu[1], sd1, lower.tail = FALSE)) < 1e-12,
        abs(r$lower_stop[2] - conditional("lower")) < 1e-10,
        abs(r$upper_stop[2] - conditional("upper")) < 1e-10,
        abs(r$continuation[2] - conditional("continuation")) < 1e-10,
        max(r$mass_error) < 1e-10)
    c(bfpwr_sim_moment_prefix(r, n, 2L), lower = sum(r$lower_stop),
      upper = sum(r$upper_stop), continuation = r$continuation[2])
}
for (alternative in c("less", "greater")) {
    baseline <- check(alternative)
    translated <- check(alternative, null = 2.2, dpm = 2.35)
    scaled <- check(alternative, null = 2, usd = 10, dpm = 3.5, dpsd = 1, psd = 3)
    stopifnot(max(abs(baseline - translated)) < 1e-8,
              max(abs(baseline - scaled)) < 1e-8)
    check(alternative, dpsd = 0)
}
less <- check("less", null = -.2, dpm = -.35)
greater <- check("greater")
stopifnot(max(abs(less[c("EN", "VarN", "continuation")] -
                     greater[c("EN", "VarN", "continuation")])) < 1e-8,
          abs(less["lower"] - greater["upper"]) < 1e-10,
          abs(less["upper"] - greater["lower"]) < 1e-10)
cat("PASS: conditional short-look probabilities, mass, both directions, point/normal designs, translation, scale, reflection\n")

## Independent joint-MVN rectangles use predpars, without the Markov transition.
integration_refinements <- numeric()
for (alternative in c("less", "greater")) {
    n <- c(10, 30, 60)
    null <- .2
    usd <- 1.3
    dpm <- .35
    dpsd <- .1
    boundary <- function(k) vapply(n, function(ni) zcrit(k, se = usd/sqrt(ni),
        null = null, tau = .3, type = "moment", alternative = alternative), numeric(1))
    lo <- pmin(boundary(.1), boundary(10))
    hi <- pmax(boundary(.1), boundary(10))
    joint <- predpars(usd/sqrt(n), null = null, dpm = dpm, dpsd = dpsd)
    probability <- function(lower, upper) {
        i <- seq_along(lower)
        value <- mvtnorm::pmvnorm(lower, upper, mean = joint$mean[i],
            sigma = joint$sigma[i, i, drop = FALSE],
            algorithm = mvtnorm::Miwa(steps = 512L))
        coarse <- mvtnorm::pmvnorm(lower, upper, mean = joint$mean[i],
            sigma = joint$sigma[i, i, drop = FALSE], algorithm = mvtnorm::Miwa(steps = 256L))
        refinement <- abs(as.numeric(value) - as.numeric(coarse))
        integration_refinements <<- c(integration_refinements, refinement)
        stopifnot(refinement <= 1e-9)
        as.numeric(value)
    }
    lower_stop <- upper_stop <- continuation <- numeric(3L)
    for (i in seq_along(n)) {
        previous <- seq_len(i - 1L)
        lower_stop[i] <- probability(c(lo[previous], -Inf), c(hi[previous], lo[i]))
        upper_stop[i] <- probability(c(lo[previous], hi[i]), c(hi[previous], Inf))
        continuation[i] <- probability(lo[seq_len(i)], hi[seq_len(i)])
    }
    actual <- bfpwr_sim_moment_markov(n, lo, hi, dpm, dpsd, null, usd,
        bfpwr_sim_moment_legendre(128L))
    ## Ten times the checked 256/512-step MVN refinement, including accumulated mass.
    tolerance <- 1e-8
    stopifnot(max(abs(actual$lower_stop - lower_stop)) < tolerance,
        max(abs(actual$upper_stop - upper_stop)) < tolerance,
        max(abs(actual$continuation - continuation)) < tolerance)
    mass <- lower_stop + upper_stop
    mass[3L] <- mass[3L] + continuation[3L]
    expected_en <- sum(n*mass)
    expected_variance <- sum(n^2*mass) - expected_en^2
    moments <- bfpwr_sim_moment_prefix(actual, n, 3L)
    ## Propagate the per-probability tolerance through the two raw moments.
    en_tolerance <- tolerance*sum(n)
    variance_tolerance <- tolerance*sum(n^2) + 2*max(n)*en_tolerance
    stopifnot(abs(moments["EN"] - expected_en) < en_tolerance,
        abs(moments["VarN"] - expected_variance) < variance_tolerance)
}
cat("PASS: independent three-look predpars/pmvnorm probabilities and stopping moments, both directions; maximum MVN refinement:",
    max(integration_refinements), "\n")

args <- commandArgs(TRUE)
if (length(args)) {
    bundle <- readRDS(args[1L])
    start <- proc.time()[["elapsed"]]
    diagnostic <- bfpwr_sim_moment_reference_check(bundle)
    print(diagnostic$summary)
    stopifnot(nrow(diagnostic$moment_rows) == nrow(bundle$moments),
        !anyNA(diagnostic$moment_rows$converged_EN), diagnostic$summary$passed,
        identical(diagnostic$moment_rows[c("case_id", "schedule_id", "evidence_threshold")],
                  bundle$moments[c("case_id", "schedule_id", "evidence_threshold")]))
    small <- bundle
    small$moments <- diagnostic$moment_rows[1L, names(bundle$moments), drop=FALSE]
    small$moments$max_observed_n <- 30
    small$moments$tolerance_EN <- 123
    small$moments$tolerance_VarN <- 456
    small$designs <- bfpwr_sim_design_cases_z()
    recorded <- bfpwr_sim_moment_reference_check(small)
    stopifnot(recorded$moment_rows$tolerance_EN == 123,
              recorded$moment_rows$tolerance_VarN == 456,
              recorded$summary$tolerance_source == "recorded bundle tolerances",
              recorded$summary$design_source == "frozen bundle design specifications")
    missing_frozen <- small
    missing_frozen$designs <- NULL
    missing_frozen$timing_schema_version <- 2L
    stopifnot(inherits(try(bfpwr_sim_moment_reference_check(missing_frozen), silent=TRUE), "try-error"))
    original_quadrature <- bfpwr_sim_moment_legendre
    bfpwr_sim_moment_legendre <- function(nodes) {
        q <- original_quadrature(nodes)
        if (nodes == 64L) q$weights <- q$weights*1.0001
        q
    }
    failed <- bfpwr_sim_moment_reference_check(small)
    bfpwr_sim_moment_legendre <- original_quadrature
    stopifnot(!failed$summary$passed, !failed$moment_rows$reference_converged)
    small$cases[[1L]]$analysis_prior$alternative <- "two.sided"
    stopifnot(inherits(try(bfpwr_sim_moment_reference_check(small), silent=TRUE), "try-error"))
    cat("PASS: frozen designs, recorded tolerances preserved, convergence failure visible, unsupported case rejected\n")
    cat("PASS: all", nrow(bundle$moments), "bundle rows; seconds:", proc.time()[["elapsed"]] - start, "\n")
    if (length(args) > 1L) saveRDS(diagnostic, args[2L])
}
