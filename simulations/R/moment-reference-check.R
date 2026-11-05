## Independent Gaussian-Markov convergence diagnostic for the one-sided
## normal-moment grid. BF boundaries always come from the package.
bfpwr_sim_moment_legendre <- function(nodes) {
    j <- seq_len(nodes - 1L)
    b <- j/sqrt(4*j^2 - 1)
    jacobi <- matrix(0, nodes, nodes)
    jacobi[cbind(j, j + 1L)] <- b
    jacobi[cbind(j + 1L, j)] <- b
    e <- eigen(jacobi, symmetric = TRUE)
    list(nodes = e$values, weights = 2*e$vectors[1L, ]^2)
}

bfpwr_sim_moment_markov <- function(n, lower, upper, dpm, dpsd, null, usd, quadrature) {
    stopifnot(length(n) == length(lower), length(n) == length(upper),
              all(is.finite(c(n, lower, upper, dpm, dpsd, null, usd))),
              all(n > 0), all(diff(n) > 0), all(lower < upper), usd > 0, dpsd >= 0)
    mu <- (dpm - null)*sqrt(n)/usd
    lower_stop <- upper_stop <- continuation <- numeric(length(n))
    for (i in seq_along(n)) {
        z <- (lower[i] + upper[i])/2 + (upper[i] - lower[i])/2*quadrature$nodes
        weights <- (upper[i] - lower[i])/2*quadrature$weights
        if (i == 1L) {
            initial_sd <- sqrt(1 + n[i]*dpsd^2/usd^2)
            lower_stop[i] <- pnorm(lower[i], mu[i], initial_sd)
            upper_stop[i] <- pnorm(upper[i], mu[i], initial_sd, lower.tail = FALSE)
            density <- dnorm(z, mu[i], initial_sd)
        } else {
            a <- (usd^2 + n[i]*dpsd^2)/(usd^2 + n[i - 1L]*dpsd^2)
            conditional_mean <- mu[i] + a*sqrt(n[i - 1L]/n[i])*(previous_z - mu[i - 1L])
            conditional_sd <- sqrt(a*(n[i] - n[i - 1L])/n[i])
            lower_stop[i] <- sum(previous_mass*pnorm(lower[i], conditional_mean, conditional_sd))
            upper_stop[i] <- sum(previous_mass*pnorm(upper[i], conditional_mean, conditional_sd,
                                                   lower.tail = FALSE))
            density <- as.vector(dnorm(outer(z, conditional_mean, "-"), sd = conditional_sd) %*%
                                 previous_mass)
        }
        previous_z <- z
        previous_mass <- density*weights
        continuation[i] <- sum(previous_mass)
    }
    list(lower_stop = lower_stop, upper_stop = upper_stop, continuation = continuation,
         mass_error = abs(cumsum(lower_stop + upper_stop) + continuation - 1))
}

bfpwr_sim_moment_prefix <- function(result, n, looks) {
    index <- seq_len(looks)
    probability <- result$lower_stop[index] + result$upper_stop[index]
    probability[looks] <- probability[looks] + result$continuation[looks]
    en <- sum(probability*n[index])
    c(EN = en, VarN = sum(probability*(n[index] - en)^2),
      mass_error = max(result$mass_error[index]))
}

bfpwr_sim_moment_reference_check <- function(bundle) {
    m <- bundle$moments
    keys <- c("case_id", "schedule_id", "evidence_threshold")
    required <- c(keys, "start_n", "increment", "n_looks", "EN", "VarN",
                  "reference_EN", "reference_VarN", "mcse_EN", "mcse_VarN")
    if (!is.data.frame(m) || !nrow(m) || !all(required %in% names(m)) ||
        anyDuplicated(m[keys])) stop("Invalid or duplicate stopping-moment rows")
    if (any(!is.finite(as.matrix(m[c("start_n", "increment", "n_looks")])) ) ||
        any(m$start_n <= 0 | m$increment <= 0 | m$n_looks < 1 | m$n_looks != floor(m$n_looks))) {
        stop("Unsupported stopping schedule")
    }
    tolerance_source <- if (all(c("tolerance_EN", "tolerance_VarN") %in% names(m))) {
        "recorded bundle tolerances"
    } else "legacy v1 observed-maximum reconstruction of the retained MC rule"
    m <- bfpwr_sim_moment_mc_fields(m, bundle$cumulative)
    designs <- bundle$designs
    design_source <- "frozen bundle design specifications"
    if (is.null(designs)) {
        if (!is.null(bundle$timing_schema_version) && bundle$timing_schema_version >= 2L) {
            stop("Frozen bundle$designs are required for refreshed bundles")
        }
        designs <- bfpwr_sim_design_cases_z()
        design_source <- "legacy v1 registry fallback (not frozen in bundle)"
    }
    case_ids <- vapply(bundle$cases, `[[`, "", "bf_prior_id")
    design_ids <- vapply(designs, `[[`, "", "design_case_id")
    if (anyDuplicated(case_ids) || anyDuplicated(design_ids)) stop("Duplicate case or design IDs")
    quadratures <- lapply(c(64L, 128L), bfpwr_sim_moment_legendre)
    for (name in c("coarse_EN", "coarse_VarN", "converged_EN", "converged_VarN", "mass_error")) {
        m[[name]] <- NA_real_
    }
    for (id in unique(m$case_id)) {
        case_index <- match(id, case_ids)
        if (is.na(case_index)) stop("Missing case: ", id)
        case <- bundle$cases[[case_index]]
        design_index <- match(case$design_case_id, design_ids)
        if (is.na(design_index)) stop("Missing design: ", case$design_case_id)
        design <- designs[[design_index]]
        if (is.null(design) || design$test_family != "z" || case$bf_type != "moment" ||
            !case$analysis_prior$alternative %in% c("less", "greater")) {
            stop("Unsupported Markov-reference case: ", id)
        }
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        usd <- design$generation$usd
        prior <- case$analysis_prior
        for (threshold in unique(m$evidence_threshold[m$case_id == id])) {
            index <- which(m$case_id == id & m$evidence_threshold == threshold)
            if (!is.finite(threshold) || threshold <= 1) stop("Unsupported evidence threshold")
            schedules <- lapply(index, function(i) seq(m$start_n[i], by = m$increment[i],
                                                       length.out = m$n_looks[i]))
            grid <- sort(unique(unlist(schedules)))
            boundary <- function(k) vapply(grid, function(ni) zcrit(k, se = usd/sqrt(ni),
                null = prior$null, tau = prior$psd, type = "moment",
                alternative = prior$alternative), numeric(1))
            b1 <- boundary(1/threshold)
            b0 <- boundary(threshold)
            lower <- pmin(b1, b0)
            upper <- pmax(b1, b0)
            ## The longest schedule supplies only exact prefixes; incompatible
            ## schedules get their own recursion.
            pending <- order(lengths(schedules), decreasing = TRUE)
            while (length(pending)) {
                longest <- pending[1L]
                n <- schedules[[longest]]
                j <- match(n, grid)
                values <- lapply(quadratures, function(q) bfpwr_sim_moment_markov(n,
                    lower[j], upper[j], dp$dpm, dp$dpsd, prior$null, usd, q))
                reuse <- pending[vapply(pending, function(h) identical(schedules[[h]],
                    n[seq_along(schedules[[h]])]), logical(1))]
                for (h in reuse) {
                    coarse <- bfpwr_sim_moment_prefix(values[[1L]], n, length(schedules[[h]]))
                    fine <- bfpwr_sim_moment_prefix(values[[2L]], n, length(schedules[[h]]))
                    i <- index[h]
                    m$coarse_EN[i] <- coarse["EN"]
                    m$coarse_VarN[i] <- coarse["VarN"]
                    m$converged_EN[i] <- fine["EN"]
                    m$converged_VarN[i] <- fine["VarN"]
                    m$mass_error[i] <- max(coarse["mass_error"], fine["mass_error"])
                }
                pending <- setdiff(pending, reuse)
            }
        }
    }
    m$refinement_EN <- abs(m$converged_EN - m$coarse_EN)
    m$refinement_VarN <- abs(m$converged_VarN - m$coarse_VarN)
    m$refinement_tolerance_EN <- 1e-8*pmax(1, abs(m$converged_EN))
    m$refinement_tolerance_VarN <- 1e-8*pmax(1, abs(m$converged_VarN))
    m$mass_tolerance <- 1e-8
    m$reference_converged <- is.finite(m$mass_error) & m$mass_error <= m$mass_tolerance &
        m$refinement_EN <= m$refinement_tolerance_EN &
        m$refinement_VarN <= m$refinement_tolerance_VarN
    relative <- function(error, value) ifelse(value != 0, error/value, NA_real_)
    standardized <- function(error, mcse) ifelse(mcse > 0, error/mcse, NA_real_)
    for (quantity in c("EN", "VarN")) {
        reference <- m[[paste0("converged_", quantity)]]
        error <- m[[paste0("reference_", quantity)]] - reference
        mc_error <- m[[quantity]] - reference
        m[[paste0("numerical_error_", quantity)]] <- error
        m[[paste0("numerical_relative_", quantity)]] <- relative(error, reference)
        m[[paste0("numerical_z_", quantity)]] <- standardized(error, m[[paste0("mcse_", quantity)]])
        m[[paste0("mc_reference_error_", quantity)]] <- mc_error
        m[[paste0("mc_reference_relative_", quantity)]] <- relative(mc_error, reference)
        m[[paste0("mc_reference_z_", quantity)]] <- standardized(mc_error, m[[paste0("mcse_", quantity)]])
    }
    list(moment_rows = m, summary = data.frame(
        diagnostic = "Reference convergence check; not a package-default accuracy pass",
        design_source = design_source, tolerance_source = tolerance_source,
        rows_checked = nrow(m), coarse_nodes = 64L,
        fine_nodes = 128L, mass_tolerance = 1e-8, refinement_relative_tolerance = 1e-8,
        failed_rows = sum(!m$reference_converged), passed = all(m$reference_converged),
        max_mass_error = max(m$mass_error), max_refinement_EN = max(m$refinement_EN),
        max_refinement_VarN = max(m$refinement_VarN)))
}
