bfpwr_sim_package_function <- function(name) {
    fn <- get(name, mode = "function", inherits = TRUE)
    if (!is.function(fn)) {
        stop("required bfpwr function is not available: ", name)
    }
    fn
}

bfpwr_sim_validate_bf_evaluation_case <- function(case) {
    if (!is.null(case$bf_prior_id)) {
        return(bfpwr_sim_validate_bf_prior_case(case))
    }
    bfpwr_sim_validate_analysis_case(case)
}

bfpwr_sim_bf_evaluation_id <- function(case) {
    if (!is.null(case$bf_prior_id)) {
        return(case$bf_prior_id)
    }
    case$analysis_case_id
}

bfpwr_sim_log_bf <- function(analysis, trajectories) {
    bfpwr_sim_validate_bf_evaluation_case(analysis)
    prior <- analysis$analysis_prior

    if (analysis$test_family == "z") {
        if (analysis$bf_type == "normal") {
            return(bfpwr_sim_package_function("bf01")(
                estimate = trajectories$estimate,
                se = trajectories$se,
                null = prior$null,
                pm = prior$pm,
                psd = prior$psd,
                log = TRUE
            ))
        }
        if (analysis$bf_type == "directional") {
            return(bfpwr_sim_package_function("dirbf01")(
                estimate = trajectories$estimate,
                se = trajectories$se,
                null = prior$null,
                pm = prior$pm,
                psd = prior$psd,
                log = TRUE
            ))
        }
        if (analysis$bf_type == "moment") {
            return(bfpwr_sim_package_function("nmbf01")(
                estimate = trajectories$estimate,
                se = trajectories$se,
                null = prior$null,
                psd = prior$psd,
                alternative = if (is.null(prior$alternative)) "two.sided" else prior$alternative,
                log = TRUE
            ))
        }
    }

    if (analysis$test_family == "t") {
        neff <- if (prior$type == "two.sample") {
            1 / (1 / trajectories$n1 + 1 / trajectories$n2)
        } else {
            trajectories$n1
        }
        t_eval <- (trajectories$estimate - prior$null) * sqrt(neff)
        return(bfpwr_sim_package_function("tbf01")(
            t = t_eval,
            n = trajectories$n1,
            n1 = trajectories$n1,
            n2 = trajectories$n2,
            plocation = prior$plocation - prior$null,
            pscale = prior$pscale,
            pdf = prior$pdf,
            type = prior$type,
            alternative = prior$alternative,
            log = TRUE
        ))
    }

    if (analysis$test_family == "binomial") {
        return(bfpwr_sim_package_function("binbf01")(
            x = trajectories$x,
            n = trajectories$n,
            p0 = prior$p0,
            type = analysis$bf_type,
            a = prior$a,
            b = prior$b,
            log = TRUE
        ))
    }

    stop("unsupported analysis family/type")
}

bfpwr_sim_bf_rows <- function(analysis, trajectories) {
    log_bf01 <- bfpwr_sim_log_bf(analysis, trajectories)
    bf_rows <- trajectories[c("design_case_id", "replicate_id", "chunk_id",
                              "look", "n")]
    if ("n1" %in% names(trajectories)) bf_rows$n1 <- trajectories$n1
    if ("n2" %in% names(trajectories)) bf_rows$n2 <- trajectories$n2
    bf_rows$bf_prior_id <- bfpwr_sim_bf_evaluation_id(analysis)
    if (!is.null(analysis$analysis_case_id)) {
        bf_rows$analysis_case_id <- analysis$analysis_case_id
    }
    bf_rows$log_bf01 <- log_bf01
    bf_rows
}
