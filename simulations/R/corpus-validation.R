bfpwr_sim_validate_design_trajectories <- function(design, trajectories) {
    stopifnot(is.data.frame(trajectories))
    required <- c("design_case_id", "replicate_id", "chunk_id", "look", "n")
    missing <- setdiff(required, names(trajectories))
    if (length(missing) > 0) {
        stop("trajectory is missing columns: ", paste(missing, collapse = ", "))
    }
    if (!all(trajectories$design_case_id == design$design_case_id)) {
        stop("trajectory contains unexpected design_case_id")
    }
    if (!all(is.finite(trajectories$replicate_id)) ||
        !all(trajectories$replicate_id == round(trajectories$replicate_id)) ||
        any(trajectories$replicate_id < 1) ||
        any(trajectories$replicate_id > design$nsim)) {
        stop("trajectory contains replicate_id outside the design range")
    }
    if (!all(is.finite(trajectories$look)) ||
        !all(trajectories$look == round(trajectories$look)) ||
        any(trajectories$look < 1) ||
        any(trajectories$look > length(design$look_grid))) {
        stop("trajectory contains invalid look indexes")
    }
    if (!all(trajectories$n == design$look_grid[trajectories$look])) {
        stop("trajectory n does not match look_grid[look]")
    }
    key <- paste(trajectories$replicate_id, trajectories$look, sep = ":")
    if (anyDuplicated(key)) {
        stop("duplicate trajectory keys found")
    }
    if (!all(sort(unique(trajectories$n)) == design$look_grid)) {
        stop("trajectory grid does not match design look_grid")
    }
    expected_rows <- length(unique(trajectories$replicate_id)) * length(design$look_grid)
    if (nrow(trajectories) != expected_rows) {
        stop("trajectory does not contain every look for every replicate")
    }
    replicate_looks <- table(trajectories$replicate_id)
    if (any(replicate_looks != length(design$look_grid))) {
        stop("one or more replicates do not contain the full look grid")
    }
    chunk_lookup <- design$chunks[match(trajectories$chunk_id, design$chunks$chunk_id), ]
    if (any(is.na(chunk_lookup$chunk_id))) {
        stop("trajectory contains unknown chunk_id")
    }
    if (any(trajectories$replicate_id < chunk_lookup$replicate_start |
            trajectories$replicate_id > chunk_lookup$replicate_end)) {
        stop("replicate_id is outside the range declared for its chunk_id")
    }
    if (design$test_family %in% c("z", "t")) {
        needed <- if (design$test_family == "z") {
            c("se", "true_effect", "estimate")
        } else {
            c("n1", "n2", "true_effect", "estimate", "t")
        }
        missing <- setdiff(needed, names(trajectories))
        if (length(missing) > 0) {
            stop("trajectory is missing columns: ", paste(missing, collapse = ", "))
        }
        if (!all(is.finite(trajectories$true_effect)) ||
            !all(is.finite(trajectories$estimate))) {
            stop("trajectory contains non-finite continuous statistics")
        }
        spread <- stats::aggregate(true_effect ~ replicate_id, trajectories,
                                   function(x) length(unique(x)))
        if (any(spread$true_effect != 1)) {
            stop("true_effect is not constant within replicate")
        }
        if (design$test_family == "z") {
            expected_se <- design$generation$usd / sqrt(trajectories$n)
            if (any(abs(trajectories$se - expected_se) > 1e-12)) {
                stop("z trajectory se does not match usd / sqrt(n)")
            }
        } else {
            if (!all(is.finite(trajectories$t))) {
                stop("t trajectory contains non-finite t statistics")
            }
            if (!all(trajectories$n1 == design$look_grid[trajectories$look])) {
                stop("t trajectory n1 does not match look grid")
            }
            expected_n2 <- if (design$generation$type == "two.sample") {
                pmax(2, round(trajectories$n1 * design$generation$n2_multiplier))
            } else {
                trajectories$n1
            }
            if (!all(trajectories$n2 == expected_n2)) {
                stop("t trajectory n2 does not match generation settings")
            }
        }
    } else {
        missing <- setdiff(c("true_p", "x"), names(trajectories))
        if (length(missing) > 0) {
            stop("trajectory is missing columns: ", paste(missing, collapse = ", "))
        }
        if (!all(is.finite(trajectories$true_p)) ||
            any(trajectories$true_p <= 0) ||
            any(trajectories$true_p >= 1)) {
            stop("binomial trajectory contains invalid true_p")
        }
        if (!all(is.finite(trajectories$x)) ||
            !all(trajectories$x == round(trajectories$x)) ||
            any(trajectories$x < 0) ||
            any(trajectories$x > trajectories$n)) {
            stop("binomial trajectory contains invalid cumulative successes")
        }
        by_rep <- split(trajectories[order(trajectories$replicate_id,
                                           trajectories$look), ],
                        trajectories$replicate_id)
        if (any(vapply(by_rep, function(x) any(diff(x$x) < 0), logical(1)))) {
            stop("binomial cumulative successes decrease within a replicate")
        }
        spread <- stats::aggregate(true_p ~ replicate_id, trajectories,
                                   function(x) length(unique(x)))
        if (any(spread$true_p != 1)) {
            stop("true_p is not constant within replicate")
        }
    }
    invisible(TRUE)
}

bfpwr_sim_validate_summary <- function(summary, mc_abs_tol = 0.002, z_tol = 4) {
    required <- c("analysis_case_id", "design_case_id", "pH1", "pH0", "pInc",
                  "mcse_pH1", "mcse_pH0", "mcse_pInc",
                  "reference_pH1", "reference_pH0", "reference_pInc")
    missing <- setdiff(required, names(summary))
    if (length(missing) > 0) {
        stop("summary is missing columns: ", paste(missing, collapse = ", "))
    }

    check_one <- function(observed, reference, mcse, nsim) {
        if (!is.finite(reference)) return(NA)
        ref_mcse <- if (reference >= 0 && reference <= 1) {
            sqrt(reference * (1 - reference) / nsim)
        } else {
            0
        }
        abs(observed - reference) <= max(mc_abs_tol, z_tol * max(mcse, ref_mcse))
    }
    data.frame(
        analysis_case_id = summary$analysis_case_id,
        has_pH1_reference = is.finite(summary$reference_pH1),
        has_pH0_reference = is.finite(summary$reference_pH0),
        has_pInc_reference = is.finite(summary$reference_pInc),
        pH1_ok = mapply(check_one, summary$pH1, summary$reference_pH1,
                        summary$mcse_pH1, summary$nsim),
        pH0_ok = mapply(check_one, summary$pH0, summary$reference_pH0,
                        summary$mcse_pH0, summary$nsim),
        pInc_ok = mapply(check_one, summary$pInc, summary$reference_pInc,
                         summary$mcse_pInc, summary$nsim)
    )
}

bfpwr_sim_validate_bf_rows <- function(bf_prior, design, bf_rows) {
    stopifnot(is.data.frame(bf_rows))
    required <- c("design_case_id", "replicate_id", "chunk_id", "look", "n",
                  "bf_prior_id", "log_bf01")
    missing <- setdiff(required, names(bf_rows))
    if (length(missing) > 0) {
        stop("BF rows are missing columns: ", paste(missing, collapse = ", "))
    }
    if (!all(bf_rows$bf_prior_id == bf_prior$bf_prior_id)) {
        stop("BF rows contain unexpected bf_prior_id")
    }
    if (!all(bf_rows$design_case_id == design$design_case_id)) {
        stop("BF rows contain unexpected design_case_id")
    }
    if (any(is.na(bf_rows$log_bf01) & !is.nan(bf_rows$log_bf01))) {
        stop("BF rows contain NA log_bf01 values")
    }

    skeleton <- bf_rows[c("design_case_id", "replicate_id", "chunk_id",
                          "look", "n")]
    if (design$test_family == "t") {
        missing_t <- setdiff(c("n1", "n2"), names(bf_rows))
        if (length(missing_t) > 0) {
            stop("t BF rows are missing columns: ",
                 paste(missing_t, collapse = ", "))
        }
        skeleton$n1 <- bf_rows$n1
        skeleton$n2 <- bf_rows$n2
        skeleton$true_effect <- 0
        skeleton$estimate <- 0
        skeleton$t <- 0
    } else if (design$test_family == "z") {
        skeleton$se <- design$generation$usd / sqrt(skeleton$n)
        skeleton$true_effect <- 0
        skeleton$estimate <- 0
    } else {
        skeleton$true_p <- 0.5
        skeleton$x <- 0
    }
    bfpwr_sim_validate_design_trajectories(design, skeleton)
    invisible(TRUE)
}
