bfpwr_sim_registry_filter <- function(cases,
                                      tiers = NULL,
                                      families = NULL,
                                      tags = NULL) {
    if (is.null(cases) || length(cases) == 0) return(list())
    keep <- rep(TRUE, length(cases))
    if (!is.null(tiers)) {
        tiers <- match.arg(tiers, bfpwr_sim_allowed_tiers(), several.ok = TRUE)
        keep <- keep & vapply(cases, function(x) x$tier %in% tiers, logical(1))
    }
    if (!is.null(families)) {
        families <- match.arg(families, c("z", "t", "binomial"), several.ok = TRUE)
        keep <- keep & vapply(cases, function(x) x$test_family %in% families, logical(1))
    }
    if (!is.null(tags)) {
        tags <- as.character(tags)
        keep <- keep & vapply(cases, function(x) all(tags %in% x$tags), logical(1))
    }
    cases[keep]
}

bfpwr_sim_case_ids <- function(cases, id_name) {
    vapply(cases, function(x) x[[id_name]], character(1))
}

bfpwr_sim_assert_unique_ids <- function(cases, id_name) {
    ids <- bfpwr_sim_case_ids(cases, id_name)
    dup <- unique(ids[duplicated(ids)])
    if (length(dup) > 0) {
        stop("duplicate case ids: ", paste(dup, collapse = ", "))
    }
    invisible(TRUE)
}

bfpwr_sim_assert_analysis_designs_exist <- function(analyses, designs) {
    design_ids <- bfpwr_sim_case_ids(designs, "design_case_id")
    missing <- setdiff(bfpwr_sim_case_ids(analyses, "design_case_id"), design_ids)
    if (length(missing) > 0) {
        stop("analysis cases reference unknown design ids: ",
             paste(missing, collapse = ", "))
    }
    invisible(TRUE)
}

bfpwr_sim_assert_bf_prior_designs_exist <- function(bf_priors, designs) {
    design_ids <- bfpwr_sim_case_ids(designs, "design_case_id")
    missing <- setdiff(bfpwr_sim_case_ids(bf_priors, "design_case_id"),
                       design_ids)
    if (length(missing) > 0) {
        stop("BF prior cases reference unknown design ids: ",
             paste(missing, collapse = ", "))
    }
    invisible(TRUE)
}

bfpwr_sim_assert_tiers_match_grids <- function(designs) {
    bad_short <- vapply(designs, function(x) {
        identical(x$tier, "short") && !identical(x$look_grid_name, "short")
    }, logical(1))
    bad_long <- vapply(designs, function(x) {
        identical(x$tier, "long") && !identical(x$look_grid_name, "long")
    }, logical(1))
    if (any(bad_short | bad_long)) {
        bad <- bfpwr_sim_case_ids(designs[bad_short | bad_long], "design_case_id")
        stop("regular short/long tiers must use matching look grids: ",
             paste(bad, collapse = ", "))
    }
    invisible(TRUE)
}

bfpwr_sim_validate_registry <- function(designs, analyses,
                                        bf_priors = list()) {
    bfpwr_sim_assert_unique_ids(designs, "design_case_id")
    bfpwr_sim_assert_unique_ids(analyses, "analysis_case_id")
    bfpwr_sim_assert_unique_ids(bf_priors, "bf_prior_id")
    bfpwr_sim_assert_analysis_designs_exist(analyses, designs)
    bfpwr_sim_assert_bf_prior_designs_exist(bf_priors, designs)
    bfpwr_sim_assert_tiers_match_grids(designs)
    invisible(TRUE)
}

bfpwr_sim_flatten_design_case <- function(case) {
    prior <- case$design_prior
    data.frame(
        design_case_id = case$design_case_id,
        family = case$test_family,
        tier = case$tier,
        look_grid = case$look_grid_name,
        nsim = case$nsim,
        chunk_size = case$chunk_size,
        n_chunks = nrow(case$chunks),
        estimated_rows = case$nsim * length(case$look_grid),
        seed = case$master_seed,
        design_prior_family = prior$family,
        dpm = if (!is.null(prior$mean)) prior$mean else NA_real_,
        dpsd = if (!is.null(prior$sd)) prior$sd else NA_real_,
        dp = if (!is.null(prior$prob)) prior$prob else NA_real_,
        da = if (!is.null(prior$shape1)) prior$shape1 else NA_real_,
        db = if (!is.null(prior$shape2)) prior$shape2 else NA_real_,
        dl = if (!is.null(prior$lower)) prior$lower else NA_real_,
        du = if (!is.null(prior$upper)) prior$upper else NA_real_,
        generation = paste(names(case$generation), unlist(case$generation),
                           sep = "=", collapse = ";"),
        tags = paste(case$tags, collapse = ";"),
        rationale = case$rationale,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_flatten_analysis_case <- function(case) {
    prior <- case$analysis_prior
    data.frame(
        analysis_case_id = case$analysis_case_id,
        design_case_id = case$design_case_id,
        family = case$test_family,
        tier = case$tier,
        bf_type = case$bf_type,
        k1 = case$k1,
        k0 = case$k0,
        analysis_prior = paste(names(prior), unlist(prior),
                               sep = "=", collapse = ";"),
        strict = case$strict,
        drange = paste(case$drange, collapse = ";"),
        tags = paste(case$tags, collapse = ";"),
        rationale = case$rationale,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_flatten_bf_prior_case <- function(case) {
    prior <- case$analysis_prior
    data.frame(
        bf_prior_id = case$bf_prior_id,
        design_case_id = case$design_case_id,
        family = case$test_family,
        tier = case$tier,
        bf_type = case$bf_type,
        prior_form = case$prior_form,
        analysis_prior = paste(names(prior), unlist(prior),
                               sep = "=", collapse = ";"),
        tags = paste(case$tags, collapse = ";"),
        rationale = case$rationale,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_design_manifest <- function(designs) {
    if (length(designs) == 0) {
        return(data.frame(
            design_case_id = character(),
            family = character(),
            tier = character(),
            look_grid = character(),
            nsim = integer(),
            chunk_size = integer(),
            n_chunks = integer(),
            estimated_rows = integer(),
            seed = integer(),
            design_prior_family = character(),
            dpm = numeric(),
            dpsd = numeric(),
            dp = numeric(),
            da = numeric(),
            db = numeric(),
            dl = numeric(),
            du = numeric(),
            generation = character(),
            tags = character(),
            rationale = character()
        ))
    }
    do.call(rbind, lapply(designs, bfpwr_sim_flatten_design_case))
}

bfpwr_sim_analysis_manifest <- function(analyses) {
    if (length(analyses) == 0) {
        return(data.frame(
            analysis_case_id = character(),
            design_case_id = character(),
            family = character(),
            tier = character(),
            bf_type = character(),
            k1 = numeric(),
            k0 = numeric(),
            analysis_prior = character(),
            strict = logical(),
            drange = character(),
            tags = character(),
            rationale = character()
        ))
    }
    do.call(rbind, lapply(analyses, bfpwr_sim_flatten_analysis_case))
}

bfpwr_sim_bf_prior_manifest <- function(bf_priors) {
    if (length(bf_priors) == 0) {
        return(data.frame(
            bf_prior_id = character(),
            design_case_id = character(),
            family = character(),
            tier = character(),
            bf_type = character(),
            prior_form = character(),
            analysis_prior = character(),
            tags = character(),
            rationale = character()
        ))
    }
    do.call(rbind, lapply(bf_priors, bfpwr_sim_flatten_bf_prior_case))
}

bfpwr_sim_run_plan <- function(designs, analyses) {
    if (length(designs) == 0) {
        return(data.frame(
            phase = character(),
            tier = character(),
            family = character(),
            design_case_id = character(),
            analysis_case_id = character(),
            chunk_id = integer(),
            replicate_start = integer(),
            replicate_end = integer(),
            seed = integer()
        ))
    }
    generation <- do.call(rbind, lapply(designs, function(design) {
        data.frame(
            phase = "generate_design",
            tier = design$tier,
            family = design$test_family,
            design_case_id = design$design_case_id,
            analysis_case_id = NA_character_,
            chunk_id = design$chunks$chunk_id,
            replicate_start = design$chunks$replicate_start,
            replicate_end = design$chunks$replicate_end,
            seed = design$chunks$seed,
            stringsAsFactors = FALSE
        )
    }))
    if (length(analyses) == 0) {
        return(generation)
    }
    materialization <- do.call(rbind, lapply(analyses, function(analysis) {
        data.frame(
            phase = "materialize_analysis",
            tier = analysis$tier,
            family = analysis$test_family,
            design_case_id = analysis$design_case_id,
            analysis_case_id = analysis$analysis_case_id,
            chunk_id = NA_integer_,
            replicate_start = NA_integer_,
            replicate_end = NA_integer_,
            seed = NA_integer_,
            stringsAsFactors = FALSE
        )
    }))
    rbind(generation, materialization)
}

bfpwr_sim_bf_prior_run_plan <- function(bf_priors, designs) {
    if (length(bf_priors) == 0) {
        return(data.frame(
            phase = character(),
            tier = character(),
            family = character(),
            design_case_id = character(),
            bf_prior_id = character(),
            chunk_id = integer(),
            replicate_start = integer(),
            replicate_end = integer(),
            seed = integer()
        ))
    }
    ids <- vapply(designs, function(x) x$design_case_id, character(1))
    do.call(rbind, lapply(bf_priors, function(bf_prior) {
        design <- designs[[match(bf_prior$design_case_id, ids)]]
        data.frame(
            phase = "materialize_bf_prior",
            tier = bf_prior$tier,
            family = bf_prior$test_family,
            design_case_id = bf_prior$design_case_id,
            bf_prior_id = bf_prior$bf_prior_id,
            chunk_id = design$chunks$chunk_id,
            replicate_start = design$chunks$replicate_start,
            replicate_end = design$chunks$replicate_end,
            seed = design$chunks$seed,
            stringsAsFactors = FALSE
        )
    }))
}

bfpwr_sim_write_registry_manifests <- function(output_dir, designs, analyses,
                                               bf_priors = list()) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    bfpwr_sim_validate_registry(designs, analyses, bf_priors)
    utils::write.csv(bfpwr_sim_design_manifest(designs),
                     file.path(output_dir, "design-cases.csv"),
                     row.names = FALSE)
    utils::write.csv(bfpwr_sim_analysis_manifest(analyses),
                     file.path(output_dir, "analysis-cases.csv"),
                     row.names = FALSE)
    utils::write.csv(bfpwr_sim_bf_prior_manifest(bf_priors),
                     file.path(output_dir, "bf-prior-cases.csv"),
                     row.names = FALSE)
    utils::write.csv(bfpwr_sim_run_plan(designs, analyses),
                     file.path(output_dir, "run-plan.csv"),
                     row.names = FALSE)
    utils::write.csv(bfpwr_sim_bf_prior_run_plan(bf_priors, designs),
                     file.path(output_dir, "bf-prior-run-plan.csv"),
                     row.names = FALSE)
    invisible(output_dir)
}
