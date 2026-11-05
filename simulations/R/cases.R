bfpwr_sim_make_chunks <- function(nsim, chunk_size, master_seed) {
    stopifnot(
        length(nsim) == 1,
        is.numeric(nsim),
        is.finite(nsim),
        nsim >= 1,
        nsim == round(nsim),
        length(chunk_size) == 1,
        is.numeric(chunk_size),
        is.finite(chunk_size),
        chunk_size >= 1,
        chunk_size == round(chunk_size),
        length(master_seed) == 1,
        is.numeric(master_seed),
        is.finite(master_seed)
    )

    chunk_count <- ceiling(nsim / chunk_size)
    chunk_id <- seq_len(chunk_count)
    replicate_start <- (chunk_id - 1) * chunk_size + 1
    replicate_end <- pmin(chunk_id * chunk_size, nsim)
    seed <- as.integer(master_seed + chunk_id * 1000003L)

    data.frame(
        chunk_id = chunk_id,
        replicate_start = replicate_start,
        replicate_end = replicate_end,
        seed = seed
    )
}

bfpwr_sim_allowed_tiers <- function() {
    c("short", "long", "adversarial", "regression")
}

bfpwr_sim_normalize_tags <- function(tags) {
    if (is.null(tags)) return(character())
    tags <- as.character(tags)
    tags <- tags[nzchar(tags)]
    unique(tags)
}

bfpwr_sim_design_case <- function(design_case_id,
                                  test_family = c("z", "t", "binomial"),
                                  look_grid_name = c("short", "long"),
                                  nsim = 10000,
                                  chunk_size = 1000,
                                  master_seed = 20260524,
                                  design_prior,
                                  generation = list(),
                                  tier = look_grid_name,
                                  tags = character(),
                                  rationale = "",
                                  corpus_version = "v1",
                                  schema_version = 1) {
    test_family <- match.arg(test_family)
    look_grid_name <- match.arg(look_grid_name)
    tier <- match.arg(tier, bfpwr_sim_allowed_tiers())
    look_grid <- bfpwr_sim_grid(look_grid_name)

    case <- list(
        schema_version = schema_version,
        corpus_version = corpus_version,
        design_case_id = design_case_id,
        test_family = test_family,
        tier = tier,
        tags = bfpwr_sim_normalize_tags(tags),
        rationale = rationale,
        look_grid_name = look_grid_name,
        look_grid = look_grid,
        nsim = as.integer(nsim),
        chunk_size = as.integer(chunk_size),
        rng_kind = "L'Ecuyer-CMRG",
        master_seed = as.integer(master_seed),
        chunks = bfpwr_sim_make_chunks(nsim, chunk_size, master_seed),
        design_prior = design_prior,
        generation = bfpwr_sim_complete_generation(test_family, generation)
    )

    bfpwr_sim_validate_design_case(case)
    class(case) <- c("bfpwr_sim_design_case", class(case))
    case
}

bfpwr_sim_analysis_case <- function(analysis_case_id,
                                    design_case_id,
                                    test_family = c("z", "t", "binomial"),
                                    bf_type,
                                    k1,
                                    k0,
                                    analysis_prior,
                                    strict = TRUE,
                                    drange = "adaptive",
                                    tier = "short",
                                    tags = character(),
                                    rationale = "",
                                    corpus_version = "v1",
                                    schema_version = 1) {
    test_family <- match.arg(test_family)
    tier <- match.arg(tier, bfpwr_sim_allowed_tiers())

    case <- list(
        schema_version = schema_version,
        corpus_version = corpus_version,
        analysis_case_id = analysis_case_id,
        design_case_id = design_case_id,
        test_family = test_family,
        tier = tier,
        tags = bfpwr_sim_normalize_tags(tags),
        rationale = rationale,
        bf_type = bf_type,
        k1 = k1,
        k0 = k0,
        analysis_prior = analysis_prior,
        strict = strict,
        drange = drange
    )

    bfpwr_sim_validate_analysis_case(case)
    class(case) <- c("bfpwr_sim_analysis_case", class(case))
    case
}

bfpwr_sim_bf_prior_case <- function(bf_prior_id,
                                    design_case_id,
                                    test_family = c("z", "t", "binomial"),
                                    bf_type,
                                    prior_form,
                                    analysis_prior,
                                    tier = "short",
                                    tags = character(),
                                    rationale = "",
                                    corpus_version = "v1",
                                    schema_version = 1) {
    test_family <- match.arg(test_family)
    tier <- match.arg(tier, bfpwr_sim_allowed_tiers())

    case <- list(
        schema_version = schema_version,
        corpus_version = corpus_version,
        bf_prior_id = bf_prior_id,
        design_case_id = design_case_id,
        test_family = test_family,
        bf_type = bf_type,
        prior_form = prior_form,
        analysis_prior = analysis_prior,
        tier = tier,
        tags = bfpwr_sim_normalize_tags(tags),
        rationale = rationale
    )

    bfpwr_sim_validate_bf_prior_case(case)
    class(case) <- c("bfpwr_sim_bf_prior_case", class(case))
    case
}

bfpwr_sim_complete_generation <- function(test_family, generation) {
    defaults <- switch(test_family,
                       z = list(usd = sqrt(2), type = "generic"),
                       t = list(type = "two.sample", n2_multiplier = 1),
                       binomial = list())
    utils::modifyList(defaults, generation)
}

bfpwr_sim_design_prior_mean_sd <- function(design) {
    prior <- design$design_prior
    if (!design$test_family %in% c("z", "t")) {
        stop("continuous design prior requested for non-continuous design")
    }
    if (identical(prior$family, "point")) {
        return(list(dpm = prior$mean, dpsd = 0))
    }
    if (identical(prior$family, "normal")) {
        return(list(dpm = prior$mean, dpsd = prior$sd))
    }
    stop("unsupported continuous design prior: ", prior$family)
}

bfpwr_sim_design_case_override <- function(design,
                                           nsim = design$nsim,
                                           chunk_size = design$chunk_size,
                                           master_seed = design$master_seed,
                                           look_grid_name = design$look_grid_name) {
    bfpwr_sim_design_case(
        design_case_id = design$design_case_id,
        test_family = design$test_family,
        look_grid_name = look_grid_name,
        nsim = nsim,
        chunk_size = chunk_size,
        master_seed = master_seed,
        design_prior = design$design_prior,
        generation = design$generation,
        tier = design$tier,
        tags = design$tags,
        rationale = design$rationale,
        corpus_version = design$corpus_version,
        schema_version = design$schema_version
    )
}

bfpwr_sim_binomial_design_args <- function(design) {
    prior <- design$design_prior
    if (identical(prior$family, "point")) {
        return(list(dp = prior$prob, da = NA_real_, db = NA_real_,
                    dl = NA_real_, du = NA_real_))
    }
    if (identical(prior$family, "beta")) {
        return(list(dp = NA_real_, da = prior$shape1, db = prior$shape2,
                    dl = prior$lower, du = prior$upper))
    }
    stop("unsupported binomial design prior: ", prior$family)
}

bfpwr_sim_validate_design_case <- function(case) {
    stopifnot(
        is.list(case),
        length(case$design_case_id) == 1,
        is.character(case$design_case_id),
        nzchar(case$design_case_id),
        case$test_family %in% c("z", "t", "binomial"),
        case$tier %in% bfpwr_sim_allowed_tiers(),
        is.character(case$tags),
        length(case$rationale) == 1,
        is.character(case$rationale),
        case$look_grid_name %in% c("short", "long"),
        case$nsim >= 1,
        case$chunk_size >= 1,
        is.data.frame(case$chunks)
    )
    bfpwr_sim_validate_look_grid(case$look_grid)

    prior <- case$design_prior
    stopifnot(is.list(prior), length(prior$family) == 1)
    if (case$test_family %in% c("z", "t")) {
        if (identical(prior$family, "point")) {
            stopifnot(is.numeric(prior$mean), length(prior$mean) == 1)
        } else if (identical(prior$family, "normal")) {
            stopifnot(is.numeric(prior$mean), length(prior$mean) == 1,
                      is.numeric(prior$sd), length(prior$sd) == 1,
                      prior$sd > 0)
        } else {
            stop("unsupported continuous design prior: ", prior$family)
        }
    } else {
        if (identical(prior$family, "point")) {
            stopifnot(is.numeric(prior$prob), prior$prob > 0, prior$prob < 1)
        } else if (identical(prior$family, "beta")) {
            stopifnot(is.numeric(prior$shape1), prior$shape1 > 0,
                      is.numeric(prior$shape2), prior$shape2 > 0,
                      is.numeric(prior$lower), is.numeric(prior$upper),
                      prior$lower >= 0, prior$lower < prior$upper,
                      prior$upper <= 1)
        } else {
            stop("unsupported binomial design prior: ", prior$family)
        }
    }
    invisible(case)
}

bfpwr_sim_validate_bf_prior_case <- function(case) {
    stopifnot(
        is.list(case),
        length(case$bf_prior_id) == 1,
        is.character(case$bf_prior_id),
        nzchar(case$bf_prior_id),
        length(case$design_case_id) == 1,
        is.character(case$design_case_id),
        nzchar(case$design_case_id),
        case$test_family %in% c("z", "t", "binomial"),
        case$tier %in% bfpwr_sim_allowed_tiers(),
        length(case$bf_type) == 1,
        is.character(case$bf_type),
        length(case$prior_form) == 1,
        is.character(case$prior_form),
        is.list(case$analysis_prior),
        is.character(case$tags),
        length(case$rationale) == 1,
        is.character(case$rationale)
    )

    prior <- case$analysis_prior
    if (case$test_family == "z") {
        stopifnot(case$bf_type %in% c("normal", "directional", "moment"))
        if (case$bf_type == "normal") {
            stopifnot(
                is.numeric(prior$null), length(prior$null) == 1,
                is.numeric(prior$pm), length(prior$pm) == 1,
                is.numeric(prior$psd), length(prior$psd) == 1,
                prior$psd >= 0,
                case$prior_form %in% c("point", "normal")
            )
            if (prior$psd == 0 && case$prior_form != "point") {
                stop("bf01 psd = 0 priors must use prior_form = 'point'")
            }
            if (prior$psd > 0 && case$prior_form != "normal") {
                stop("bf01 psd > 0 priors must use prior_form = 'normal'")
            }
        } else if (case$bf_type == "directional") {
            stopifnot(
                is.numeric(prior$null), length(prior$null) == 1,
                is.numeric(prior$pm), length(prior$pm) == 1,
                is.numeric(prior$psd), length(prior$psd) == 1,
                prior$psd > 0,
                identical(case$prior_form, "directional-normal")
            )
        } else {
            stopifnot(
                is.numeric(prior$null), length(prior$null) == 1,
                is.numeric(prior$psd), length(prior$psd) == 1,
                prior$psd > 0,
                identical(case$prior_form, "moment")
            )
        }
    } else if (case$test_family == "t") {
        stopifnot(
            identical(case$bf_type, "t"),
            case$prior_form %in% c("cauchy", "student-t"),
            is.numeric(prior$null), length(prior$null) == 1,
            is.numeric(prior$plocation), length(prior$plocation) == 1,
            is.numeric(prior$pscale), length(prior$pscale) == 1,
            prior$pscale > 0,
            is.numeric(prior$pdf), length(prior$pdf) == 1,
            prior$pdf > 0,
            prior$type %in% c("two.sample", "one.sample", "paired"),
            prior$alternative %in% c("two.sided", "less", "greater")
        )
        if (prior$pdf == 1 && case$prior_form != "cauchy") {
            stop("t priors with pdf = 1 must use prior_form = 'cauchy'")
        }
        if (prior$pdf != 1 && case$prior_form != "student-t") {
            stop("t priors with pdf != 1 must use prior_form = 'student-t'")
        }
    } else {
        stopifnot(
            case$bf_type %in% c("point", "direction"),
            case$prior_form %in% c("beta", "directional-beta"),
            is.numeric(prior$p0), length(prior$p0) == 1,
            prior$p0 > 0, prior$p0 < 1,
            is.numeric(prior$a), length(prior$a) == 1,
            prior$a > 0,
            is.numeric(prior$b), length(prior$b) == 1,
            prior$b > 0
        )
        if (case$bf_type == "point" && case$prior_form != "beta") {
            stop("binomial point BF priors must use prior_form = 'beta'")
        }
        if (case$bf_type == "direction" &&
            case$prior_form != "directional-beta") {
            stop("binomial directional BF priors must use prior_form = 'directional-beta'")
        }
    }
    invisible(case)
}

bfpwr_sim_validate_analysis_case <- function(case) {
    stopifnot(
        is.list(case),
        length(case$analysis_case_id) == 1,
        is.character(case$analysis_case_id),
        nzchar(case$analysis_case_id),
        length(case$design_case_id) == 1,
        is.character(case$design_case_id),
        nzchar(case$design_case_id),
        case$test_family %in% c("z", "t", "binomial"),
        case$tier %in% bfpwr_sim_allowed_tiers(),
        is.character(case$tags),
        length(case$rationale) == 1,
        is.character(case$rationale),
        length(case$bf_type) == 1,
        is.character(case$bf_type),
        is.numeric(case$k1),
        length(case$k1) == 1,
        case$k1 > 0,
        case$k1 <= 1,
        is.numeric(case$k0),
        length(case$k0) == 1,
        case$k0 >= 1,
        is.list(case$analysis_prior),
        is.logical(case$strict),
        length(case$strict) == 1
    )
    if (case$test_family == "z") {
        stopifnot(case$bf_type %in% c("normal", "directional", "moment"))
        prior <- case$analysis_prior
        if (case$bf_type == "normal") {
            stopifnot(
                is.numeric(prior$null), length(prior$null) == 1,
                is.numeric(prior$pm), length(prior$pm) == 1,
                is.numeric(prior$psd), length(prior$psd) == 1,
                prior$psd >= 0
            )
        } else if (case$bf_type == "directional") {
            stopifnot(
                is.numeric(prior$null), length(prior$null) == 1,
                is.numeric(prior$pm), length(prior$pm) == 1,
                is.numeric(prior$psd), length(prior$psd) == 1,
                prior$psd > 0
            )
        } else {
            stopifnot(
                is.numeric(prior$null), length(prior$null) == 1,
                is.numeric(prior$psd), length(prior$psd) == 1,
                prior$psd > 0
            )
        }
    } else if (case$test_family == "t") {
        stopifnot(identical(case$bf_type, "t"))
        prior <- case$analysis_prior
        stopifnot(
            is.numeric(prior$null), length(prior$null) == 1,
            is.numeric(prior$plocation), length(prior$plocation) == 1,
            is.numeric(prior$pscale), length(prior$pscale) == 1,
            prior$pscale > 0,
            is.numeric(prior$pdf), length(prior$pdf) == 1,
            prior$pdf > 0,
            prior$type %in% c("two.sample", "one.sample", "paired"),
            prior$alternative %in% c("two.sided", "less", "greater")
        )
    } else {
        stopifnot(case$bf_type %in% c("point", "direction"))
        prior <- case$analysis_prior
        stopifnot(
            is.numeric(prior$p0), length(prior$p0) == 1,
            prior$p0 > 0, prior$p0 < 1,
            is.numeric(prior$a), length(prior$a) == 1,
            prior$a > 0,
            is.numeric(prior$b), length(prior$b) == 1,
            prior$b > 0
        )
    }
    invisible(case)
}
