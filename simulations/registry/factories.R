bfpwr_sim_decimal_id <- function(x) {
    x <- as.character(x)
    x <- sub("^-", "m", x)
    x <- gsub("[.]", "p", x)
    x <- gsub("[/]", "over", x)
    x
}

bfpwr_sim_threshold_id <- function(k) {
    if (k < 1) {
        paste0("k1over", bfpwr_sim_decimal_id(round(1 / k, 8)))
    } else {
        paste0("k", bfpwr_sim_decimal_id(round(k, 8)))
    }
}

bfpwr_sim_c <- function(...) {
    bfpwr_sim_normalize_tags(unlist(list(...), use.names = FALSE))
}

bfpwr_sim_design_z <- function(id,
                               tier = c("short", "long", "adversarial", "regression"),
                               look_grid = c("short", "long"),
                               nsim = 10000,
                               chunk_size = 1000,
                               seed,
                               dpm,
                               dpsd,
                               usd = sqrt(2),
                               tags = character(),
                               rationale = "") {
    tier <- match.arg(tier, bfpwr_sim_allowed_tiers())
    look_grid <- match.arg(look_grid)
    prior <- if (dpsd == 0) {
        list(family = "point", mean = dpm, sd = 0)
    } else {
        list(family = "normal", mean = dpm, sd = dpsd)
    }
    bfpwr_sim_design_case(
        design_case_id = id,
        test_family = "z",
        look_grid_name = look_grid,
        nsim = nsim,
        chunk_size = chunk_size,
        master_seed = seed,
        design_prior = prior,
        generation = list(usd = usd),
        tier = tier,
        tags = bfpwr_sim_c("z", if (dpsd == 0) "point-design" else "normal-design", tags),
        rationale = rationale
    )
}

bfpwr_sim_design_t <- function(id,
                               tier = c("short", "long", "adversarial", "regression"),
                               look_grid = c("short", "long"),
                               nsim = 10000,
                               chunk_size = 1000,
                               seed,
                               dpm,
                               dpsd,
                               type = c("two.sample", "one.sample", "paired"),
                               n2_multiplier = 1,
                               tags = character(),
                               rationale = "") {
    tier <- match.arg(tier, bfpwr_sim_allowed_tiers())
    look_grid <- match.arg(look_grid)
    type <- match.arg(type)
    prior <- if (dpsd == 0) {
        list(family = "point", mean = dpm, sd = 0)
    } else {
        list(family = "normal", mean = dpm, sd = dpsd)
    }
    bfpwr_sim_design_case(
        design_case_id = id,
        test_family = "t",
        look_grid_name = look_grid,
        nsim = nsim,
        chunk_size = chunk_size,
        master_seed = seed,
        design_prior = prior,
        generation = list(type = type, n2_multiplier = n2_multiplier),
        tier = tier,
        tags = bfpwr_sim_c("t", type,
                           if (dpsd == 0) "point-design" else "normal-design",
                           tags),
        rationale = rationale
    )
}

bfpwr_sim_design_binomial <- function(id,
                                      tier = c("short", "long", "adversarial", "regression"),
                                      look_grid = c("short", "long"),
                                      nsim = 10000,
                                      chunk_size = 1000,
                                      seed,
                                      design_prior,
                                      tags = character(),
                                      rationale = "") {
    tier <- match.arg(tier, bfpwr_sim_allowed_tiers())
    look_grid <- match.arg(look_grid)
    bfpwr_sim_design_case(
        design_case_id = id,
        test_family = "binomial",
        look_grid_name = look_grid,
        nsim = nsim,
        chunk_size = chunk_size,
        master_seed = seed,
        design_prior = design_prior,
        generation = list(),
        tier = tier,
        tags = bfpwr_sim_c("binomial", paste0(design_prior$family, "-design"), tags),
        rationale = rationale
    )
}

bfpwr_sim_analysis_z_normal <- function(design,
                                        id = NULL,
                                        tier = design$tier,
                                        null = 0,
                                        pm,
                                        psd,
                                        k1 = 1 / 10,
                                        k0 = 10,
                                        strict = FALSE,
                                        tags = character(),
                                        rationale = "") {
    if (is.null(id)) {
        id <- paste("z-normal",
                    paste0("pm", bfpwr_sim_decimal_id(pm)),
                    paste0("psd", bfpwr_sim_decimal_id(psd)),
                    bfpwr_sim_threshold_id(k1),
                    bfpwr_sim_threshold_id(k0),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_analysis_case(
        analysis_case_id = id,
        design_case_id = design$design_case_id,
        test_family = "z",
        bf_type = "normal",
        k1 = k1,
        k0 = k0,
        analysis_prior = list(null = null, pm = pm, psd = psd),
        strict = strict,
        tier = tier,
        tags = bfpwr_sim_c("z", "normal-prior", tags),
        rationale = rationale
    )
}

bfpwr_sim_analysis_z_moment <- function(design,
                                        id = NULL,
                                        tier = design$tier,
                                        null = 0,
                                        psd,
                                        k1 = 1 / 10,
                                        k0 = 10,
                                        strict = FALSE,
                                        tags = character(),
                                        rationale = "") {
    if (is.null(id)) {
        id <- paste("z-moment",
                    paste0("psd", bfpwr_sim_decimal_id(psd)),
                    bfpwr_sim_threshold_id(k1),
                    bfpwr_sim_threshold_id(k0),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_analysis_case(
        analysis_case_id = id,
        design_case_id = design$design_case_id,
        test_family = "z",
        bf_type = "moment",
        k1 = k1,
        k0 = k0,
        analysis_prior = list(null = null, psd = psd),
        strict = strict,
        tier = tier,
        tags = bfpwr_sim_c("z", "moment-prior", tags),
        rationale = rationale
    )
}

bfpwr_sim_analysis_z_directional <- function(design,
                                             id = NULL,
                                             tier = design$tier,
                                             null = 0,
                                             pm,
                                             psd,
                                             k1 = 1 / 10,
                                             k0 = 10,
                                             strict = FALSE,
                                             tags = character(),
                                             rationale = "") {
    if (is.null(id)) {
        id <- paste("z-directional",
                    paste0("pm", bfpwr_sim_decimal_id(pm)),
                    paste0("psd", bfpwr_sim_decimal_id(psd)),
                    bfpwr_sim_threshold_id(k1),
                    bfpwr_sim_threshold_id(k0),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_analysis_case(
        analysis_case_id = id,
        design_case_id = design$design_case_id,
        test_family = "z",
        bf_type = "directional",
        k1 = k1,
        k0 = k0,
        analysis_prior = list(null = null, pm = pm, psd = psd),
        strict = strict,
        tier = tier,
        tags = bfpwr_sim_c("z", "directional-prior", tags),
        rationale = rationale
    )
}

bfpwr_sim_analysis_t <- function(design,
                                 id = NULL,
                                 tier = design$tier,
                                 null = 0,
                                 plocation = 0,
                                 pscale = 1 / sqrt(2),
                                 pdf = 1,
                                 alternative = c("two.sided", "less", "greater"),
                                 k1 = 1 / 10,
                                 k0 = 10,
                                 strict = FALSE,
                                 drange = "adaptive",
                                 tags = character(),
                                 rationale = "") {
    alternative <- match.arg(alternative)
    type <- design$generation$type
    if (is.null(id)) {
        id <- paste("t", type, alternative,
                    bfpwr_sim_threshold_id(k1),
                    bfpwr_sim_threshold_id(k0),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_analysis_case(
        analysis_case_id = id,
        design_case_id = design$design_case_id,
        test_family = "t",
        bf_type = "t",
        k1 = k1,
        k0 = k0,
        analysis_prior = list(null = null, plocation = plocation,
                              pscale = pscale, pdf = pdf, type = type,
                              alternative = alternative),
        strict = strict,
        drange = drange,
        tier = tier,
        tags = bfpwr_sim_c("t", type, alternative, tags),
        rationale = rationale
    )
}

bfpwr_sim_analysis_binomial <- function(design,
                                        id = NULL,
                                        tier = design$tier,
                                        bf_type = c("point", "direction"),
                                        p0 = 0.5,
                                        a = 1,
                                        b = 1,
                                        k1 = 1 / 10,
                                        k0 = 10,
                                        tags = character(),
                                        rationale = "") {
    bf_type <- match.arg(bf_type)
    if (is.null(id)) {
        id <- paste("binomial", bf_type,
                    paste0("p0", bfpwr_sim_decimal_id(p0)),
                    bfpwr_sim_threshold_id(k1),
                    bfpwr_sim_threshold_id(k0),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_analysis_case(
        analysis_case_id = id,
        design_case_id = design$design_case_id,
        test_family = "binomial",
        bf_type = bf_type,
        k1 = k1,
        k0 = k0,
        analysis_prior = list(p0 = p0, a = a, b = b),
        strict = FALSE,
        tier = tier,
        tags = bfpwr_sim_c("binomial", paste0(bf_type, "-bf"), tags),
        rationale = rationale
    )
}

bfpwr_sim_bf_prior_z_bf01 <- function(design,
                                      id = NULL,
                                      tier = design$tier,
                                      null = 0,
                                      pm,
                                      psd,
                                      tags = character(),
                                      rationale = "") {
    prior_form <- if (psd == 0) "point" else "normal"
    if (is.null(id)) {
        id <- paste("z-bf01", prior_form,
                    paste0("null", bfpwr_sim_decimal_id(null)),
                    paste0("pm", bfpwr_sim_decimal_id(pm)),
                    paste0("psd", bfpwr_sim_decimal_id(round(psd, 8))),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_bf_prior_case(
        bf_prior_id = id,
        design_case_id = design$design_case_id,
        test_family = "z",
        bf_type = "normal",
        prior_form = prior_form,
        analysis_prior = list(null = null, pm = pm, psd = psd),
        tier = tier,
        tags = bfpwr_sim_c("z", "bf01", paste0(prior_form, "-prior"), tags),
        rationale = rationale
    )
}

bfpwr_sim_bf_prior_z_directional <- function(design,
                                             id = NULL,
                                             tier = design$tier,
                                             null = 0,
                                             pm,
                                             psd,
                                             tags = character(),
                                             rationale = "") {
    if (is.null(id)) {
        id <- paste("z-dirbf01-directional-normal",
                    paste0("null", bfpwr_sim_decimal_id(null)),
                    paste0("pm", bfpwr_sim_decimal_id(pm)),
                    paste0("psd", bfpwr_sim_decimal_id(round(psd, 8))),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_bf_prior_case(
        bf_prior_id = id,
        design_case_id = design$design_case_id,
        test_family = "z",
        bf_type = "directional",
        prior_form = "directional-normal",
        analysis_prior = list(null = null, pm = pm, psd = psd),
        tier = tier,
        tags = bfpwr_sim_c("z", "dirbf01", "directional-normal-prior", tags),
        rationale = rationale
    )
}

bfpwr_sim_bf_prior_z_moment <- function(design,
                                        id = NULL,
                                        tier = design$tier,
                                        null = 0,
                                        psd,
                                        tags = character(),
                                        rationale = "") {
    if (is.null(id)) {
        id <- paste("z-nmbf01-moment",
                    paste0("null", bfpwr_sim_decimal_id(null)),
                    paste0("psd", bfpwr_sim_decimal_id(round(psd, 8))),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_bf_prior_case(
        bf_prior_id = id,
        design_case_id = design$design_case_id,
        test_family = "z",
        bf_type = "moment",
        prior_form = "moment",
        analysis_prior = list(null = null, psd = psd),
        tier = tier,
        tags = bfpwr_sim_c("z", "nmbf01", "moment-prior", tags),
        rationale = rationale
    )
}

bfpwr_sim_bf_prior_t <- function(design,
                                 id = NULL,
                                 tier = design$tier,
                                 null = 0,
                                 plocation = 0,
                                 pscale = 1 / sqrt(2),
                                 pdf = 1,
                                 alternative = c("two.sided", "less", "greater"),
                                 tags = character(),
                                 rationale = "") {
    alternative <- match.arg(alternative)
    prior_form <- if (pdf == 1) "cauchy" else "student-t"
    type <- design$generation$type
    if (is.null(id)) {
        id <- paste("t-tbf01", prior_form,
                    paste0("null", bfpwr_sim_decimal_id(null)),
                    paste0("loc", bfpwr_sim_decimal_id(plocation)),
                    paste0("scale", bfpwr_sim_decimal_id(round(pscale, 8))),
                    paste0("df", bfpwr_sim_decimal_id(pdf)),
                    gsub("[.]", "-", alternative),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_bf_prior_case(
        bf_prior_id = id,
        design_case_id = design$design_case_id,
        test_family = "t",
        bf_type = "t",
        prior_form = prior_form,
        analysis_prior = list(null = null, plocation = plocation,
                              pscale = pscale, pdf = pdf, type = type,
                              alternative = alternative),
        tier = tier,
        tags = bfpwr_sim_c("t", "tbf01", type, prior_form, alternative, tags),
        rationale = rationale
    )
}

bfpwr_sim_bf_prior_binomial <- function(design,
                                        id = NULL,
                                        tier = design$tier,
                                        bf_type = c("point", "direction"),
                                        p0 = 0.5,
                                        a = 1,
                                        b = 1,
                                        tags = character(),
                                        rationale = "") {
    bf_type <- match.arg(bf_type)
    prior_form <- if (bf_type == "point") "beta" else "directional-beta"
    if (is.null(id)) {
        id <- paste("binom-binbf01", bf_type,
                    paste0("p0", bfpwr_sim_decimal_id(p0)),
                    paste0("a", bfpwr_sim_decimal_id(a)),
                    paste0("b", bfpwr_sim_decimal_id(b)),
                    "on", design$design_case_id, sep = "-")
    }
    bfpwr_sim_bf_prior_case(
        bf_prior_id = id,
        design_case_id = design$design_case_id,
        test_family = "binomial",
        bf_type = bf_type,
        prior_form = prior_form,
        analysis_prior = list(p0 = p0, a = a, b = b),
        tier = tier,
        tags = bfpwr_sim_c("binomial", "binbf01",
                           paste0(bf_type, "-bf"), prior_form, tags),
        rationale = rationale
    )
}
