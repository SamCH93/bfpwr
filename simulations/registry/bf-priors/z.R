bfpwr_sim_bf_prior_cases_z <- function(designs) {
    z_designs <- designs[vapply(designs, function(x) x$test_family == "z",
                                logical(1))]
    ids <- vapply(z_designs, function(x) x$design_case_id, character(1))
    by_id <- function(id) {
        hit <- match(id, ids)
        if (is.na(hit)) stop("z design case not found: ", id)
        z_designs[[hit]]
    }

    cases <- list()
    add <- function(case) {
        cases[[length(cases) + 1L]] <<- case
    }

    standard <- z_designs[vapply(z_designs, function(x) {
        identical(x$generation$usd, 1)
    }, logical(1))]
    for (design in standard) {
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 0.2, psd = 0,
            rationale = "local point alternative"))
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 0.5, psd = 0,
            rationale = "moderate point alternative"))
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 0, psd = 1 / sqrt(2),
            rationale = "centered default-scale normal alternative"))
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 0.2, psd = 0.3,
            rationale = "uncertain local normal alternative"))
        add(bfpwr_sim_bf_prior_z_directional(
            design, null = 0, pm = 0, psd = 1 / sqrt(2),
            rationale = "directional normal alternative"))
        add(bfpwr_sim_bf_prior_z_moment(
            design, null = 0, psd = 0.5 / sqrt(2),
            rationale = "narrower moment alternative"))
        add(bfpwr_sim_bf_prior_z_moment(
            design, null = 0, psd = 1 / sqrt(2),
            rationale = "broader moment alternative"))
    }

    extreme <- z_designs[vapply(z_designs, function(x) {
        identical(x$generation$usd, 50)
    }, logical(1))]
    for (design in extreme) {
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 10, psd = 0,
            rationale = "matched extreme-scale point alternative"))
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 10, psd = 20,
            rationale = "matched uncertain extreme-scale normal alternative"))
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0, pm = 0, psd = 20,
            rationale = "diffuse extreme-scale normal alternative"))
        add(bfpwr_sim_bf_prior_z_directional(
            design, null = 0, pm = 0, psd = 20,
            rationale = "extreme-scale directional normal alternative"))
        add(bfpwr_sim_bf_prior_z_moment(
            design, null = 0, psd = 10 / sqrt(2),
            rationale = "extreme-scale moment alternative"))
    }

    shifted_ids <- c(
        "z-dpoint-0p2-usd1-short",
        "z-dpoint-0p5-usd1-short",
        "z-dnorm-0p5-s0p1-usd1-short",
        "z-dnorm-0p2-s0p3-usd1-short",
        "z-dnorm-0p2-s0p5-usd1-short",
        "z-dpoint-0p2-usd1-long",
        "z-dnorm-0p2-s0p1-usd1-long",
        "z-dnorm-0p2-s0p5-usd1-long"
    )
    for (id in shifted_ids) {
        design <- by_id(id)
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0.2, pm = 0.5, psd = 0,
            rationale = "shifted-null point alternative"))
        add(bfpwr_sim_bf_prior_z_bf01(
            design, null = 0.2, pm = 0.2, psd = 1 / sqrt(2),
            rationale = "shifted-null centered normal alternative"))
        add(bfpwr_sim_bf_prior_z_directional(
            design, null = 0.2, pm = 0.2, psd = 1 / sqrt(2),
            rationale = "shifted-null directional normal alternative"))
        add(bfpwr_sim_bf_prior_z_moment(
            design, null = 0.2, psd = 0.5 / sqrt(2),
            rationale = "shifted-null moment alternative"))
    }

    cases
}
