bfpwr_sim_bf_prior_cases_t <- function(designs) {
    t_designs <- designs[vapply(designs, function(x) x$test_family == "t",
                                logical(1))]
    ids <- vapply(t_designs, function(x) x$design_case_id, character(1))
    by_id <- function(id) {
        hit <- match(id, ids)
        if (is.na(hit)) stop("t design case not found: ", id)
        t_designs[[hit]]
    }

    cases <- list()
    add <- function(case) {
        cases[[length(cases) + 1L]] <<- case
    }

    for (design in t_designs) {
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 1 / sqrt(2),
            pdf = 1, alternative = "two.sided",
            rationale = "default Cauchy two-sided t-test prior"))
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 1 / sqrt(2),
            pdf = 1, alternative = "greater",
            rationale = "default Cauchy greater one-sided t-test prior"))
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 1 / sqrt(2),
            pdf = 1, alternative = "less",
            rationale = "default Cauchy less one-sided t-test prior"))
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 0.35,
            pdf = 1, alternative = "two.sided",
            rationale = "narrow Cauchy two-sided t-test prior"))
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 1,
            pdf = 1, alternative = "two.sided",
            rationale = "wide Cauchy two-sided t-test prior"))
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 1 / sqrt(2),
            pdf = 3, alternative = "two.sided",
            rationale = "Student-t df=3 two-sided t-test prior"))
        add(bfpwr_sim_bf_prior_t(
            design, null = 0, plocation = 0, pscale = 1 / sqrt(2),
            pdf = 30, alternative = "two.sided",
            rationale = "high-df Student-t two-sided t-test prior"))
    }

    positive_ids <- c(
        "t-two-dpoint-0p5-short",
        "t-two-dnorm-0p5-s0p1-short",
        "t-two-dnorm-0p5-s0p5-short",
        "t-one-dpoint-0p5-short",
        "t-paired-dpoint-0p5-short",
        "t-two-dpoint-0p35-n2x1p1-short",
        "t-two-dpoint-0p2-long",
        "t-two-dnorm-0p2-s0p1-long",
        "t-two-dnorm-0p2-s0p5-long"
    )
    for (id in positive_ids) {
        add(bfpwr_sim_bf_prior_t(
            by_id(id), null = 0, plocation = 0.5, pscale = 0.1,
            pdf = 3, alternative = "greater",
            rationale = "localized positive Student-t prior"))
    }

    local_shift_ids <- c(
        "t-two-dpoint-0p35-n2x1p1-short",
        "t-two-dpoint-0p2-long",
        "t-two-dnorm-0p2-s0p1-long",
        "t-two-dnorm-0p2-s0p5-long"
    )
    for (id in local_shift_ids) {
        add(bfpwr_sim_bf_prior_t(
            by_id(id), null = 0.2, plocation = 0.5, pscale = 0.35,
            pdf = 3, alternative = "greater",
            rationale = "shifted-null local positive Student-t prior"))
    }

    null_shift_ids <- c(
        "t-two-dpoint-0-short",
        "t-one-dpoint-0-short",
        "t-paired-dpoint-0-short",
        "t-two-dpoint-0-n2x1p1-short",
        "t-two-dpoint-0-long"
    )
    for (id in null_shift_ids) {
        add(bfpwr_sim_bf_prior_t(
            by_id(id), null = 0.2, plocation = 0.5, pscale = 0.35,
            pdf = 30, alternative = "greater",
            rationale = "shifted-null high-df positive Student-t prior"))
    }

    cases
}
