bfpwr_sim_bf_prior_cases_binomial <- function(designs) {
    binom_designs <- designs[vapply(designs, function(x) {
        x$test_family == "binomial"
    }, logical(1))]
    ids <- vapply(binom_designs, function(x) x$design_case_id, character(1))
    by_id <- function(id) {
        hit <- match(id, ids)
        if (is.na(hit)) stop("binomial design case not found: ", id)
        binom_designs[[hit]]
    }

    cases <- list()
    add <- function(case) {
        cases[[length(cases) + 1L]] <<- case
    }

    for (design in binom_designs) {
        add(bfpwr_sim_bf_prior_binomial(
            design, bf_type = "point", p0 = 0.5, a = 1, b = 1,
            rationale = "default point-null beta(1,1) alternative"))
        add(bfpwr_sim_bf_prior_binomial(
            design, bf_type = "direction", p0 = 0.5, a = 1, b = 1,
            rationale = "default directional beta(1,1) alternative"))
        add(bfpwr_sim_bf_prior_binomial(
            design, bf_type = "point", p0 = 0.5, a = 60, b = 40,
            rationale = "concentrated beta alternative around p=0.6"))
        add(bfpwr_sim_bf_prior_binomial(
            design, bf_type = "point", p0 = 0.5, a = 40, b = 60,
            rationale = "concentrated opposite beta alternative around p=0.4"))
    }

    for (type in c("point", "direction")) {
        add(bfpwr_sim_bf_prior_binomial(
            by_id("binom-dpoint-0p2-short"), bf_type = type,
            p0 = 0.2, a = 1, b = 1,
            rationale = "low-proportion null check"))
        add(bfpwr_sim_bf_prior_binomial(
            by_id("binom-dpoint-0p75-short"), bf_type = type,
            p0 = 0.75, a = 1, b = 1,
            rationale = "high-proportion null check"))
        add(bfpwr_sim_bf_prior_binomial(
            by_id("binom-dpoint-0p6-short"), bf_type = type,
            p0 = 0.6, a = 1, b = 1,
            rationale = "p=0.6 null check"))
        add(bfpwr_sim_bf_prior_binomial(
            by_id("binom-dbeta-60-40-short"), bf_type = type,
            p0 = 0.6, a = 1, b = 1,
            rationale = "p=0.6 null check against beta design prior"))
        add(bfpwr_sim_bf_prior_binomial(
            by_id("binom-dpoint-0p55-long"), bf_type = type,
            p0 = 0.55, a = 1, b = 1,
            rationale = "long-grid local null check"))
    }

    for (id in c("binom-dpoint-0p5-long", "binom-dpoint-0p55-long")) {
        add(bfpwr_sim_bf_prior_binomial(
            by_id(id), bf_type = "point", p0 = 0.5,
            a = 5100, b = 4900,
            rationale = "large-shape beta numeric stress check"))
    }

    cases
}
