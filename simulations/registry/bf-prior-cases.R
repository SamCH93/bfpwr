bfpwr_sim_bf_prior_cases <- function(designs = bfpwr_sim_design_cases(),
                                     tiers = NULL,
                                     families = NULL,
                                     tags = NULL) {
    cases <- c(bfpwr_sim_bf_prior_cases_z(designs),
               bfpwr_sim_bf_prior_cases_t(designs),
               bfpwr_sim_bf_prior_cases_binomial(designs))
    bfpwr_sim_registry_filter(cases, tiers = tiers, families = families,
                              tags = tags)
}

bfpwr_sim_find_bf_prior_case <- function(bf_prior_id,
                                         cases = bfpwr_sim_bf_prior_cases()) {
    ids <- vapply(cases, function(x) x$bf_prior_id, character(1))
    hit <- which(ids == bf_prior_id)
    if (length(hit) != 1) {
        stop("BF prior case not found or not unique: ", bf_prior_id)
    }
    cases[[hit]]
}

bfpwr_sim_bf_prior_case_set <- function(name = c("production", "smoke")) {
    name <- match.arg(name)
    switch(name,
           production = bfpwr_sim_bf_prior_cases(),
           smoke = bfpwr_sim_smoke_bf_prior_cases())
}
