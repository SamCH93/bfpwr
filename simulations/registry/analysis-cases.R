bfpwr_sim_analysis_cases <- function(tiers = NULL, families = NULL, tags = NULL,
                                     designs = bfpwr_sim_design_cases()) {
    cases <- c(bfpwr_sim_analysis_cases_z(designs),
               bfpwr_sim_analysis_cases_t(designs),
               bfpwr_sim_analysis_cases_binomial(designs))
    bfpwr_sim_registry_filter(cases, tiers = tiers, families = families, tags = tags)
}

bfpwr_sim_find_analysis_case <- function(analysis_case_id,
                                         cases = bfpwr_sim_analysis_cases()) {
    ids <- vapply(cases, function(x) x$analysis_case_id, character(1))
    hit <- which(ids == analysis_case_id)
    if (length(hit) != 1) {
        stop("analysis case not found or not unique: ", analysis_case_id)
    }
    cases[[hit]]
}

bfpwr_sim_analysis_case_set <- function(name = c("production", "smoke")) {
    name <- match.arg(name)
    switch(name,
           production = bfpwr_sim_analysis_cases(),
           smoke = bfpwr_sim_smoke_analysis_cases())
}
