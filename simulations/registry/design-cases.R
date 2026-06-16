bfpwr_sim_design_cases <- function(tiers = NULL, families = NULL, tags = NULL) {
    cases <- c(bfpwr_sim_design_cases_z(),
               bfpwr_sim_design_cases_t(),
               bfpwr_sim_design_cases_binomial())
    bfpwr_sim_registry_filter(cases, tiers = tiers, families = families, tags = tags)
}

bfpwr_sim_find_design_case <- function(design_case_id, cases = bfpwr_sim_design_cases()) {
    ids <- vapply(cases, function(x) x$design_case_id, character(1))
    hit <- which(ids == design_case_id)
    if (length(hit) != 1) {
        stop("design case not found or not unique: ", design_case_id)
    }
    cases[[hit]]
}

bfpwr_sim_design_case_set <- function(name = c("production", "smoke")) {
    name <- match.arg(name)
    switch(name,
           production = bfpwr_sim_design_cases(),
           smoke = bfpwr_sim_smoke_design_cases())
}
