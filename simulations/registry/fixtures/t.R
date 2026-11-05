bfpwr_sim_t_fixture_schedules_core <- function() {
    bfpwr_sim_z_fixture_schedules_core()
}

bfpwr_sim_fixture_specs_t <- function() {
    evidence_thresholds <- bfpwr_sim_fixture_evidence_thresholds_core()
    search_targets <- bfpwr_sim_fixture_search_targets_core()
    sequential_search_targets <- bfpwr_sim_search_targets(
        evidence_thresholds = evidence_thresholds,
        target_prob = numeric(0),
        evidence = c("H1", "H0"))
    schedules <- bfpwr_sim_t_fixture_schedules_core()

    list(
        bfpwr_sim_fixed_fixture_spec(
            fixture_set_id = "t-tbf01-fixed-core-v1",
            family = "t",
            bf_types = "t",
            evidence_thresholds = evidence_thresholds,
            search_targets = search_targets,
            tags = c("t", "tbf01", "fixed"),
            rationale = paste(
                "Fixed-look t-test BF01 fixture summaries for threshold,",
                "sample-size search, sampling-type, prior-family, and",
                "one-sided analysis checks.")),
        bfpwr_sim_sequential_fixture_spec(
            fixture_set_id = "t-tbf01-sequential-core-v1",
            family = "t",
            bf_types = "t",
            evidence_thresholds = evidence_thresholds,
            schedules = schedules,
            search_targets = sequential_search_targets,
            tags = c("t", "tbf01", "sequential"),
            rationale = paste(
                "Sequential t-test BF01 fixture summaries for dense and",
                "stress look schedules. Strict package-reference checks are",
                "curated separately because long exact ptbf01seq() runs are",
                "expensive."))
    )
}
