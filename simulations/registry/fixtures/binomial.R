bfpwr_sim_binomial_fixture_schedules_core <- function() {
    bfpwr_sim_z_fixture_schedules_core()
}

bfpwr_sim_fixture_specs_binomial <- function() {
    evidence_thresholds <- bfpwr_sim_fixture_evidence_thresholds_core()
    search_targets <- bfpwr_sim_fixture_search_targets_core()
    sequential_search_targets <- bfpwr_sim_search_targets(
        evidence_thresholds = evidence_thresholds,
        target_prob = numeric(0),
        evidence = c("H1", "H0"))
    schedules <- bfpwr_sim_binomial_fixture_schedules_core()

    list(
        bfpwr_sim_fixed_fixture_spec(
            fixture_set_id = "binom-binbf01-point-fixed-core-v1",
            family = "binomial",
            bf_types = "point",
            evidence_thresholds = evidence_thresholds,
            search_targets = search_targets,
            tags = c("binomial", "binbf01", "point", "fixed"),
            rationale = paste(
                "Fixed-look binomial point-null BF01 fixture summaries for",
                "threshold and sample-size search checks.")),
        bfpwr_sim_fixed_fixture_spec(
            fixture_set_id = "binom-binbf01-direction-fixed-core-v1",
            family = "binomial",
            bf_types = "direction",
            evidence_thresholds = evidence_thresholds,
            search_targets = search_targets,
            tags = c("binomial", "binbf01", "direction", "fixed"),
            rationale = paste(
                "Fixed-look binomial directional BF01 fixture summaries for",
                "threshold and sample-size search checks.")),
        bfpwr_sim_sequential_fixture_spec(
            fixture_set_id = "binom-binbf01-point-sequential-core-v1",
            family = "binomial",
            bf_types = "point",
            evidence_thresholds = evidence_thresholds,
            schedules = schedules,
            search_targets = sequential_search_targets,
            tags = c("binomial", "binbf01", "point", "sequential"),
            rationale = paste(
                "Sequential binomial point-null BF01 fixture summaries for",
                "dense and stress look schedules.")),
        bfpwr_sim_sequential_fixture_spec(
            fixture_set_id = "binom-binbf01-direction-sequential-core-v1",
            family = "binomial",
            bf_types = "direction",
            evidence_thresholds = evidence_thresholds,
            schedules = schedules,
            search_targets = sequential_search_targets,
            tags = c("binomial", "binbf01", "direction", "sequential"),
            rationale = paste(
                "Sequential binomial directional BF01 fixture summaries for",
                "dense and stress look schedules."))
    )
}
