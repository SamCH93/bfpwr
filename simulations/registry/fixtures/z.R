bfpwr_sim_z_fixture_schedules_core <- function() {
    comparison_schedules <- lapply(2:20, function(n_looks) {
        bfpwr_sim_schedule(
            paste0("start20-by10-looks",
                   sprintf("%02d", n_looks)),
            start_n = 20,
            increment = 10,
            n_looks = n_looks,
            schedule_family_id = "start20-by10",
            tags = c("short-compatible", "comparison-grid"))
    })
    stress_schedules <- list(
        bfpwr_sim_schedule(
            "start10-by10-looks50",
            start_n = 10, increment = 10, n_looks = 50,
            schedule_family_id = "start10-by10",
            tags = c("short-compatible", "stress")),
        bfpwr_sim_schedule(
            "start10-by10-looks100",
            start_n = 10, increment = 10, n_looks = 100,
            schedule_family_id = "start10-by10",
            tags = c("long-only", "stress")),
        bfpwr_sim_schedule(
            "start10-by10-looks200",
            start_n = 10, increment = 10, n_looks = 200,
            schedule_family_id = "start10-by10",
            tags = c("long-only", "stress"))
    )
    c(comparison_schedules, stress_schedules)
}

bfpwr_sim_fixture_specs_z <- function() {
    evidence_thresholds <- bfpwr_sim_fixture_evidence_thresholds_core()
    search_targets <- bfpwr_sim_fixture_search_targets_core()
    sequential_search_targets <- bfpwr_sim_search_targets(
        evidence_thresholds = evidence_thresholds,
        target_prob = numeric(0),
        evidence = c("H1", "H0"))
    schedules <- bfpwr_sim_z_fixture_schedules_core()

    list(
        bfpwr_sim_fixed_fixture_spec(
            fixture_set_id = "z-bf01-fixed-core-v1",
            family = "z",
            bf_types = "normal",
            evidence_thresholds = evidence_thresholds,
            search_targets = search_targets,
            tags = c("z", "bf01", "fixed"),
            rationale = paste(
                "Fixed-look z-test BF01 fixture summaries for symmetric",
                "evidence thresholds and sample-size search checks.")),
        bfpwr_sim_fixed_fixture_spec(
            fixture_set_id = "z-nmbf01-fixed-core-v1",
            family = "z",
            bf_types = "moment",
            evidence_thresholds = evidence_thresholds,
            search_targets = search_targets,
            tags = c("z", "nmbf01", "fixed"),
            rationale = paste(
                "Fixed-look z-test moment-prior fixture summaries for",
                "threshold and sample-size search checks.")),
        bfpwr_sim_fixed_fixture_spec(
            fixture_set_id = "z-dirbf01-fixed-core-v1",
            family = "z",
            bf_types = "directional",
            evidence_thresholds = evidence_thresholds,
            search_targets = search_targets,
            tags = c("z", "dirbf01", "fixed"),
            rationale = paste(
                "Fixed-look z-test directional BF fixture summaries for",
                "threshold and sample-size search checks.")),
        bfpwr_sim_sequential_fixture_spec(
            fixture_set_id = "z-bf01-sequential-core-v1",
            family = "z",
            bf_types = "normal",
            evidence_thresholds = evidence_thresholds,
            schedules = schedules,
            search_targets = sequential_search_targets,
            tags = c("z", "bf01", "sequential"),
            rationale = paste(
                "Sequential z-test BF01 fixture summaries for dense and",
                "coarse look schedules.")),
        bfpwr_sim_sequential_fixture_spec(
            fixture_set_id = "z-nmbf01-sequential-core-v1",
            family = "z",
            bf_types = "moment",
            evidence_thresholds = evidence_thresholds,
            schedules = schedules,
            search_targets = sequential_search_targets,
            tags = c("z", "nmbf01", "sequential"),
            rationale = paste(
                "Sequential z-test moment-prior fixture summaries for dense",
                "and coarse look schedules.")),
        bfpwr_sim_sequential_fixture_spec(
            fixture_set_id = "z-dirbf01-sequential-core-v1",
            family = "z",
            bf_types = "directional",
            evidence_thresholds = evidence_thresholds,
            schedules = schedules,
            search_targets = sequential_search_targets,
            tags = c("z", "dirbf01", "sequential"),
            rationale = paste(
                "Sequential z-test directional BF fixture summaries for dense",
                "and coarse look schedules."))
    )
}
