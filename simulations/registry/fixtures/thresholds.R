bfpwr_sim_fixture_evidence_thresholds_core <- function() {
    c(3, 10, 30)
}

bfpwr_sim_fixture_target_probabilities_core <- function() {
    c(0.1, 0.3, 0.8, 0.9, 0.95)
}

bfpwr_sim_fixture_search_targets_core <- function() {
    bfpwr_sim_search_targets(
        evidence_thresholds = bfpwr_sim_fixture_evidence_thresholds_core(),
        target_prob = bfpwr_sim_fixture_target_probabilities_core(),
        evidence = c("H1", "H0"))
}
