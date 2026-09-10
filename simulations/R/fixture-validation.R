bfpwr_sim_near <- function(x, y, tol = 1e-12) {
    both_na <- is.na(x) & is.na(y)
    out <- rep(FALSE, length(x))
    finite <- is.finite(x) & is.finite(y)
    out[finite] <- abs(x[finite] - y[finite]) <= tol
    out[both_na] <- TRUE
    out
}

bfpwr_sim_integerish <- function(x, tol = 1e-8) {
    is.finite(x) & abs(x - round(x)) <= tol
}

bfpwr_sim_row_key <- function(x, cols) {
    cols <- cols[cols %in% names(x)]
    if (length(cols) == 0) {
        return(rep("", nrow(x)))
    }
    do.call(paste, c(x[cols], sep = "\r"))
}

bfpwr_sim_failure <- function(rule,
                              n_failed,
                              details,
                              severity = "error") {
    data.frame(
        severity = severity,
        rule = rule,
        n_failed = as.integer(n_failed),
        details = details,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_bind_failures <- function(...) {
    x <- list(...)
    x <- x[vapply(x, function(y) is.data.frame(y) && nrow(y) > 0,
                  logical(1))]
    if (length(x) == 0) {
        return(data.frame(severity = character(),
                          rule = character(),
                          n_failed = integer(),
                          details = character()))
    }
    do.call(rbind, x)
}

bfpwr_sim_check_required_columns <- function(x, required, table_name) {
    missing <- setdiff(required, names(x))
    if (length(missing) == 0) {
        return(bfpwr_sim_bind_failures())
    }
    bfpwr_sim_failure(
        rule = paste0(table_name, ":required_columns"),
        n_failed = length(missing),
        details = paste("missing columns:", paste(missing, collapse = ", "))
    )
}

bfpwr_sim_fixture_manifest_status <- function(fixture_dir) {
    manifest_file <- file.path(fixture_dir, "manifest.csv")
    if (!file.exists(manifest_file)) {
        return(data.frame(file = "manifest.csv",
                          exists = FALSE,
                          expected_bytes = NA_real_,
                          actual_bytes = NA_real_,
                          bytes_ok = FALSE,
                          expected_sha256 = NA_character_,
                          actual_sha256 = NA_character_,
                          sha256_ok = FALSE,
                          stringsAsFactors = FALSE))
    }
    manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
    required <- c("file", "bytes", "sha256")
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0) {
        stop("fixture manifest is missing columns: ",
             paste(missing, collapse = ", "))
    }

    paths <- file.path(fixture_dir, manifest$file)
    exists <- file.exists(paths)
    actual_bytes <- rep(NA_real_, length(paths))
    actual_sha256 <- rep(NA_character_, length(paths))
    if (any(exists)) {
        actual_bytes[exists] <- file.info(paths[exists])$size
        actual_sha256[exists] <- unname(tools::sha256sum(paths[exists]))
    }

    data.frame(
        file = manifest$file,
        exists = exists,
        expected_bytes = manifest$bytes,
        actual_bytes = actual_bytes,
        bytes_ok = exists & actual_bytes == manifest$bytes,
        expected_sha256 = manifest$sha256,
        actual_sha256 = actual_sha256,
        sha256_ok = exists & actual_sha256 == manifest$sha256,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_validate_fixture_manifest <- function(fixture_dir,
                                                expected_files = NULL) {
    status <- bfpwr_sim_fixture_manifest_status(fixture_dir)
    manifest_file <- file.path(fixture_dir, "manifest.csv")
    manifest <- if (file.exists(manifest_file)) {
        utils::read.csv(manifest_file, stringsAsFactors = FALSE)
    } else {
        data.frame(file = character(), bytes = numeric(),
                   sha256 = character())
    }
    unsafe <- grepl("^([A-Za-z]:)?[/\\\\]", manifest$file) |
        grepl("(^|[/\\\\])[.][.]([/\\\\]|$)", manifest$file)
    duplicated_files <- duplicated(manifest$file)
    missing_expected <- if (is.null(expected_files)) {
        character()
    } else {
        setdiff(expected_files, manifest$file)
    }
    extra_files <- if (is.null(expected_files)) {
        character()
    } else {
        setdiff(manifest$file, expected_files)
    }
    bfpwr_sim_bind_failures(
        if (any(unsafe)) {
            bfpwr_sim_failure(
                "manifest:relative_safe_paths",
                sum(unsafe),
                "manifest file paths must be relative and stay inside the fixture directory"
            )
        },
        if (any(duplicated_files)) {
            bfpwr_sim_failure(
                "manifest:unique_files",
                sum(duplicated_files),
                "manifest file paths must be unique"
            )
        },
        if (length(missing_expected) > 0) {
            bfpwr_sim_failure(
                "manifest:expected_files_present",
                length(missing_expected),
                paste("missing expected manifest entries:",
                      paste(missing_expected, collapse = ", "))
            )
        },
        if (length(extra_files) > 0) {
            bfpwr_sim_failure(
                "manifest:no_extra_files",
                length(extra_files),
                paste("unexpected manifest entries:",
                      paste(extra_files, collapse = ", "))
            )
        },
        if (any(!status$exists)) {
            bfpwr_sim_failure(
                "manifest:files_exist",
                sum(!status$exists),
                paste("missing files:",
                      paste(status$file[!status$exists], collapse = ", "))
            )
        },
        if (any(status$exists & !status$bytes_ok)) {
            bfpwr_sim_failure(
                "manifest:byte_sizes_match",
                sum(status$exists & !status$bytes_ok),
                paste("byte-size mismatches:",
                      paste(status$file[status$exists & !status$bytes_ok],
                            collapse = ", "))
            )
        },
        if (any(status$exists & !status$sha256_ok)) {
            bfpwr_sim_failure(
                "manifest:sha256_hashes_match",
                sum(status$exists & !status$sha256_ok),
                paste("SHA256 mismatches:",
                      paste(status$file[status$exists & !status$sha256_ok],
                            collapse = ", "))
            )
        }
    )
}

bfpwr_sim_expected_fixture_files <- function(mode = c("fixed",
                                                      "sequential")) {
    mode <- match.arg(mode)
    if (mode == "fixed") {
        return(c("spec.rds", "tail-summary.rds", "decision-summary.rds",
                 "search-summary.rds", "overview.csv"))
    }
    c("spec.rds", "cumulative-summary.rds", "final-summary.rds",
      "search-summary.rds", "overview.csv")
}

bfpwr_sim_validate_fixture_overview <- function(fixture, mode) {
    overview_file <- file.path(fixture$fixture_dir, "overview.csv")
    if (!file.exists(overview_file)) {
        return(bfpwr_sim_failure(
            "overview:exists", 1L,
            "overview.csv must exist in each fixture summary directory"))
    }
    overview <- utils::read.csv(overview_file, stringsAsFactors = FALSE)
    if (nrow(overview) != 1) {
        return(bfpwr_sim_failure(
            "overview:single_row", nrow(overview),
            "overview.csv must contain exactly one row"))
    }
    required <- if (mode == "fixed") {
        c("fixture_set_id", "family", "mode", "bf_priors", "tail_rows",
          "decision_rows", "search_rows")
    } else {
        c("fixture_set_id", "family", "mode", "bf_priors",
          "cumulative_rows", "final_rows", "search_rows")
    }
    failures <- bfpwr_sim_check_required_columns(overview, required,
                                                 "overview")
    if (nrow(failures) > 0) return(failures)

    common_bad <- !identical(overview$fixture_set_id[[1]],
                             fixture$spec$fixture_set_id) ||
        !identical(overview$family[[1]], fixture$spec$family) ||
        !identical(overview$mode[[1]], fixture$spec$mode)
    if (mode == "fixed") {
        row_bad <- overview$bf_priors[[1]] !=
            length(unique(fixture$tail_summary$bf_prior_id)) ||
            overview$tail_rows[[1]] != nrow(fixture$tail_summary) ||
            overview$decision_rows[[1]] != nrow(fixture$decision_summary) ||
            overview$search_rows[[1]] != nrow(fixture$search_summary)
    } else {
        row_bad <- overview$bf_priors[[1]] !=
            length(unique(fixture$cumulative_summary$bf_prior_id)) ||
            overview$cumulative_rows[[1]] !=
                nrow(fixture$cumulative_summary) ||
            overview$final_rows[[1]] != nrow(fixture$final_summary) ||
            overview$search_rows[[1]] != nrow(fixture$search_summary)
    }
    bfpwr_sim_bind_failures(
        if (common_bad) {
            bfpwr_sim_failure(
                "overview:fixture_metadata",
                1L,
                "overview fixture_set_id, family, and mode must match spec.rds"
            )
        },
        if (row_bad) {
            bfpwr_sim_failure(
                "overview:row_counts",
                1L,
                "overview row counts must match loaded summary tables"
            )
        }
    )
}

bfpwr_sim_validate_probability_counts <- function(x,
                                                  table_name,
                                                  prob_col,
                                                  count_col,
                                                  nsim_col = "nsim",
                                                  mcse_col = NULL,
                                                  tol = 1e-12) {
    required <- c(prob_col, count_col, nsim_col)
    if (!is.null(mcse_col)) required <- c(required, mcse_col)
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
        return(bfpwr_sim_check_required_columns(x, required, table_name))
    }

    prob <- x[[prob_col]]
    count <- x[[count_col]]
    nsim <- x[[nsim_col]]
    expected_prob <- count / nsim
    failures <- list()

    bad_nsim <- !bfpwr_sim_integerish(nsim) | nsim <= 0
    if (any(bad_nsim, na.rm = TRUE)) {
        failures[[length(failures) + 1L]] <- bfpwr_sim_failure(
            paste0(table_name, ":", nsim_col, "_positive_integer"),
            sum(bad_nsim, na.rm = TRUE),
            paste(nsim_col, "must be a positive integer")
        )
    }

    bad_count <- !bfpwr_sim_integerish(count) | count < 0 | count > nsim
    if (any(bad_count, na.rm = TRUE)) {
        failures[[length(failures) + 1L]] <- bfpwr_sim_failure(
            paste0(table_name, ":", count_col, "_valid_count"),
            sum(bad_count, na.rm = TRUE),
            paste(count_col, "must be an integer between 0 and", nsim_col)
        )
    }

    bad_prob <- !is.finite(prob) | prob < 0 | prob > 1
    if (any(bad_prob, na.rm = TRUE)) {
        failures[[length(failures) + 1L]] <- bfpwr_sim_failure(
            paste0(table_name, ":", prob_col, "_probability"),
            sum(bad_prob, na.rm = TRUE),
            paste(prob_col, "must be finite and in [0, 1]")
        )
    }

    bad_prob_count <- !bfpwr_sim_near(prob, expected_prob, tol = tol)
    if (any(bad_prob_count, na.rm = TRUE)) {
        failures[[length(failures) + 1L]] <- bfpwr_sim_failure(
            paste0(table_name, ":", prob_col, "_matches_", count_col),
            sum(bad_prob_count, na.rm = TRUE),
            paste(prob_col, "must equal", count_col, "/", nsim_col)
        )
    }

    if (!is.null(mcse_col)) {
        mcse <- x[[mcse_col]]
        expected_mcse <- bfpwr_sim_mcse(prob, nsim)
        bad_mcse <- !is.finite(mcse) |
            !bfpwr_sim_near(mcse, expected_mcse, tol = tol)
        if (any(bad_mcse, na.rm = TRUE)) {
            failures[[length(failures) + 1L]] <- bfpwr_sim_failure(
                paste0(table_name, ":", mcse_col, "_matches_binomial_mcse"),
                sum(bad_mcse, na.rm = TRUE),
                paste(mcse_col, "must equal sqrt(p * (1 - p) / nsim)")
            )
        }
    }

    bfpwr_sim_bind_failures(failures)
}

bfpwr_sim_validate_log_bf_status_counts <- function(x, table_name,
                                                    multiplier_col = NULL) {
    required <- c("nsim", "n_finite_log_bf01", "n_infinite_log_bf01",
                  "n_nan_log_bf01", "n_na_log_bf01")
    failures <- bfpwr_sim_check_required_columns(x, required, table_name)
    if (nrow(failures) > 0) return(failures)

    expected <- x$nsim
    if (!is.null(multiplier_col) && multiplier_col %in% names(x)) {
        expected <- expected * x[[multiplier_col]]
    }
    total <- x$n_finite_log_bf01 + x$n_infinite_log_bf01 +
        x$n_nan_log_bf01 + x$n_na_log_bf01

    bfpwr_sim_bind_failures(
        if (any(total != expected, na.rm = TRUE)) {
            bfpwr_sim_failure(
                paste0(table_name, ":log_bf_status_counts_sum"),
                sum(total != expected, na.rm = TRUE),
                "finite/infinite/NaN/NA log-BF counts must sum to the expected number of evaluations"
            )
        },
        if (any(x$n_nan_log_bf01 != 0, na.rm = TRUE)) {
            bfpwr_sim_failure(
                paste0(table_name, ":no_nan_log_bf01"),
                sum(x$n_nan_log_bf01 != 0, na.rm = TRUE),
                "NaN log_bf01 values are not allowed in finalized fixtures"
            )
        },
        if (any(x$n_na_log_bf01 != 0, na.rm = TRUE)) {
            bfpwr_sim_failure(
                paste0(table_name, ":no_na_log_bf01"),
                sum(x$n_na_log_bf01 != 0, na.rm = TRUE),
                "NA log_bf01 values are not allowed in finalized fixtures"
            )
        }
    )
}

bfpwr_sim_validate_fixed_tail_summary <- function(tail_summary) {
    required <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "threshold_id",
                  "evidence_threshold", "tail", "threshold",
                  "log_threshold", "n", "nsim", "n_event", "prob", "mcse",
                  "n_finite_log_bf01", "n_infinite_log_bf01",
                  "n_nan_log_bf01", "n_na_log_bf01", "source_chunk_count")
    failures <- bfpwr_sim_check_required_columns(tail_summary, required,
                                                 "tail_summary")
    if (nrow(failures) > 0) return(failures)

    expected_threshold <- ifelse(tail_summary$tail == "H1",
                                 1 / tail_summary$evidence_threshold,
                                 tail_summary$evidence_threshold)
    key_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "threshold_id",
                  "tail", "n")
    duplicate_keys <- duplicated(bfpwr_sim_row_key(tail_summary, key_cols))
    bad_threshold <- !tail_summary$tail %in% c("H1", "H0") |
        !is.finite(tail_summary$evidence_threshold) |
        tail_summary$evidence_threshold <= 1 |
        !bfpwr_sim_near(tail_summary$threshold, expected_threshold) |
        !bfpwr_sim_near(tail_summary$log_threshold,
                        log(tail_summary$threshold))

    log_count_failures <- bfpwr_sim_validate_log_bf_status_counts(
        tail_summary, "tail_summary")
    count_failures <- bfpwr_sim_validate_probability_counts(
        tail_summary, "tail_summary", "prob", "n_event",
        mcse_col = "mcse")

    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "tail", "n")
    group_key <- bfpwr_sim_row_key(tail_summary, group_cols)
    ord <- order(group_key, tail_summary$evidence_threshold)
    ordered_group <- group_key[ord]
    ordered_prob <- tail_summary$prob[ord]
    same_group <- ordered_group[-1L] == ordered_group[-length(ordered_group)]
    monotone_bad_rows <- same_group & diff(ordered_prob) > 1e-12
    monotone_bad_groups <- unique(ordered_group[-1L][monotone_bad_rows])

    bfpwr_sim_bind_failures(
        if (any(duplicate_keys)) {
            bfpwr_sim_failure(
                "tail_summary:unique_keys",
                sum(duplicate_keys),
                "tail summary key rows must be unique"
            )
        },
        if (any(bad_threshold, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "tail_summary:threshold_metadata",
                sum(bad_threshold, na.rm = TRUE),
                "tail, threshold, evidence_threshold, and log_threshold must agree"
            )
        },
        log_count_failures,
        count_failures,
        if (length(monotone_bad_groups) > 0) {
            bfpwr_sim_failure(
                "tail_summary:threshold_monotonicity",
                length(monotone_bad_groups),
                "event probability must not increase as the evidence threshold becomes more stringent"
            )
        }
    )
}

bfpwr_sim_validate_fixed_decision_summary <- function(decision_summary,
                                                      tail_summary) {
    required <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "n",
                  "threshold_pair_id", "evidence_threshold", "k1", "k0",
                  "pH1", "pH0", "pInc", "nH1", "nH0", "nInc", "nsim",
                  "mcse_pH1", "mcse_pH0", "mcse_pInc")
    failures <- bfpwr_sim_check_required_columns(decision_summary, required,
                                                 "decision_summary")
    if (nrow(failures) > 0) return(failures)

    count_failures <- bfpwr_sim_bind_failures(
        bfpwr_sim_validate_probability_counts(
            decision_summary, "decision_summary", "pH1", "nH1",
            mcse_col = "mcse_pH1"),
        bfpwr_sim_validate_probability_counts(
            decision_summary, "decision_summary", "pH0", "nH0",
            mcse_col = "mcse_pH0"),
        bfpwr_sim_validate_probability_counts(
            decision_summary, "decision_summary", "pInc", "nInc",
            mcse_col = "mcse_pInc")
    )

    bad_threshold <- !bfpwr_sim_near(decision_summary$k1,
                                     1 / decision_summary$evidence_threshold) |
        !bfpwr_sim_near(decision_summary$k0,
                        decision_summary$evidence_threshold)
    key_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "n",
                  "threshold_pair_id")
    duplicate_keys <- duplicated(bfpwr_sim_row_key(decision_summary,
                                                   key_cols))
    bad_counts <- decision_summary$nH1 + decision_summary$nH0 +
        decision_summary$nInc != decision_summary$nsim
    bad_probs <- !bfpwr_sim_near(decision_summary$pH1 +
                                     decision_summary$pH0 +
                                     decision_summary$pInc,
                                 rep(1, nrow(decision_summary)))

    relation_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                       "design_case_id", "look_grid_name", "n",
                       "evidence_threshold")
    dkey <- bfpwr_sim_row_key(decision_summary, relation_cols)
    h1 <- tail_summary[tail_summary$tail == "H1", , drop = FALSE]
    h0 <- tail_summary[tail_summary$tail == "H0", , drop = FALSE]
    h1_match <- match(dkey, bfpwr_sim_row_key(h1, relation_cols))
    h0_match <- match(dkey, bfpwr_sim_row_key(h0, relation_cols))
    missing_tail <- is.na(h1_match) | is.na(h0_match)
    relation_bad <- rep(FALSE, nrow(decision_summary))
    relation_bad[!missing_tail] <-
        !bfpwr_sim_near(decision_summary$pH1[!missing_tail],
                        h1$prob[h1_match[!missing_tail]]) |
        !bfpwr_sim_near(decision_summary$pH0[!missing_tail],
                        h0$prob[h0_match[!missing_tail]]) |
        decision_summary$nH1[!missing_tail] !=
            h1$n_event[h1_match[!missing_tail]] |
        decision_summary$nH0[!missing_tail] !=
            h0$n_event[h0_match[!missing_tail]]

    bfpwr_sim_bind_failures(
        count_failures,
        if (any(duplicate_keys)) {
            bfpwr_sim_failure(
                "decision_summary:unique_keys",
                sum(duplicate_keys),
                "decision summary key rows must be unique"
            )
        },
        if (any(bad_threshold, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "decision_summary:threshold_pair_metadata",
                sum(bad_threshold, na.rm = TRUE),
                "k1 and k0 must equal 1 / evidence_threshold and evidence_threshold"
            )
        },
        if (any(bad_counts, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "decision_summary:counts_sum_to_nsim",
                sum(bad_counts, na.rm = TRUE),
                "nH1 + nH0 + nInc must equal nsim"
            )
        },
        if (any(bad_probs, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "decision_summary:probabilities_sum_to_one",
                sum(bad_probs, na.rm = TRUE),
                "pH1 + pH0 + pInc must equal 1"
            )
        },
        if (any(missing_tail)) {
            bfpwr_sim_failure(
                "decision_summary:matching_tail_rows_exist",
                sum(missing_tail),
                "each decision row must have matching H1 and H0 tail rows"
            )
        },
        if (any(relation_bad)) {
            bfpwr_sim_failure(
                "decision_summary:matches_tail_summary",
                sum(relation_bad),
                "decision H1/H0 probabilities and counts must match tail_summary"
            )
        }
    )
}

bfpwr_sim_validate_fixed_search_summary <- function(search_summary,
                                                    tail_summary) {
    required <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "tail",
                  "evidence_threshold", "threshold", "log_threshold",
                  "search_target_id", "target_prob", "criterion",
                  "achieved", "n_found", "prob_found", "mcse_found",
                  "n_previous", "prob_previous", "best_n", "best_prob",
                  "max_n_searched")
    failures <- bfpwr_sim_check_required_columns(search_summary, required,
                                                 "search_summary")
    if (nrow(failures) > 0 || nrow(search_summary) == 0) return(failures)

    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "tail",
                    "evidence_threshold", "threshold", "log_threshold")
    groups <- split(tail_summary, bfpwr_sim_row_key(tail_summary, group_cols))
    key_cols <- c(group_cols, "search_target_id")
    duplicate_keys <- duplicated(bfpwr_sim_row_key(search_summary, key_cols))
    bad <- rep(FALSE, nrow(search_summary))
    missing_curve <- rep(FALSE, nrow(search_summary))
    for (i in seq_len(nrow(search_summary))) {
        key <- bfpwr_sim_row_key(search_summary[i, , drop = FALSE],
                                 group_cols)
        curve <- groups[[key]]
        if (is.null(curve)) {
            missing_curve[[i]] <- TRUE
            next
        }
        curve <- curve[order(curve$n), , drop = FALSE]
        hit <- which(curve$prob >= search_summary$target_prob[[i]])
        achieved <- length(hit) > 0
        first <- if (achieved) hit[[1]] else NA_integer_
        previous <- if (achieved && first > 1) first - 1L else NA_integer_
        best <- which.max(curve$prob)

        expected <- list(
            criterion = "first_n_with_prob_ge_target",
            achieved = achieved,
            n_found = if (achieved) curve$n[[first]] else NA_integer_,
            prob_found = if (achieved) curve$prob[[first]] else NA_real_,
            mcse_found = if (achieved) curve$mcse[[first]] else NA_real_,
            n_previous = if (is.na(previous)) NA_integer_ else curve$n[[previous]],
            prob_previous = if (is.na(previous)) NA_real_ else curve$prob[[previous]],
            best_n = curve$n[[best]],
            best_prob = curve$prob[[best]],
            max_n_searched = max(curve$n)
        )
        bad[[i]] <- !identical(search_summary$criterion[[i]],
                               expected$criterion) ||
            !identical(as.logical(search_summary$achieved[[i]]),
                       expected$achieved) ||
            !bfpwr_sim_near(search_summary$n_found[[i]], expected$n_found) ||
            !bfpwr_sim_near(search_summary$prob_found[[i]],
                            expected$prob_found) ||
            !bfpwr_sim_near(search_summary$mcse_found[[i]],
                            expected$mcse_found) ||
            !bfpwr_sim_near(search_summary$n_previous[[i]],
                            expected$n_previous) ||
            !bfpwr_sim_near(search_summary$prob_previous[[i]],
                            expected$prob_previous) ||
            !bfpwr_sim_near(search_summary$best_n[[i]], expected$best_n) ||
            !bfpwr_sim_near(search_summary$best_prob[[i]],
                            expected$best_prob) ||
            !bfpwr_sim_near(search_summary$max_n_searched[[i]],
                            expected$max_n_searched)
    }

    bfpwr_sim_bind_failures(
        if (any(duplicate_keys)) {
            bfpwr_sim_failure(
                "search_summary:unique_keys",
                sum(duplicate_keys),
                "fixed search summary key rows must be unique"
            )
        },
        if (any(missing_curve)) {
            bfpwr_sim_failure(
                "search_summary:matching_tail_curve_exists",
                sum(missing_curve),
                "each fixed search row must match a tail-summary curve"
            )
        },
        if (any(bad)) {
            bfpwr_sim_failure(
                "search_summary:first_crossing_matches_tail_curve",
                sum(bad),
                "fixed search rows must encode the first probability crossing"
            )
        }
    )
}

bfpwr_sim_validate_fixed_fixture_deterministic <- function(fixture) {
    bfpwr_sim_bind_failures(
        bfpwr_sim_validate_fixed_tail_summary(fixture$tail_summary),
        bfpwr_sim_validate_fixed_decision_summary(fixture$decision_summary,
                                                  fixture$tail_summary),
        bfpwr_sim_validate_fixed_search_summary(fixture$search_summary,
                                                fixture$tail_summary)
    )
}

bfpwr_sim_validate_sequential_cumulative_summary <- function(cumulative) {
    required <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "schedule_id",
                  "schedule_family_id", "start_n", "increment", "n_looks",
                  "max_n", "threshold_pair_id", "evidence_threshold",
                  "k1", "k0", "look", "n", "nsim",
                  "stop_H1_at_look", "stop_H0_at_look", "cum_pH1",
                  "cum_pH0", "cum_pInc", "mcse_cum_pH1",
                  "mcse_cum_pH0", "mcse_cum_pInc", "EN_to_look",
                  "n_finite_log_bf01", "n_infinite_log_bf01",
                  "n_nan_log_bf01", "n_na_log_bf01", "source_chunk_count")
    failures <- bfpwr_sim_check_required_columns(cumulative, required,
                                                 "cumulative_summary")
    if (nrow(failures) > 0 || nrow(cumulative) == 0) return(failures)

    count_failures <- bfpwr_sim_bind_failures()

    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "schedule_id",
                    "threshold_pair_id", "evidence_threshold")
    groups <- split(cumulative, bfpwr_sim_row_key(cumulative, group_cols))
    duplicate_keys <- duplicated(bfpwr_sim_row_key(
        cumulative, c(group_cols, "look")))
    group_bad <- list()
    add_group_failure <- function(rule, bad, details) {
        if (any(bad)) {
            group_bad[[length(group_bad) + 1L]] <<-
                bfpwr_sim_failure(rule, sum(bad), details)
        }
    }

    bad_counts <- bad_probs <- bad_mcse <- bad_monotone <- bad_en <- logical()
    bad_log_counts <- logical()
    for (g in groups) {
        g <- g[order(g$look), , drop = FALSE]
        cum_h1 <- cumsum(g$stop_H1_at_look)
        cum_h0 <- cumsum(g$stop_H0_at_look)
        cum_inc <- g$nsim[[1]] - cum_h1 - cum_h0
        expected_eval_count <- g$nsim * g$n_looks
        status_total <- g$n_finite_log_bf01 + g$n_infinite_log_bf01 +
            g$n_nan_log_bf01 + g$n_na_log_bf01
        bad_counts <- c(bad_counts,
                        any(cum_inc < 0) ||
                            any(g$stop_H1_at_look < 0) ||
                            any(g$stop_H0_at_look < 0))
        bad_probs <- c(bad_probs,
                       !all(bfpwr_sim_near(g$cum_pH1, cum_h1 / g$nsim)) ||
                           !all(bfpwr_sim_near(g$cum_pH0, cum_h0 / g$nsim)) ||
                           !all(bfpwr_sim_near(g$cum_pInc, cum_inc / g$nsim)))
        bad_mcse <- c(bad_mcse,
                      !all(bfpwr_sim_near(g$mcse_cum_pH1,
                                          bfpwr_sim_mcse(g$cum_pH1, g$nsim))) ||
                          !all(bfpwr_sim_near(g$mcse_cum_pH0,
                                              bfpwr_sim_mcse(g$cum_pH0, g$nsim))) ||
                          !all(bfpwr_sim_near(g$mcse_cum_pInc,
                                              bfpwr_sim_mcse(g$cum_pInc, g$nsim))))
        bad_monotone <- c(bad_monotone,
                          any(diff(g$cum_pH1) < -1e-12) ||
                              any(diff(g$cum_pH0) < -1e-12) ||
                              any(diff(g$cum_pInc) > 1e-12))
        bad_en <- c(bad_en, any(!is.finite(g$EN_to_look)) ||
                        any(diff(g$EN_to_look) < -1e-12))
        bad_log_counts <- c(bad_log_counts,
                            any(status_total != expected_eval_count) ||
                                any(g$n_nan_log_bf01 != 0) ||
                                any(g$n_na_log_bf01 != 0))
    }
    add_group_failure("cumulative_summary:valid_cumulative_counts",
                      bad_counts,
                      "sequential cumulative counts must be non-negative")
    add_group_failure("cumulative_summary:probabilities_match_counts",
                      bad_probs,
                      "cum_pH1, cum_pH0, and cum_pInc must match cumulative counts")
    add_group_failure("cumulative_summary:mcse_matches_probabilities",
                      bad_mcse,
                      "cumulative MCSE columns must match binomial MCSEs")
    add_group_failure("cumulative_summary:cumulative_monotonicity",
                      bad_monotone,
                      "cum_pH1/cum_pH0 must be non-decreasing and cum_pInc non-increasing")
    add_group_failure("cumulative_summary:EN_to_look_valid",
                      bad_en,
                      "EN_to_look must be finite and non-decreasing")
    add_group_failure("cumulative_summary:log_bf_status_counts",
                      bad_log_counts,
                      "log-BF status counts must match nsim * n_looks with no NaN/NA values")

    bad_threshold <- !bfpwr_sim_near(cumulative$k1,
                                     1 / cumulative$evidence_threshold) |
        !bfpwr_sim_near(cumulative$k0, cumulative$evidence_threshold)
    bfpwr_sim_bind_failures(
        count_failures,
        if (any(duplicate_keys)) {
            bfpwr_sim_failure(
                "cumulative_summary:unique_keys",
                sum(duplicate_keys),
                "cumulative summary key rows must be unique"
            )
        },
        bfpwr_sim_bind_failures(group_bad),
        if (any(bad_threshold, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "cumulative_summary:threshold_pair_metadata",
                sum(bad_threshold, na.rm = TRUE),
                "k1 and k0 must equal 1 / evidence_threshold and evidence_threshold"
            )
        }
    )
}

bfpwr_sim_validate_sequential_final_summary <- function(final, cumulative) {
    required <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "schedule_id",
                  "schedule_family_id", "start_n", "increment", "n_looks",
                  "max_n", "threshold_pair_id", "evidence_threshold",
                  "k1", "k0", "nsim", "pH1", "pH0", "pInc",
                  "nH1", "nH0", "nInc", "mcse_pH1", "mcse_pH0",
                  "mcse_pInc", "EN", "VarN", "q25_N", "median_N",
                  "q75_N", "p_stop_by_final", "n_finite_log_bf01",
                  "n_infinite_log_bf01", "n_nan_log_bf01", "n_na_log_bf01",
                  "source_chunk_count")
    failures <- bfpwr_sim_check_required_columns(final, required,
                                                 "final_summary")
    if (nrow(failures) > 0 || nrow(final) == 0) return(failures)

    count_failures <- bfpwr_sim_bind_failures(
        bfpwr_sim_validate_probability_counts(
            final, "final_summary", "pH1", "nH1", mcse_col = "mcse_pH1"),
        bfpwr_sim_validate_probability_counts(
            final, "final_summary", "pH0", "nH0", mcse_col = "mcse_pH0"),
        bfpwr_sim_validate_probability_counts(
            final, "final_summary", "pInc", "nInc", mcse_col = "mcse_pInc")
    )
    log_count_failures <- bfpwr_sim_validate_log_bf_status_counts(
        final, "final_summary", multiplier_col = "n_looks")

    bad_counts <- final$nH1 + final$nH0 + final$nInc != final$nsim
    key_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "schedule_id",
                  "threshold_pair_id", "evidence_threshold")
    duplicate_keys <- duplicated(bfpwr_sim_row_key(final, key_cols))
    bad_probs <- !bfpwr_sim_near(final$pH1 + final$pH0 + final$pInc,
                                 rep(1, nrow(final)))
    bad_threshold <- !bfpwr_sim_near(final$k1,
                                     1 / final$evidence_threshold) |
        !bfpwr_sim_near(final$k0, final$evidence_threshold)
    bad_stop <- !bfpwr_sim_near(final$p_stop_by_final, 1 - final$pInc)
    bad_n <- !is.finite(final$EN) | final$EN < final$start_n |
        final$EN > final$max_n | !is.finite(final$VarN) | final$VarN < -1e-8 |
        final$q25_N > final$median_N | final$median_N > final$q75_N

    cgroups <- split(cumulative, bfpwr_sim_row_key(cumulative, key_cols))
    missing_cumulative <- relation_bad <- rep(FALSE, nrow(final))
    for (i in seq_len(nrow(final))) {
        key <- bfpwr_sim_row_key(final[i, , drop = FALSE], key_cols)
        g <- cgroups[[key]]
        if (is.null(g)) {
            missing_cumulative[[i]] <- TRUE
            next
        }
        g <- g[order(g$look), , drop = FALSE]
        last <- g[nrow(g), , drop = FALSE]
        stop_counts <- g$stop_H1_at_look + g$stop_H0_at_look
        n_inc <- final$nsim[[i]] - sum(stop_counts)
        stop_n_sum <- sum(stop_counts * g$n) + n_inc * final$max_n[[i]]
        stop_n2_sum <- sum(stop_counts * g$n^2) + n_inc * final$max_n[[i]]^2
        expected_en <- stop_n_sum / final$nsim[[i]]
        expected_var <- stop_n2_sum / final$nsim[[i]] - expected_en^2
        relation_bad[[i]] <- !bfpwr_sim_near(final$pH1[[i]],
                                             last$cum_pH1[[1]]) ||
            !bfpwr_sim_near(final$pH0[[i]], last$cum_pH0[[1]]) ||
            !bfpwr_sim_near(final$pInc[[i]], last$cum_pInc[[1]]) ||
            final$nH1[[i]] != sum(g$stop_H1_at_look) ||
            final$nH0[[i]] != sum(g$stop_H0_at_look) ||
            final$nInc[[i]] != n_inc ||
            !bfpwr_sim_near(final$EN[[i]], expected_en) ||
            !bfpwr_sim_near(final$VarN[[i]], expected_var, tol = 1e-8)
    }

    bfpwr_sim_bind_failures(
        count_failures,
        log_count_failures,
        if (any(duplicate_keys)) {
            bfpwr_sim_failure(
                "final_summary:unique_keys",
                sum(duplicate_keys),
                "final summary key rows must be unique"
            )
        },
        if (any(bad_counts, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "final_summary:counts_sum_to_nsim",
                sum(bad_counts, na.rm = TRUE),
                "nH1 + nH0 + nInc must equal nsim"
            )
        },
        if (any(bad_probs, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "final_summary:probabilities_sum_to_one",
                sum(bad_probs, na.rm = TRUE),
                "pH1 + pH0 + pInc must equal 1"
            )
        },
        if (any(bad_threshold, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "final_summary:threshold_pair_metadata",
                sum(bad_threshold, na.rm = TRUE),
                "k1 and k0 must equal 1 / evidence_threshold and evidence_threshold"
            )
        },
        if (any(bad_stop, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "final_summary:p_stop_by_final",
                sum(bad_stop, na.rm = TRUE),
                "p_stop_by_final must equal 1 - pInc"
            )
        },
        if (any(bad_n, na.rm = TRUE)) {
            bfpwr_sim_failure(
                "final_summary:sample_size_summaries_valid",
                sum(bad_n, na.rm = TRUE),
                "EN, VarN, and N quantiles must be internally valid"
            )
        },
        if (any(missing_cumulative)) {
            bfpwr_sim_failure(
                "final_summary:matching_cumulative_rows_exist",
                sum(missing_cumulative),
                "each final row must have matching cumulative rows"
            )
        },
        if (any(relation_bad)) {
            bfpwr_sim_failure(
                "final_summary:matches_cumulative_summary",
                sum(relation_bad),
                "final probabilities/counts/EN/VarN must match cumulative stopping counts"
            )
        }
    )
}

bfpwr_sim_validate_sequential_search_summary <- function(search_summary,
                                                         cumulative) {
    required <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                  "design_case_id", "look_grid_name", "schedule_id",
                  "schedule_family_id", "start_n", "increment", "n_looks",
                  "max_n", "threshold_pair_id", "evidence_threshold",
                  "k1", "k0", "search_target_id", "evidence",
                  "target_prob", "criterion", "achieved", "look_found",
                  "n_found", "prob_found", "mcse_found", "look_previous",
                  "n_previous", "prob_previous", "best_look", "best_n",
                  "best_prob", "max_n_searched")
    failures <- bfpwr_sim_check_required_columns(search_summary, required,
                                                 "search_summary")
    if (nrow(failures) > 0 || nrow(search_summary) == 0) return(failures)

    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "schedule_id",
                    "threshold_pair_id", "evidence_threshold")
    groups <- split(cumulative, bfpwr_sim_row_key(cumulative, group_cols))
    duplicate_keys <- duplicated(bfpwr_sim_row_key(
        search_summary, c(group_cols, "search_target_id", "evidence")))
    bad <- missing_curve <- rep(FALSE, nrow(search_summary))
    for (i in seq_len(nrow(search_summary))) {
        key <- bfpwr_sim_row_key(search_summary[i, , drop = FALSE],
                                 group_cols)
        curve <- groups[[key]]
        if (is.null(curve)) {
            missing_curve[[i]] <- TRUE
            next
        }
        curve <- curve[order(curve$look), , drop = FALSE]
        prob_col <- if (search_summary$evidence[[i]] == "H1") {
            "cum_pH1"
        } else {
            "cum_pH0"
        }
        mcse_col <- if (search_summary$evidence[[i]] == "H1") {
            "mcse_cum_pH1"
        } else {
            "mcse_cum_pH0"
        }
        prob <- curve[[prob_col]]
        hit <- which(prob >= search_summary$target_prob[[i]])
        achieved <- length(hit) > 0
        first <- if (achieved) hit[[1]] else NA_integer_
        previous <- if (achieved && first > 1) first - 1L else NA_integer_
        best <- which.max(prob)
        bad[[i]] <- !identical(search_summary$criterion[[i]],
                               "first_look_with_cumulative_prob_ge_target") ||
            !identical(as.logical(search_summary$achieved[[i]]), achieved) ||
            !bfpwr_sim_near(search_summary$look_found[[i]],
                            if (achieved) curve$look[[first]] else NA_integer_) ||
            !bfpwr_sim_near(search_summary$n_found[[i]],
                            if (achieved) curve$n[[first]] else NA_integer_) ||
            !bfpwr_sim_near(search_summary$prob_found[[i]],
                            if (achieved) prob[[first]] else NA_real_) ||
            !bfpwr_sim_near(search_summary$mcse_found[[i]],
                            if (achieved) curve[[mcse_col]][[first]] else NA_real_) ||
            !bfpwr_sim_near(search_summary$look_previous[[i]],
                            if (is.na(previous)) NA_integer_ else curve$look[[previous]]) ||
            !bfpwr_sim_near(search_summary$n_previous[[i]],
                            if (is.na(previous)) NA_integer_ else curve$n[[previous]]) ||
            !bfpwr_sim_near(search_summary$prob_previous[[i]],
                            if (is.na(previous)) NA_real_ else prob[[previous]]) ||
            !bfpwr_sim_near(search_summary$best_look[[i]],
                            curve$look[[best]]) ||
            !bfpwr_sim_near(search_summary$best_n[[i]], curve$n[[best]]) ||
            !bfpwr_sim_near(search_summary$best_prob[[i]], prob[[best]]) ||
            !bfpwr_sim_near(search_summary$max_n_searched[[i]], max(curve$n))
    }

    bfpwr_sim_bind_failures(
        if (any(duplicate_keys)) {
            bfpwr_sim_failure(
                "search_summary:unique_keys",
                sum(duplicate_keys),
                "sequential search summary key rows must be unique"
            )
        },
        if (any(missing_curve)) {
            bfpwr_sim_failure(
                "search_summary:matching_cumulative_curve_exists",
                sum(missing_curve),
                "each sequential search row must match a cumulative-summary curve"
            )
        },
        if (any(bad)) {
            bfpwr_sim_failure(
                "search_summary:first_crossing_matches_cumulative_curve",
                sum(bad),
                "sequential search rows must encode the first cumulative crossing"
            )
        }
    )
}

bfpwr_sim_validate_sequential_fixture_deterministic <- function(fixture) {
    bfpwr_sim_bind_failures(
        bfpwr_sim_validate_sequential_cumulative_summary(
            fixture$cumulative_summary),
        bfpwr_sim_validate_sequential_final_summary(
            fixture$final_summary, fixture$cumulative_summary),
        bfpwr_sim_validate_sequential_search_summary(
            fixture$search_summary, fixture$cumulative_summary)
    )
}

bfpwr_sim_find_fixture_schedule <- function(fixture, schedule_id) {
    hits <- vapply(fixture$spec$schedules, function(x) {
        identical(x$schedule_id, schedule_id)
    }, logical(1))
    if (sum(hits) != 1L) {
        stop("schedule not found or not unique: ", schedule_id)
    }
    fixture$spec$schedules[[which(hits)]]
}

bfpwr_sim_z_sequential_reference <- function(fixture,
                                             row,
                                             bf_prior,
                                             design,
                                             strict = TRUE) {
    schedule <- bfpwr_sim_find_fixture_schedule(
        fixture, row$schedule_id[[1]])
    n <- schedule$n
    se <- design$generation$usd / sqrt(n)
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    prior <- bf_prior$analysis_prior

    if (identical(bf_prior$bf_type, "normal")) {
        pm <- prior$pm - prior$null
        psd <- prior$psd
        type <- "normal"
    } else if (identical(bf_prior$bf_type, "moment")) {
        pm <- NULL
        psd <- prior$psd
        type <- "moment"
    } else if (identical(bf_prior$bf_type, "directional")) {
        pm <- prior$pm - prior$null
        psd <- prior$psd
        type <- "directional"
    } else {
        stop("unsupported z sequential BF type: ", bf_prior$bf_type)
    }

    pbf01seq(
        k1 = row$k1[[1]],
        k0 = row$k0[[1]],
        se = se,
        n = n,
        pm = pm,
        psd = psd,
        dpm = dp$dpm - prior$null,
        dpsd = dp$dpsd,
        type = type,
        strict = strict)
}

bfpwr_sim_z_sequential_reference_case_grid <- function(fixture,
                                                       profile = c("curated",
                                                                   "none")) {
    profile <- match.arg(profile)
    if (profile == "none") {
        return(data.frame(bf_prior_id = character(),
                          schedule_id = character(),
                          evidence_threshold = numeric(),
                          stringsAsFactors = FALSE))
    }

    fixture_set_id <- fixture$spec$fixture_set_id
    base <- switch(
        fixture_set_id,
        "z-bf01-sequential-core-v1" = data.frame(
            bf_prior_id = c(
                "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-short",
                "z-bf01-normal-null0-pm0p2-psd0p3-on-z-dnorm-0p2-s0p5-usd1-short",
                "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-short",
                "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0-usd1-long",
                "z-bf01-point-null0-pm0p2-psd0-on-z-dpoint-0p2-usd1-long",
                "z-bf01-normal-null0p2-pm0p2-psd0p70710678-on-z-dnorm-0p2-s0p5-usd1-short"
            ),
            schedule_id = c(
                "start20-by10-looks02",
                "start20-by10-looks02",
                "start10-by10-looks50",
                "start10-by10-looks100",
                "start10-by10-looks200",
                "start20-by10-looks02"
            ),
            thresholds = c("all", "all", "10", "10", "10", "all"),
            stringsAsFactors = FALSE),
        "z-nmbf01-sequential-core-v1" = data.frame(
            bf_prior_id = c(
                "z-nmbf01-moment-null0-psd0p35355339-on-z-dpoint-0-usd1-short",
                "z-nmbf01-moment-null0-psd0p35355339-on-z-dpoint-0p2-usd1-short",
                "z-nmbf01-moment-null0-psd0p70710678-on-z-dpoint-0p2-usd1-short",
                "z-nmbf01-moment-null0p2-psd0p35355339-on-z-dnorm-0p2-s0p5-usd1-short"
            ),
            schedule_id = c(
                "start20-by10-looks02",
                "start20-by10-looks05",
                "start20-by10-looks05",
                "start20-by10-looks05"
            ),
            thresholds = c("all", "all", "all", "all"),
            stringsAsFactors = FALSE),
        "z-dirbf01-sequential-core-v1" = data.frame(
            bf_prior_id = c(
                "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short",
                "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0p2-usd1-short",
                "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0-usd1-short",
                "z-dirbf01-directional-normal-null0-pm0-psd0p70710678-on-z-dpoint-0p2-usd1-long",
                "z-dirbf01-directional-normal-null0p2-pm0p2-psd0p70710678-on-z-dnorm-0p2-s0p5-usd1-long",
                "z-dirbf01-directional-normal-null0p2-pm0p2-psd0p70710678-on-z-dnorm-0p2-s0p5-usd1-short"
            ),
            schedule_id = c(
                "start20-by10-looks02",
                "start20-by10-looks20",
                "start10-by10-looks50",
                "start10-by10-looks100",
                "start10-by10-looks200",
                "start20-by10-looks20"
            ),
            thresholds = c("all", "all", "10", "10", "10", "all"),
            stringsAsFactors = FALSE),
        data.frame(bf_prior_id = character(),
                   schedule_id = character(),
                   stringsAsFactors = FALSE))

    if (nrow(base) == 0) {
        return(data.frame(bf_prior_id = character(),
                          schedule_id = character(),
                          evidence_threshold = numeric(),
                          stringsAsFactors = FALSE))
    }

    cases <- do.call(rbind, lapply(seq_len(nrow(base)), function(i) {
        thresholds <- if (identical(base$thresholds[[i]], "all")) {
            fixture$spec$evidence_thresholds
        } else {
            as.numeric(strsplit(base$thresholds[[i]], ",",
                                fixed = TRUE)[[1]])
        }
        data.frame(
            bf_prior_id = base$bf_prior_id[[i]],
            schedule_id = base$schedule_id[[i]],
            evidence_threshold = thresholds,
            stringsAsFactors = FALSE)
    }))
    cases <- unique(cases)
    rownames(cases) <- NULL

    available <- unique(fixture$final_summary[
        c("bf_prior_id", "schedule_id", "evidence_threshold")])
    key <- bfpwr_sim_row_key(cases, c("bf_prior_id", "schedule_id",
                                      "evidence_threshold"))
    available_key <- bfpwr_sim_row_key(
        available, c("bf_prior_id", "schedule_id", "evidence_threshold"))
    cases[key %in% available_key, , drop = FALSE]
}

bfpwr_sim_z_sequential_reference_table <- function(fixture,
                                                   designs,
                                                   bf_priors,
                                                   profile = c("curated",
                                                               "none"),
                                                   strict = TRUE,
                                                   cases = NULL) {
    profile <- match.arg(profile)
    if (is.null(cases)) {
        cases <- bfpwr_sim_z_sequential_reference_case_grid(
            fixture, profile = profile)
    }
    if (nrow(cases) == 0) {
        return(list(reference_rows = data.frame(),
                    en_rows = data.frame(),
                    timings = data.frame()))
    }

    out <- list()
    en_out <- list()
    timing_out <- list()
    idx <- 0L
    en_idx <- 0L
    timing_idx <- 0L
    for (i in seq_len(nrow(cases))) {
        case <- cases[i, , drop = FALSE]
        validation_case_id <- if ("package_verification_case_id" %in%
                                  names(case)) {
            case$package_verification_case_id[[1]]
        } else NA_character_
        validation_role <- if ("validation_role" %in% names(case)) {
            case$validation_role[[1]]
        } else "sequential_package_reference"
        include_probability <- !identical(
            validation_role, "sequential_expected_sample_size_reference")
        final_row <- fixture$final_summary[
            fixture$final_summary$bf_prior_id == case$bf_prior_id[[1]] &
                fixture$final_summary$schedule_id == case$schedule_id[[1]] &
                fixture$final_summary$evidence_threshold ==
                    case$evidence_threshold[[1]],
            , drop = FALSE]
        curve <- fixture$cumulative_summary[
            fixture$cumulative_summary$bf_prior_id == case$bf_prior_id[[1]] &
                fixture$cumulative_summary$schedule_id ==
                    case$schedule_id[[1]] &
                fixture$cumulative_summary$evidence_threshold ==
                    case$evidence_threshold[[1]],
            , drop = FALSE]
        curve <- curve[order(curve$look), , drop = FALSE]

        if (nrow(final_row) != 1L || nrow(curve) != final_row$n_looks[[1]]) {
            next
        }

        bf_prior <- bfpwr_sim_find_bf_prior_case(case$bf_prior_id[[1]],
                                                 bf_priors)
        design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                             designs)
        timing_error <- ""
        t0 <- proc.time()[["elapsed"]]
        reference <- tryCatch({
            bfpwr_sim_z_sequential_reference(
                fixture = fixture,
                row = final_row,
                bf_prior = bf_prior,
                design = design,
                strict = strict)
        }, error = function(e) {
            timing_error <<- conditionMessage(e)
            warning("could not compute sequential reference for ",
                    case$bf_prior_id[[1]], " / ",
                    case$schedule_id[[1]], " / BF ",
                    case$evidence_threshold[[1]], ": ",
                    conditionMessage(e))
            NULL
        })
        timing_idx <- timing_idx + 1L
        timing_out[[timing_idx]] <- bfpwr_sim_reference_timing_row(
            rows = final_row,
            package_function = "pbf01seq",
            mode = "sequential",
            elapsed_seconds = proc.time()[["elapsed"]] - t0,
            n_rows = 1L,
            status = if (nzchar(timing_error)) "error" else "ok",
            error = timing_error,
            package_verification_case_id = validation_case_id,
            validation_role = validation_role)

        if (include_probability) {
            tails <- list(
                H1 = list(prob = curve$cum_pH1,
                          reference = if (is.null(reference)) {
                              rep(NA_real_, nrow(curve))
                          } else reference$cumpH1),
                H0 = list(prob = curve$cum_pH0,
                          reference = if (is.null(reference)) {
                              rep(NA_real_, nrow(curve))
                          } else reference$cumpH0),
                Inc = list(prob = curve$cum_pInc,
                           reference = if (is.null(reference)) {
                               rep(NA_real_, nrow(curve))
                           } else reference$cumpInc)
            )
            for (tail in names(tails)) {
                idx <- idx + 1L
                prob <- tails[[tail]]$prob
                out[[idx]] <- data.frame(
                    package_verification_case_id = validation_case_id,
                    validation_role = validation_role,
                    fixture_set_id = curve$fixture_set_id,
                    family = curve$family,
                    bf_type = curve$bf_type,
                    bf_prior_id = curve$bf_prior_id,
                    design_case_id = curve$design_case_id,
                    look_grid_name = curve$look_grid_name,
                    schedule_id = curve$schedule_id,
                    schedule_family_id = curve$schedule_family_id,
                    threshold_pair_id = curve$threshold_pair_id,
                    evidence_threshold = curve$evidence_threshold,
                    tail = tail,
                    look = curve$look,
                    n = curve$n,
                    nsim = curve$nsim,
                    n_event = as.integer(round(prob * curve$nsim)),
                    prob = prob,
                    reference_prob = tails[[tail]]$reference,
                    reference_status =
                        ifelse(is.finite(tails[[tail]]$reference),
                               "available", "unavailable"),
                    stringsAsFactors = FALSE)
            }
        }

        if (!is.null(reference)) {
            en_idx <- en_idx + 1L
            se_en <- sqrt(max(final_row$VarN[[1]], 0) / final_row$nsim[[1]])
            en_out[[en_idx]] <- data.frame(
                package_verification_case_id = validation_case_id,
                validation_role = validation_role,
                fixture_set_id = final_row$fixture_set_id[[1]],
                family = final_row$family[[1]],
                bf_type = final_row$bf_type[[1]],
                bf_prior_id = final_row$bf_prior_id[[1]],
                design_case_id = final_row$design_case_id[[1]],
                look_grid_name = final_row$look_grid_name[[1]],
                schedule_id = final_row$schedule_id[[1]],
                schedule_family_id = final_row$schedule_family_id[[1]],
                threshold_pair_id = final_row$threshold_pair_id[[1]],
                evidence_threshold = final_row$evidence_threshold[[1]],
                n_looks = final_row$n_looks[[1]],
                nsim = final_row$nsim[[1]],
                EN = final_row$EN[[1]],
                reference_EN = reference$EN,
                abs_error = abs(final_row$EN[[1]] - reference$EN),
                tolerance = max(1, 4 * se_en),
                stringsAsFactors = FALSE)
        }
    }

    reference_rows <- if (length(out) == 0) data.frame() else do.call(rbind, out)
    en_rows <- if (length(en_out) == 0) data.frame() else do.call(rbind, en_out)
    timings <- if (length(timing_out) == 0) {
        data.frame()
    } else do.call(rbind, timing_out)
    list(reference_rows = reference_rows,
         en_rows = en_rows,
         timings = timings)
}

bfpwr_sim_t_n2 <- function(n1, design) {
    if (identical(design$generation$type, "two.sample")) {
        pmax(2, round(n1 * design$generation$n2_multiplier))
    } else {
        n1
    }
}

bfpwr_sim_t_sequential_reference <- function(fixture,
                                             row,
                                             bf_prior,
                                             design,
                                             strict = TRUE) {
    schedule <- bfpwr_sim_find_fixture_schedule(
        fixture, row$schedule_id[[1]])
    n1 <- schedule$n
    n2 <- bfpwr_sim_t_n2(n1, design)
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)

    ptbf01seq(
        k1 = row$k1[[1]],
        k0 = row$k0[[1]],
        n = n1,
        n1 = n1,
        n2 = n2,
        plocation = prior$plocation - prior$null,
        pscale = prior$pscale,
        pdf = prior$pdf,
        dpm = dp$dpm - prior$null,
        dpsd = dp$dpsd,
        type = prior$type,
        alternative = prior$alternative,
        strict = strict)
}

bfpwr_sim_t_sequential_reference_case_grid <- function(fixture,
                                                       profile = c("curated",
                                                                   "none")) {
    profile <- match.arg(profile)
    if (profile == "none") {
        return(data.frame(bf_prior_id = character(),
                          schedule_id = character(),
                          evidence_threshold = numeric(),
                          stringsAsFactors = FALSE))
    }

    if (!identical(fixture$spec$fixture_set_id,
                   "t-tbf01-sequential-core-v1")) {
        return(data.frame(bf_prior_id = character(),
                          schedule_id = character(),
                          evidence_threshold = numeric(),
                          stringsAsFactors = FALSE))
    }

    cases <- data.frame(
        bf_prior_id = c(
            "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-0p5-short",
            "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-less-on-t-two-dpoint-0p5-short",
            "t-tbf01-student-t-null0-loc0p5-scale0p1-df3-greater-on-t-one-dpoint-0p5-short",
            "t-tbf01-cauchy-null0-loc0-scale0p70710678-df1-greater-on-t-two-dpoint-m0p4-short"
        ),
        schedule_id = rep("start20-by10-looks20", 4),
        evidence_threshold = c(10, 10, 10, 3),
        stringsAsFactors = FALSE)

    available <- unique(fixture$final_summary[
        c("bf_prior_id", "schedule_id", "evidence_threshold")])
    key <- bfpwr_sim_row_key(cases, c("bf_prior_id", "schedule_id",
                                      "evidence_threshold"))
    available_key <- bfpwr_sim_row_key(
        available, c("bf_prior_id", "schedule_id", "evidence_threshold"))
    cases[key %in% available_key, , drop = FALSE]
}

bfpwr_sim_t_sequential_reference_table <- function(fixture,
                                                   designs,
                                                   bf_priors,
                                                   profile = c("curated",
                                                               "none"),
                                                   strict = TRUE,
                                                   cases = NULL) {
    profile <- match.arg(profile)
    if (is.null(cases)) {
        cases <- bfpwr_sim_t_sequential_reference_case_grid(
            fixture, profile = profile)
    }
    if (nrow(cases) == 0) {
        return(list(reference_rows = data.frame(),
                    en_rows = data.frame(),
                    timings = data.frame()))
    }

    out <- list()
    en_out <- list()
    timing_out <- list()
    idx <- 0L
    en_idx <- 0L
    timing_idx <- 0L
    for (i in seq_len(nrow(cases))) {
        case <- cases[i, , drop = FALSE]
        validation_case_id <- if ("package_verification_case_id" %in%
                                  names(case)) {
            case$package_verification_case_id[[1]]
        } else NA_character_
        validation_role <- if ("validation_role" %in% names(case)) {
            case$validation_role[[1]]
        } else "sequential_package_reference"
        include_probability <- !identical(
            validation_role, "sequential_expected_sample_size_reference")
        final_row <- fixture$final_summary[
            fixture$final_summary$bf_prior_id == case$bf_prior_id[[1]] &
                fixture$final_summary$schedule_id == case$schedule_id[[1]] &
                fixture$final_summary$evidence_threshold ==
                    case$evidence_threshold[[1]],
            , drop = FALSE]
        curve <- fixture$cumulative_summary[
            fixture$cumulative_summary$bf_prior_id == case$bf_prior_id[[1]] &
                fixture$cumulative_summary$schedule_id ==
                    case$schedule_id[[1]] &
                fixture$cumulative_summary$evidence_threshold ==
                    case$evidence_threshold[[1]],
            , drop = FALSE]
        curve <- curve[order(curve$look), , drop = FALSE]

        if (nrow(final_row) != 1L || nrow(curve) != final_row$n_looks[[1]]) {
            next
        }

        bf_prior <- bfpwr_sim_find_bf_prior_case(case$bf_prior_id[[1]],
                                                 bf_priors)
        design <- bfpwr_sim_find_design_case(final_row$design_case_id[[1]],
                                             designs)
        timing_error <- ""
        t0 <- proc.time()[["elapsed"]]
        reference <- tryCatch({
            bfpwr_sim_t_sequential_reference(
                fixture = fixture,
                row = final_row,
                bf_prior = bf_prior,
                design = design,
                strict = strict)
        }, error = function(e) {
            timing_error <<- conditionMessage(e)
            warning("could not compute t sequential reference for ",
                    case$bf_prior_id[[1]], " / ",
                    case$schedule_id[[1]], " / BF ",
                    case$evidence_threshold[[1]], ": ",
                    conditionMessage(e))
            NULL
        })
        timing_idx <- timing_idx + 1L
        timing_out[[timing_idx]] <- bfpwr_sim_reference_timing_row(
            rows = final_row,
            package_function = "ptbf01seq",
            mode = "sequential",
            elapsed_seconds = proc.time()[["elapsed"]] - t0,
            n_rows = 1L,
            status = if (nzchar(timing_error)) "error" else "ok",
            error = timing_error,
            package_verification_case_id = validation_case_id,
            validation_role = validation_role)

        if (include_probability) {
            tails <- list(
                H1 = list(prob = curve$cum_pH1,
                          reference = if (is.null(reference)) {
                              rep(NA_real_, nrow(curve))
                          } else reference$cumpH1),
                H0 = list(prob = curve$cum_pH0,
                          reference = if (is.null(reference)) {
                              rep(NA_real_, nrow(curve))
                          } else reference$cumpH0),
                Inc = list(prob = curve$cum_pInc,
                           reference = if (is.null(reference)) {
                               rep(NA_real_, nrow(curve))
                           } else reference$cumpInc)
            )
            for (tail in names(tails)) {
                idx <- idx + 1L
                prob <- tails[[tail]]$prob
                out[[idx]] <- data.frame(
                    package_verification_case_id = validation_case_id,
                    validation_role = validation_role,
                    fixture_set_id = curve$fixture_set_id,
                    family = curve$family,
                    bf_type = curve$bf_type,
                    bf_prior_id = curve$bf_prior_id,
                    design_case_id = curve$design_case_id,
                    look_grid_name = curve$look_grid_name,
                    schedule_id = curve$schedule_id,
                    schedule_family_id = curve$schedule_family_id,
                    threshold_pair_id = curve$threshold_pair_id,
                    evidence_threshold = curve$evidence_threshold,
                    tail = tail,
                    look = curve$look,
                    n = curve$n,
                    nsim = curve$nsim,
                    n_event = as.integer(round(prob * curve$nsim)),
                    prob = prob,
                    reference_prob = tails[[tail]]$reference,
                    reference_status =
                        ifelse(is.finite(tails[[tail]]$reference),
                               "available", "unavailable"),
                    stringsAsFactors = FALSE)
            }
        }

        if (!is.null(reference)) {
            en_idx <- en_idx + 1L
            se_en <- sqrt(max(final_row$VarN[[1]], 0) / final_row$nsim[[1]])
            ## Sequential t EN is a mean over a discrete look schedule. The
            ## expanded package manifest includes high-threshold cases where
            ## Monte Carlo and integration agree within about one grid step but
            ## can exceed a one-sample absolute floor. Keep the tolerance tied
            ## to both MCSE and the schedule resolution.
            grid_floor <- 0.15 * final_row$increment[[1]]
            en_out[[en_idx]] <- data.frame(
                package_verification_case_id = validation_case_id,
                validation_role = validation_role,
                fixture_set_id = final_row$fixture_set_id[[1]],
                family = final_row$family[[1]],
                bf_type = final_row$bf_type[[1]],
                bf_prior_id = final_row$bf_prior_id[[1]],
                design_case_id = final_row$design_case_id[[1]],
                look_grid_name = final_row$look_grid_name[[1]],
                schedule_id = final_row$schedule_id[[1]],
                schedule_family_id = final_row$schedule_family_id[[1]],
                threshold_pair_id = final_row$threshold_pair_id[[1]],
                evidence_threshold = final_row$evidence_threshold[[1]],
                n_looks = final_row$n_looks[[1]],
                nsim = final_row$nsim[[1]],
                EN = final_row$EN[[1]],
                reference_EN = reference$EN1,
                abs_error = abs(final_row$EN[[1]] - reference$EN1),
                tolerance = max(1, 4 * se_en, grid_floor),
                stringsAsFactors = FALSE)
        }
    }

    reference_rows <- if (length(out) == 0) data.frame() else do.call(rbind, out)
    en_rows <- if (length(en_out) == 0) data.frame() else do.call(rbind, en_out)
    timings <- if (length(timing_out) == 0) {
        data.frame()
    } else do.call(rbind, timing_out)
    list(reference_rows = reference_rows,
         en_rows = en_rows,
         timings = timings)
}

bfpwr_sim_z_directional_reference <- function(k1,
                                              k0,
                                              n,
                                              usd,
                                              prior,
                                              dpm,
                                              dpsd,
                                              tail = c("H1", "H0"),
                                              strict = TRUE) {
    tail <- match.arg(tail)
    ## pbf01seq() encodes directional z references with the boundary at zero.
    ## Shift analysis and design means so nonzero dirbf01() nulls agree.
    res <- pbf01seq(
        k1 = k1,
        k0 = k0,
        se = usd / sqrt(n),
        n = n,
        pm = prior$pm - prior$null,
        psd = prior$psd,
        dpm = dpm - prior$null,
        dpsd = dpsd,
        type = "directional",
        strict = strict)
    if (tail == "H1") {
        tail(res$cumpH1, 1)
    } else {
        tail(res$cumpH0, 1)
    }
}

bfpwr_sim_fixed_reference_function_name <- function(bf_prior) {
    if (identical(bf_prior$test_family, "z")) {
        if (identical(bf_prior$bf_type, "moment")) return("pnmbf01")
        if (identical(bf_prior$bf_type, "directional")) return("pbf01seq")
        return("pbf01")
    }
    if (identical(bf_prior$test_family, "binomial")) return("pbinbf01")
    if (identical(bf_prior$test_family, "t")) return("ptbf01")
    NA_character_
}

bfpwr_sim_reference_timing_row <- function(rows,
                                           package_function,
                                           mode,
                                           elapsed_seconds,
                                           n_rows,
                                           status = "ok",
                                           error = "",
                                           package_verification_case_id =
                                               NA_character_,
                                           validation_role = NA_character_) {
    first <- if (nrow(rows) > 0) rows[1, , drop = FALSE] else data.frame()
    value <- function(name, default = NA) {
        if (name %in% names(first)) first[[name]][[1]] else default
    }
    if (is.na(validation_role) && "validation_role" %in% names(first)) {
        validation_role <- value("validation_role", NA_character_)
    }
    data.frame(
        package_verification_case_id = package_verification_case_id,
        validation_role = validation_role,
        fixture_set_id = value("fixture_set_id", NA_character_),
        family = value("family", NA_character_),
        mode = mode,
        bf_type = value("bf_type", NA_character_),
        package_function = package_function,
        bf_prior_id = value("bf_prior_id", NA_character_),
        design_case_id = value("design_case_id", NA_character_),
        look_grid_name = value("look_grid_name", NA_character_),
        schedule_id = value("schedule_id", NA_character_),
        evidence_threshold = value("evidence_threshold", NA_real_),
        threshold_id = value("threshold_id", NA_character_),
        tail = value("tail", NA_character_),
        n_rows = as.integer(n_rows),
        n_looks = value("n_looks", NA_integer_),
        elapsed_seconds = as.numeric(elapsed_seconds),
        timing_status = status,
        timing_error = error,
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_fixed_reference_group <- function(rows, bf_prior, design) {
    prior <- bf_prior$analysis_prior
    lower_tail <- rows$tail == "H1"

    if (identical(bf_prior$test_family, "z")) {
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        if (identical(bf_prior$bf_type, "normal")) {
            return(pbf01(
                k = rows$threshold,
                n = rows$n,
                usd = design$generation$usd,
                null = prior$null,
                pm = prior$pm,
                psd = prior$psd,
                dpm = dp$dpm,
                dpsd = dp$dpsd,
                lower.tail = lower_tail))
        }
        if (identical(bf_prior$bf_type, "moment")) {
            return(pnmbf01(
                k = rows$threshold,
                n = rows$n,
                usd = design$generation$usd,
                null = prior$null,
                psd = prior$psd,
                dpm = dp$dpm,
                dpsd = dp$dpsd,
                lower.tail = lower_tail))
        }
        if (identical(bf_prior$bf_type, "directional")) {
            return(vapply(seq_len(nrow(rows)), function(i) {
                bfpwr_sim_z_directional_reference(
                    k1 = 1 / rows$evidence_threshold[[i]],
                    k0 = rows$evidence_threshold[[i]],
                    n = rows$n[[i]],
                    usd = design$generation$usd,
                    prior = prior,
                    dpm = dp$dpm,
                    dpsd = dp$dpsd,
                    tail = rows$tail[[i]])
            }, numeric(1)))
        }
        return(rep(NA_real_, nrow(rows)))
    }

    if (identical(bf_prior$test_family, "binomial")) {
        design_args <- bfpwr_sim_binomial_design_args(design)
        return(pbinbf01(
            k = rows$threshold,
            n = rows$n,
            p0 = prior$p0,
            type = bf_prior$bf_type,
            a = prior$a,
            b = prior$b,
            dp = design_args$dp,
            da = design_args$da,
            db = design_args$db,
            dl = design_args$dl,
            du = design_args$du,
            lower.tail = lower_tail))
    }

    if (identical(bf_prior$test_family, "t")) {
        dp <- bfpwr_sim_design_prior_mean_sd(design)
        n1 <- rows$n
        n2 <- if (identical(design$generation$type, "two.sample")) {
            pmax(2, round(n1 * design$generation$n2_multiplier))
        } else {
            n1
        }
        return(ptbf01(
            k = rows$threshold,
            n = n1,
            n1 = n1,
            n2 = n2,
            null = prior$null,
            plocation = prior$plocation,
            pscale = prior$pscale,
            pdf = prior$pdf,
            dpm = dp$dpm,
            dpsd = dp$dpsd,
            type = prior$type,
            alternative = prior$alternative,
            lower.tail = lower_tail))
    }

    rep(NA_real_, nrow(rows))
}

bfpwr_sim_fixed_reference_table <- function(tail_summary,
                                            designs,
                                            bf_priors,
                                            manifest_cases = NULL) {
    if (!is.null(manifest_cases)) {
        keep <- manifest_cases$validation_role ==
            "fixed_probability_reference"
        manifest_cases <- manifest_cases[keep, , drop = FALSE]
        if (nrow(manifest_cases) == 0) {
            tail_summary <- tail_summary[0, , drop = FALSE]
        } else {
            exact_cols <- c("fixture_set_id", "bf_prior_id", "threshold_id",
                            "tail", "n")
            has_exact_cells <- all(exact_cols %in% names(manifest_cases)) &&
                all(exact_cols %in% names(tail_summary)) &&
                any(nzchar(as.character(manifest_cases$threshold_id))) &&
                any(nzchar(as.character(manifest_cases$tail))) &&
                any(is.finite(suppressWarnings(as.numeric(manifest_cases$n))))
            if (has_exact_cells) {
                case_key <- bfpwr_sim_row_key(manifest_cases, exact_cols)
                if (any(duplicated(case_key))) {
                    stop("fixed package manifest has duplicate exact cells",
                         call. = FALSE)
                }
                tail_key <- bfpwr_sim_row_key(tail_summary, exact_cols)
                idx <- match(tail_key, case_key)
                keep_tail <- !is.na(idx)
                tail_summary <- tail_summary[keep_tail, , drop = FALSE]
                idx <- idx[keep_tail]
                tail_summary$package_verification_case_id <-
                    manifest_cases$package_verification_case_id[idx]
                tail_summary$validation_role <-
                    manifest_cases$validation_role[idx]
            } else if ("bf_prior_id" %in% names(manifest_cases)) {
                tail_summary <- tail_summary[
                    tail_summary$bf_prior_id %in% manifest_cases$bf_prior_id,
                    , drop = FALSE]
            }
        }
    }
    ids <- unique(tail_summary$bf_prior_id)
    out <- vector("list", length(ids))
    timings <- vector("list", length(ids))
    for (i in seq_along(ids)) {
        id <- ids[[i]]
        rows <- tail_summary[tail_summary$bf_prior_id == id, , drop = FALSE]
        bf_prior <- bfpwr_sim_find_bf_prior_case(id, bf_priors)
        design <- bfpwr_sim_find_design_case(rows$design_case_id[[1]],
                                             designs)
        timing_error <- ""
        t0 <- proc.time()[["elapsed"]]
        reference <- tryCatch({
            bfpwr_sim_fixed_reference_group(rows, bf_prior, design)
        }, error = function(e) {
            timing_error <<- conditionMessage(e)
            warning("could not compute fixed reference for ", id, ": ",
                    conditionMessage(e))
            rep(NA_real_, nrow(rows))
        })
        elapsed <- proc.time()[["elapsed"]] - t0
        timings[[i]] <- bfpwr_sim_reference_timing_row(
            rows = rows,
            package_function = bfpwr_sim_fixed_reference_function_name(
                bf_prior),
            mode = "fixed",
            elapsed_seconds = elapsed,
            n_rows = nrow(rows),
            status = if (nzchar(timing_error)) "error" else "ok",
            error = timing_error)
        rows$reference_prob <- reference
        rows$reference_status <- ifelse(is.finite(reference),
                                        "available", "unavailable")
        if (!is.null(manifest_cases) && nrow(manifest_cases) > 0 &&
            !"package_verification_case_id" %in% names(rows)) {
            case_row <- manifest_cases[
                manifest_cases$bf_prior_id == id, , drop = FALSE]
            if (nrow(case_row) > 0) {
                rows$package_verification_case_id <-
                    case_row$package_verification_case_id[[1]]
                rows$validation_role <- case_row$validation_role[[1]]
            }
        }
        out[[i]] <- rows
    }
    if (length(out) == 0) {
        empty <- tail_summary[0, , drop = FALSE]
        attr(empty, "timings") <- data.frame()
        return(empty)
    }
    reference_rows <- do.call(rbind, out)
    attr(reference_rows, "timings") <- do.call(rbind, timings)
    reference_rows
}

bfpwr_sim_hash_uniform <- function(key) {
    vapply(as.character(key), function(x) {
        bytes <- as.integer(charToRaw(enc2utf8(x)))
        hash <- 1
        for (byte in bytes) {
            hash <- (hash * 131 + byte + 1) %% 2147483647
        }
        u <- sin(hash * 12.9898 + 78.233) * 43758.5453123
        u - floor(u)
    }, numeric(1))
}

bfpwr_sim_randomized_quantile_residual <- function(event,
                                                   probability,
                                                   key,
                                                   eps = 1e-12) {
    event <- as.logical(event)
    probability <- as.numeric(probability)
    valid <- !is.na(event) & is.finite(probability) &
        probability > 0 & probability < 1
    out <- rep(NA_real_, length(probability))
    if (!any(valid)) {
        return(out)
    }

    p <- probability[valid]
    lower <- ifelse(event[valid], 1 - p, 0)
    upper <- ifelse(event[valid], 1, 1 - p)
    u <- lower + bfpwr_sim_hash_uniform(key[valid]) * (upper - lower)
    u <- pmin(1 - eps, pmax(eps, u))
    out[valid] <- stats::qnorm(u)
    out
}

bfpwr_sim_fixed_cell_metadata <- function(rows, designs, bf_priors) {
    value <- function(x, name, default = NA_real_) {
        if (!is.null(x[[name]])) x[[name]] else default
    }
    design_meta <- do.call(rbind, lapply(designs, function(design) {
        data.frame(
            design_case_id = design$design_case_id,
            design_prior_form = design$design_prior$family,
            design_tier = design$tier,
            stringsAsFactors = FALSE
        )
    }))
    prior_meta <- do.call(rbind, lapply(bf_priors, function(prior_case) {
        data.frame(
            bf_prior_id = prior_case$bf_prior_id,
            analysis_prior_form = prior_case$prior_form,
            analysis_null = value(prior_case$analysis_prior, "null"),
            analysis_pm = value(prior_case$analysis_prior, "pm"),
            analysis_psd = value(prior_case$analysis_prior, "psd"),
            stringsAsFactors = FALSE
        )
    }))
    rows <- merge(rows, design_meta, by = "design_case_id", all.x = TRUE,
                  sort = FALSE)
    rows <- merge(rows, prior_meta, by = "bf_prior_id", all.x = TRUE,
                  sort = FALSE)
    rows$look_fraction <- vapply(seq_len(nrow(rows)), function(i) {
        design <- bfpwr_sim_find_design_case(rows$design_case_id[[i]],
                                             designs)
        index <- match(rows$n[[i]], design$look_grid)
        if (is.na(index)) return(NA_real_)
        index / length(design$look_grid)
    }, numeric(1))
    rows
}

bfpwr_sim_select_independent_fixed_cells <- function(reference_rows,
                                                     designs,
                                                     bf_priors,
                                                     max_cells = 12L) {
    rows <- reference_rows[is.finite(reference_rows$reference_prob) &
                               reference_rows$reference_prob > 0 &
                               reference_rows$reference_prob < 1 &
                               reference_rows$nsim > 0, , drop = FALSE]
    if (nrow(rows) == 0) {
        return(data.frame())
    }
    rows <- bfpwr_sim_fixed_cell_metadata(rows, designs, bf_priors)
    rows <- rows[is.finite(rows$look_fraction), , drop = FALSE]
    if (nrow(rows) == 0) {
        return(data.frame())
    }

    design_ids <- sort(unique(rows$design_case_id))
    if (is.finite(max_cells)) {
        design_ids <- head(design_ids, max(0L, as.integer(max_cells)))
    }
    target_tail <- rep(c("H1", "H0"), length.out = length(design_ids))
    target_threshold <- rep(c(3, 10, 30), length.out = length(design_ids))
    target_prior <- rep(c("point", "normal"), length.out = length(design_ids))
    target_look <- rep(c(0.25, 0.50, 0.75), length.out = length(design_ids))

    selected <- vector("list", length(design_ids))
    for (i in seq_along(design_ids)) {
        group <- rows[rows$design_case_id == design_ids[[i]], , drop = FALSE]
        score <- 100 * (group$tail != target_tail[[i]]) +
            30 * (abs(group$evidence_threshold -
                          target_threshold[[i]]) > 1e-12) +
            20 * (group$analysis_prior_form != target_prior[[i]]) +
            abs(group$look_fraction - target_look[[i]]) +
            0.01 * abs(group$reference_prob - 0.5) +
            ifelse(group$reference_prob < 0.01 |
                       group$reference_prob > 0.99, 25,
                   ifelse(group$reference_prob < 0.05 |
                              group$reference_prob > 0.95, 5, 0))
        selected[[i]] <- group[which.min(score), , drop = FALSE]
    }
    out <- do.call(rbind, selected)
    rownames(out) <- NULL
    out$diagnostic_cell_id <- sprintf("cell-%02d", seq_len(nrow(out)))
    out
}

bfpwr_sim_fixed_event_probabilities <- function(cell,
                                                bf_prior,
                                                design,
                                                trajectories,
                                                reference = c("conditional",
                                                              "unconditional")) {
    reference <- match.arg(reference)
    if (identical(reference, "unconditional")) {
        return(rep(cell$reference_prob[[1]], nrow(trajectories)))
    }

    prior <- bf_prior$analysis_prior
    lower_tail <- cell$tail[[1]] == "H1"
    if (identical(bf_prior$test_family, "z")) {
        if (identical(bf_prior$bf_type, "normal")) {
            return(bfpwr_sim_package_function("pbf01")(
                k = cell$threshold[[1]],
                n = trajectories$n,
                usd = design$generation$usd,
                null = prior$null,
                pm = prior$pm,
                psd = prior$psd,
                dpm = trajectories$true_effect,
                dpsd = 0,
                lower.tail = lower_tail))
        }
        if (identical(bf_prior$bf_type, "moment")) {
            return(bfpwr_sim_package_function("pnmbf01")(
                k = cell$threshold[[1]],
                n = trajectories$n,
                usd = design$generation$usd,
                null = prior$null,
                psd = prior$psd,
                dpm = trajectories$true_effect,
                dpsd = 0,
                lower.tail = lower_tail))
        }
        if (identical(bf_prior$bf_type, "directional")) {
            return(vapply(seq_len(nrow(trajectories)), function(i) {
                bfpwr_sim_z_directional_reference(
                    k1 = 1 / cell$evidence_threshold[[1]],
                    k0 = cell$evidence_threshold[[1]],
                    n = trajectories$n[[i]],
                    usd = design$generation$usd,
                    prior = prior,
                    dpm = trajectories$true_effect[[i]],
                    dpsd = 0,
                    tail = cell$tail[[1]])
            }, numeric(1)))
        }
    }

    if (identical(bf_prior$test_family, "t") &&
        "true_effect" %in% names(trajectories)) {
        return(bfpwr_sim_package_function("ptbf01")(
            k = cell$threshold[[1]],
            n = trajectories$n,
            n1 = trajectories$n1,
            n2 = trajectories$n2,
            null = prior$null,
            plocation = prior$plocation,
            pscale = prior$pscale,
            pdf = prior$pdf,
            dpm = trajectories$true_effect,
            dpsd = 0,
            type = prior$type,
            alternative = prior$alternative,
            lower.tail = lower_tail))
    }

    if (identical(bf_prior$test_family, "binomial") &&
        "true_p" %in% names(trajectories)) {
        return(bfpwr_sim_package_function("pbinbf01")(
            k = cell$threshold[[1]],
            n = trajectories$n,
            p0 = prior$p0,
            type = bf_prior$bf_type,
            a = prior$a,
            b = prior$b,
            dp = trajectories$true_p,
            lower.tail = lower_tail))
    }

    rep(cell$reference_prob[[1]], nrow(trajectories))
}

bfpwr_sim_fixed_independent_residuals <- function(corpus_root,
                                                  fixture,
                                                  designs,
                                                  bf_priors,
                                                  reference_rows = NULL,
                                                  reference = c("conditional",
                                                                "unconditional"),
                                                  max_cells = 12L,
                                                  max_chunks_per_cell = Inf) {
    reference <- match.arg(reference)
    if (is.null(reference_rows)) {
        reference_rows <- bfpwr_sim_fixed_reference_table(
            fixture$tail_summary, designs, bf_priors)
    }
    cells <- bfpwr_sim_select_independent_fixed_cells(
        reference_rows, designs, bf_priors, max_cells = max_cells)
    if (nrow(cells) == 0) {
        return(list(
            residual_rows = data.frame(),
            selected_cells = cells,
            summary = data.frame(
                reference = reference,
                n_cells = 0L,
                n_residuals = 0L,
                mean = NA_real_,
                sd = NA_real_,
                p_abs_gt_2 = NA_real_,
                stringsAsFactors = FALSE)))
    }

    residuals <- list()
    row_id <- 0L
    for (i in seq_len(nrow(cells))) {
        cell <- cells[i, , drop = FALSE]
        design <- bfpwr_sim_find_design_case(cell$design_case_id[[1]],
                                             designs)
        bf_prior <- bfpwr_sim_find_bf_prior_case(cell$bf_prior_id[[1]],
                                                 bf_priors)
        bf_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior, design)
        design_files <- bfpwr_sim_design_chunk_files(corpus_root, design)
        chunk_indices <- seq_along(bf_files)
        if (is.finite(max_chunks_per_cell)) {
            max_chunks <- max(0L, as.integer(max_chunks_per_cell))
            if (max_chunks < length(chunk_indices)) {
                chunk_indices <- unique(round(seq(1, length(chunk_indices),
                                                  length.out = max_chunks)))
            }
        }
        for (chunk_i in chunk_indices) {
            bf_rows <- readRDS(bf_files[[chunk_i]])
            bf_rows <- bf_rows[bf_rows$n == cell$n[[1]], , drop = FALSE]
            if (nrow(bf_rows) == 0) next

            trajectories <- readRDS(design_files[[chunk_i]])
            trajectories <- trajectories[trajectories$n == cell$n[[1]],
                                         , drop = FALSE]
            trajectories <- trajectories[
                match(bf_rows$replicate_id, trajectories$replicate_id),
                , drop = FALSE]
            if (any(is.na(trajectories$replicate_id))) {
                stop("could not align design and BF rows for ",
                     cell$bf_prior_id[[1]], " chunk ", chunk_i)
            }

            hit <- if (cell$tail[[1]] == "H1") {
                bf_rows$log_bf01 <= cell$log_threshold[[1]]
            } else {
                bf_rows$log_bf01 >= cell$log_threshold[[1]]
            }
            hit[is.na(hit)] <- FALSE
            event_prob <- bfpwr_sim_fixed_event_probabilities(
                cell, bf_prior, design, trajectories, reference = reference)
            key <- paste(cell$diagnostic_cell_id[[1]],
                         bf_rows$chunk_id,
                         bf_rows$replicate_id,
                         sep = "\r")
            residual <- bfpwr_sim_randomized_quantile_residual(
                hit, event_prob, key)
            keep <- is.finite(residual)
            if (!any(keep)) next

            row_id <- row_id + 1L
            residuals[[row_id]] <- data.frame(
                fixture_set_id = cell$fixture_set_id[[1]],
                diagnostic_cell_id = cell$diagnostic_cell_id[[1]],
                reference = reference,
                family = cell$family[[1]],
                bf_type = cell$bf_type[[1]],
                bf_prior_id = cell$bf_prior_id[[1]],
                design_case_id = cell$design_case_id[[1]],
                design_prior_form = cell$design_prior_form[[1]],
                analysis_prior_form = cell$analysis_prior_form[[1]],
                look_grid_name = cell$look_grid_name[[1]],
                tail = cell$tail[[1]],
                evidence_threshold = cell$evidence_threshold[[1]],
                n = cell$n[[1]],
                chunk_id = bf_rows$chunk_id[keep],
                replicate_id = bf_rows$replicate_id[keep],
                event = as.integer(hit[keep]),
                event_probability = event_prob[keep],
                residual = residual[keep],
                stringsAsFactors = FALSE
            )
        }
    }

    residual_rows <- if (length(residuals) == 0) {
        data.frame()
    } else {
        do.call(rbind, residuals)
    }
    z <- residual_rows$residual
    summary <- data.frame(
        reference = reference,
        n_cells = nrow(cells),
        n_residuals = length(z),
        mean = if (length(z) > 0) mean(z) else NA_real_,
        sd = if (length(z) > 1) stats::sd(z) else NA_real_,
        p_abs_gt_2 = if (length(z) > 0) mean(abs(z) > 2) else NA_real_,
        p_abs_gt_3 = if (length(z) > 0) mean(abs(z) > 3) else NA_real_,
        stringsAsFactors = FALSE
    )
    list(residual_rows = residual_rows,
         selected_cells = cells,
         summary = summary)
}

bfpwr_sim_mc_reference_diagnostics <- function(rows,
                                               label = "fixed_tail",
                                               family_alpha = 0.01,
                                               severe_alpha = 1e-8,
                                               group_cols = c("family",
                                                              "bf_type",
                                                              "tail",
                                                              "evidence_threshold",
                                                              "look_grid_name"),
                                               min_group_rows = 30L) {
    rows <- rows[is.finite(rows$reference_prob) &
                     rows$reference_prob >= 0 &
                     rows$reference_prob <= 1, , drop = FALSE]
    if (nrow(rows) == 0) {
        return(list(
            summary = data.frame(label = label, rows_checked = 0L,
                                 rows_with_z = 0L,
                                 unsupported_rows = NA_integer_,
                                 exact_p_min = NA_real_,
                                 holm_p_min = NA_real_,
                                 holm_failures = NA_integer_,
                                 exact_p_lt_0p05 = NA_integer_,
                                 exact_p_lt_0p01 = NA_integer_,
                                 exact_p_lt_0p001 = NA_integer_,
                                 exact_p_lt_0p0001 = NA_integer_,
                                 mean_z = NA_real_, sd_z = NA_real_,
                                 median_abs_error = NA_real_,
                                 max_abs_error = NA_real_,
                                 coverage_1se = NA_real_,
                                 coverage_2se = NA_real_,
                                 coverage_3se = NA_real_,
                                 outliers_abs_z_gt_4 = NA_integer_,
                                 outliers_abs_z_gt_4_bound = NA_integer_,
                                 outliers_abs_z_gt_5 = NA_integer_,
                                 outliers_abs_z_gt_5_bound = NA_integer_,
                                 severe_interval_failures = NA_integer_,
                                 stringsAsFactors = FALSE),
            group_summary = data.frame(),
            outliers = data.frame(),
            failures = bfpwr_sim_bind_failures()))
    }

    rows$abs_error <- abs(rows$prob - rows$reference_prob)
    ref_mcse <- bfpwr_sim_mcse(rows$reference_prob, rows$nsim)
    rows$reference_mcse <- ref_mcse
    rows$z <- NA_real_
    interior_for_z <- rows$nsim * rows$reference_prob >= 10 &
        rows$nsim * (1 - rows$reference_prob) >= 10
    with_z <- ref_mcse > 0 & interior_for_z
    rows$z[with_z] <- (rows$prob[with_z] -
                           rows$reference_prob[with_z]) / ref_mcse[with_z]

    exact_p <- rep(NA_real_, nrow(rows))
    p0 <- rows$reference_prob == 0
    p1 <- rows$reference_prob == 1
    interior <- !(p0 | p1)
    exact_p[p0] <- ifelse(rows$n_event[p0] == 0, 1, 0)
    exact_p[p1] <- ifelse(rows$n_event[p1] == rows$nsim[p1], 1, 0)
    if (any(interior)) {
        lower <- stats::pbinom(rows$n_event[interior],
                               rows$nsim[interior],
                               rows$reference_prob[interior])
        upper <- stats::pbinom(rows$n_event[interior] - 1,
                               rows$nsim[interior],
                               rows$reference_prob[interior],
                               lower.tail = FALSE)
        exact_p[interior] <- pmin(1, 2 * pmin(lower, upper))
    }
    rows$exact_p <- exact_p
    rows$holm_p <- stats::p.adjust(exact_p, method = "holm")
    holm_bad <- rows$holm_p < family_alpha

    lower <- stats::qbinom(severe_alpha / 2, rows$nsim,
                           rows$reference_prob)
    upper <- stats::qbinom(1 - severe_alpha / 2, rows$nsim,
                           rows$reference_prob)
    severe_bad <- rows$n_event < lower | rows$n_event > upper

    z <- rows$z[is.finite(rows$z)]
    m <- length(z)
    p4 <- 2 * stats::pnorm(-4)
    p5 <- 2 * stats::pnorm(-5)
    out4 <- sum(abs(z) > 4)
    out5 <- sum(abs(z) > 5)
    expected4 <- m * p4
    expected5 <- m * p5

    summary <- data.frame(
        label = label,
        rows_checked = nrow(rows),
        rows_with_z = m,
        unsupported_rows = NA_integer_,
        exact_p_min = min(rows$exact_p, na.rm = TRUE),
        holm_p_min = min(rows$holm_p, na.rm = TRUE),
        holm_failures = sum(holm_bad, na.rm = TRUE),
        exact_p_lt_0p05 = sum(rows$exact_p < 0.05, na.rm = TRUE),
        exact_p_lt_0p01 = sum(rows$exact_p < 0.01, na.rm = TRUE),
        exact_p_lt_0p001 = sum(rows$exact_p < 0.001, na.rm = TRUE),
        exact_p_lt_0p0001 = sum(rows$exact_p < 0.0001, na.rm = TRUE),
        mean_z = if (m > 0) mean(z) else NA_real_,
        sd_z = if (m > 1) stats::sd(z) else NA_real_,
        median_abs_error = stats::median(rows$abs_error),
        max_abs_error = max(rows$abs_error),
        coverage_1se = if (m > 0) mean(abs(z) <= 1) else NA_real_,
        coverage_2se = if (m > 0) mean(abs(z) <= 2) else NA_real_,
        coverage_3se = if (m > 0) mean(abs(z) <= 3) else NA_real_,
        outliers_abs_z_gt_4 = out4,
        expected_abs_z_gt_4 = expected4,
        outliers_abs_z_gt_5 = out5,
        expected_abs_z_gt_5 = expected5,
        severe_interval_failures = sum(severe_bad),
        stringsAsFactors = FALSE
    )

    groups <- data.frame()
    group_cols <- group_cols[group_cols %in% names(rows)]
    if (length(group_cols) > 0 && m > 0) {
        z_rows <- rows[is.finite(rows$z), , drop = FALSE]
        split_rows <- split(z_rows, bfpwr_sim_row_key(z_rows, group_cols))
        groups <- do.call(rbind, lapply(split_rows, function(g) {
            data.frame(
                group = paste(g[1, group_cols], collapse = " | "),
                n_rows = nrow(g),
                mean_z = mean(g$z),
                sd_z = if (nrow(g) > 1) stats::sd(g$z) else NA_real_,
                coverage_2se = mean(abs(g$z) <= 2),
                max_abs_z = max(abs(g$z)),
                stringsAsFactors = FALSE
            )
        }))
        rownames(groups) <- NULL
    }

    outliers <- rows[order(-abs(rows$z)), , drop = FALSE]
    outliers <- outliers[is.finite(outliers$z), , drop = FALSE]
    keep_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                   "design_case_id", "look_grid_name", "schedule_id",
                   "schedule_family_id", "threshold_pair_id", "tail",
                   "evidence_threshold", "look", "n", "nsim", "n_event",
                   "prob", "reference_prob", "reference_mcse", "exact_p",
                   "holm_p", "z", "abs_error")
    outliers <- rows[order(rows$holm_p, -abs(rows$z)), , drop = FALSE]
    outliers <- outliers[seq_len(min(20L, nrow(outliers))),
                         keep_cols[keep_cols %in% names(outliers)],
                         drop = FALSE]

    group_bad <- if (nrow(groups) > 0) {
        groups$n_rows >= min_group_rows &
            (abs(groups$mean_z) > 0.35 |
                 (!is.na(groups$sd_z) & (groups$sd_z < 0.5 |
                                             groups$sd_z > 1.6)) |
                 groups$coverage_2se < 0.90)
    } else {
        logical()
    }

    failures <- bfpwr_sim_bind_failures(
        if (any(severe_bad, na.rm = TRUE)) {
            bfpwr_sim_failure(
                paste0(label, ":severe_prediction_interval"),
                sum(severe_bad, na.rm = TRUE),
                paste("rows must fall inside the severe binomial prediction interval at alpha",
                      severe_alpha)
            )
        }
    )

    list(summary = summary,
         group_summary = groups,
         outliers = outliers,
         diagnostic_rows = rows,
         failures = failures)
}

bfpwr_sim_validate_fixed_references <- function(fixture,
                                                designs,
                                                bf_priors,
                                                manifest_cases = NULL) {
    reference_rows <- bfpwr_sim_fixed_reference_table(
        fixture$tail_summary, designs, bf_priors,
        manifest_cases = manifest_cases)
    unsupported <- sum(!is.finite(reference_rows$reference_prob))
    diagnostics <- bfpwr_sim_mc_reference_diagnostics(
        reference_rows, label = "fixed_tail_package_reference")
    diagnostics$summary$unsupported_rows <- unsupported
    diagnostics$reference_rows <- reference_rows
    diagnostics$timings <- attr(reference_rows, "timings", exact = TRUE)
    if (is.null(diagnostics$timings)) {
        diagnostics$timings <- data.frame()
    }
    diagnostics
}

bfpwr_sim_select_fixed_reference_plot_rows <- function(tail_summary,
                                                       n_fractions = c(0.05,
                                                                       0.25,
                                                                       0.50,
                                                                       0.75,
                                                                       1.00)) {
    if (nrow(tail_summary) == 0) {
        return(tail_summary)
    }
    n_fractions <- sort(unique(as.numeric(n_fractions)))
    n_fractions <- n_fractions[is.finite(n_fractions)]
    if (length(n_fractions) == 0) {
        stop("n_fractions must contain at least one finite value")
    }
    n_fractions <- pmin(1, pmax(0, n_fractions))

    group_cols <- c("bf_prior_id", "threshold_id")
    missing_cols <- setdiff(c(group_cols, "n"), names(tail_summary))
    if (length(missing_cols) > 0) {
        stop("tail_summary is missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }

    keys <- unique(tail_summary[group_cols])
    selected <- vector("list", nrow(keys))
    for (i in seq_len(nrow(keys))) {
        keep <- tail_summary$bf_prior_id == keys$bf_prior_id[[i]] &
            tail_summary$threshold_id == keys$threshold_id[[i]]
        group <- tail_summary[keep, , drop = FALSE]
        group <- group[order(group$n), , drop = FALSE]
        n_values <- sort(unique(group$n))
        positions <- unique(pmax(
            1L,
            pmin(length(n_values),
                 round(1 + (length(n_values) - 1) * n_fractions))))
        selected_n <- n_values[positions]
        selected[[i]] <- group[group$n %in% selected_n, , drop = FALSE]
    }
    out <- do.call(rbind, selected)
    out <- out[order(out$family, out$bf_type, out$design_case_id,
                     out$bf_prior_id, out$threshold_id, out$n),
               , drop = FALSE]
    rownames(out) <- NULL
    out
}

bfpwr_sim_validate_sequential_references <- function(fixture,
                                                     designs,
                                                     bf_priors,
                                                     profile = c("curated",
                                                                 "none"),
                                                     strict = TRUE,
                                                     manifest_cases = NULL) {
    profile <- match.arg(profile)
    cases <- NULL
    if (!is.null(manifest_cases) && nrow(manifest_cases) > 0) {
        keep <- manifest_cases$validation_role %in%
            c("sequential_package_reference",
              "sequential_expected_sample_size_reference")
        cases <- manifest_cases[keep, , drop = FALSE]
    } else if (!is.null(manifest_cases)) {
        cases <- data.frame()
    }
    reference <- switch(
        fixture$spec$family,
        z = bfpwr_sim_z_sequential_reference_table(
            fixture = fixture,
            designs = designs,
            bf_priors = bf_priors,
            profile = profile,
            strict = strict,
            cases = cases),
        t = bfpwr_sim_t_sequential_reference_table(
            fixture = fixture,
            designs = designs,
            bf_priors = bf_priors,
            profile = profile,
            strict = strict,
            cases = cases),
        list(reference_rows = data.frame(),
             en_rows = data.frame(),
             timings = data.frame()))
    reference_rows <- reference$reference_rows
    if (nrow(reference_rows) == 0 &&
        !"reference_prob" %in% names(reference_rows)) {
        reference_rows <- data.frame(
            prob = numeric(),
            reference_prob = numeric(),
            nsim = integer(),
            n_event = integer())
    }
    unsupported <- if (nrow(reference_rows) == 0) {
        0L
    } else {
        sum(!is.finite(reference_rows$reference_prob))
    }
    diagnostics <- bfpwr_sim_mc_reference_diagnostics(
        reference_rows,
        label = "sequential_cumulative_package_reference",
        group_cols = c("family", "bf_type", "tail", "evidence_threshold",
                       "look_grid_name", "schedule_id"))
    diagnostics$summary$unsupported_rows <- unsupported

    en_rows <- reference$en_rows
    en_bad <- if (nrow(en_rows) == 0) {
        logical()
    } else {
        en_rows$abs_error > en_rows$tolerance
    }
    en_summary <- data.frame(
        label = "sequential_expected_sample_size_package_reference",
        rows_checked = nrow(en_rows),
        max_abs_error = if (nrow(en_rows) > 0) max(en_rows$abs_error) else NA_real_,
        max_abs_error_over_tolerance = if (nrow(en_rows) > 0) {
            max(en_rows$abs_error / en_rows$tolerance)
        } else NA_real_,
        failures = sum(en_bad),
        stringsAsFactors = FALSE)

    diagnostics$failures <- bfpwr_sim_bind_failures(
        diagnostics$failures,
        if (any(en_bad)) {
            bfpwr_sim_failure(
                "sequential_expected_sample_size_package_reference:tolerance",
                sum(en_bad),
                "expected sample size must agree with the package sequential reference within the Monte Carlo tolerance"
            )
        }
    )
    diagnostics$reference_rows <- reference_rows
    diagnostics$en_rows <- en_rows
    diagnostics$en_summary <- en_summary
    diagnostics$timings <- if (is.null(reference$timings)) {
        data.frame()
    } else reference$timings
    diagnostics
}

bfpwr_sim_validate_fixture_summary <- function(corpus_root,
                                               fixture_set_id,
                                               family = "z",
                                               mode = c("fixed", "sequential"),
                                               deterministic = TRUE,
                                               references = TRUE,
                                               case_set = "production",
                                               sequential_reference_profile = c("curated",
                                                                                "none"),
                                               manifest_cases = NULL,
                                               output_dir = NULL) {
    mode <- match.arg(mode)
    sequential_reference_profile <- match.arg(sequential_reference_profile)
    fixture <- bfpwr_sim_read_fixture_summary(
        corpus_root = corpus_root,
        fixture_set_id = fixture_set_id,
        family = family,
        mode = mode)
    manifest_failures <- bfpwr_sim_validate_fixture_manifest(
        fixture$fixture_dir,
        expected_files = bfpwr_sim_expected_fixture_files(mode))
    overview_failures <- bfpwr_sim_validate_fixture_overview(fixture, mode)
    deterministic_failures <- if (deterministic) {
        if (mode == "fixed") {
            bfpwr_sim_validate_fixed_fixture_deterministic(fixture)
        } else {
            bfpwr_sim_validate_sequential_fixture_deterministic(fixture)
        }
    } else {
        bfpwr_sim_bind_failures()
    }

    reference <- NULL
    if (references && mode == "fixed") {
        reference <- bfpwr_sim_validate_fixed_references(
            fixture,
            designs = bfpwr_sim_design_case_set(case_set),
            bf_priors = bfpwr_sim_bf_prior_case_set(case_set),
            manifest_cases = manifest_cases)
    } else if (references && mode == "sequential" &&
               family %in% c("z", "t") &&
               !identical(sequential_reference_profile, "none")) {
        reference <- bfpwr_sim_validate_sequential_references(
            fixture,
            designs = bfpwr_sim_design_case_set(case_set),
            bf_priors = bfpwr_sim_bf_prior_case_set(case_set),
            profile = sequential_reference_profile,
            strict = TRUE,
            manifest_cases = manifest_cases)
    }

    failures <- bfpwr_sim_bind_failures(
        manifest_failures,
        overview_failures,
        deterministic_failures,
        if (!is.null(reference)) reference$failures else bfpwr_sim_bind_failures()
    )
    summary <- data.frame(
        fixture_set_id = fixture_set_id,
        family = family,
        mode = mode,
        deterministic_failures = nrow(manifest_failures) +
            nrow(overview_failures) +
            nrow(deterministic_failures),
        reference_failures = if (is.null(reference)) 0L else
            nrow(reference$failures),
        reference_rows_checked = if (is.null(reference)) 0L else
            reference$summary$rows_checked,
        reference_rows_unsupported = if (is.null(reference)) NA_integer_ else
            reference$summary$unsupported_rows,
        passed = nrow(failures) == 0,
        stringsAsFactors = FALSE
    )

    result <- list(fixture = fixture,
                   summary = summary,
                   failures = failures,
                   reference = reference)
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        prefix <- file.path(output_dir, fixture_set_id)
        utils::write.csv(summary, paste0(prefix, "-summary.csv"),
                         row.names = FALSE)
        utils::write.csv(failures, paste0(prefix, "-failures.csv"),
                         row.names = FALSE)
        if (!is.null(reference)) {
            utils::write.csv(reference$summary,
                             paste0(prefix, "-reference-summary.csv"),
                             row.names = FALSE)
            utils::write.csv(reference$group_summary,
                             paste0(prefix, "-reference-groups.csv"),
                             row.names = FALSE)
            utils::write.csv(reference$outliers,
                             paste0(prefix, "-reference-outliers.csv"),
                             row.names = FALSE)
            optional <- c("-reference-diagnostic-rows.csv",
                          "-reference-en-summary.csv",
                          "-reference-en-rows.csv",
                          "-reference-timings.csv")
            for (suffix in optional) {
                path <- paste0(prefix, suffix)
                if (file.exists(path)) {
                    unlink(path)
                }
            }
            if (!is.null(reference$diagnostic_rows)) {
                utils::write.csv(reference$diagnostic_rows,
                                 paste0(prefix,
                                        "-reference-diagnostic-rows.csv"),
                                 row.names = FALSE)
            }
            if (!is.null(reference$en_summary)) {
                utils::write.csv(reference$en_summary,
                                 paste0(prefix, "-reference-en-summary.csv"),
                                 row.names = FALSE)
            }
            if (!is.null(reference$en_rows)) {
                utils::write.csv(reference$en_rows,
                                 paste0(prefix, "-reference-en-rows.csv"),
                                 row.names = FALSE)
            }
            if (!is.null(reference$timings)) {
                utils::write.csv(reference$timings,
                                 paste0(prefix, "-reference-timings.csv"),
                                 row.names = FALSE)
            }
        }
    }
    result
}
