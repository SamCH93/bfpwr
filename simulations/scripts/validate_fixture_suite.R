script_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    if (length(hit) > 0) {
        return(normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                             winslash = "/", mustWork = TRUE))
    }
    normalizePath("simulations/scripts/validate_fixture_suite.R",
                  winslash = "/", mustWork = FALSE)
}

source(file.path(dirname(script_path()), "common.R"))

validate_numeric_counts <- function(df, failures, prefix) {
    if ("prob" %in% names(df)) {
        failures <- add_failure(
            failures, "error", paste0(prefix, ":probability_bounds"),
            count_bad_prob(df$prob), "probabilities must be finite and in [0, 1]"
        )
    }
    if (all(c("pH1", "pH0", "pInc") %in% names(df))) {
        failures <- add_failure(
            failures, "error", paste0(prefix, ":decision_sum"),
            count_not_close(df$pH1 + df$pH0 + df$pInc, 1),
            "H1, H0, and inconclusive probabilities must sum to one"
        )
    }
    if (all(c("cum_pH1", "cum_pH0", "cum_pInc") %in% names(df))) {
        failures <- add_failure(
            failures, "error", paste0(prefix, ":cumulative_sum"),
            count_not_close(df$cum_pH1 + df$cum_pH0 + df$cum_pInc, 1),
            "cumulative H1, H0, and continuation probabilities must sum to one"
        )
    }
    if (all(c("nH1", "nH0", "nInc", "nsim") %in% names(df))) {
        failures <- add_failure(
            failures, "error", paste0(prefix, ":decision_counts"),
            count_not_close(df$nH1 + df$nH0 + df$nInc, df$nsim),
            "decision counts must sum to nsim"
        )
    }
    if (all(c("n_nan_log_bf01", "n_na_log_bf01") %in% names(df))) {
        failures <- add_failure(
            failures, "error", paste0(prefix, ":no_nan_log_bf01"),
            safe_sum(df$n_nan_log_bf01 > 0),
            "NaN log_bf01 values are not allowed in finalized fixtures"
        )
        failures <- add_failure(
            failures, "error", paste0(prefix, ":no_na_log_bf01"),
            safe_sum(df$n_na_log_bf01 > 0),
            "NA log_bf01 values are not allowed in finalized fixtures"
        )
    }
    failures
}

validate_search_summary <- function(df, failures, prefix) {
    if (!is.data.frame(df) || nrow(df) == 0) {
        return(failures)
    }
    if (all(c("achieved", "n_found", "prob_found", "target_prob") %in%
            names(df))) {
        reached <- isTRUE(df$achieved) | (!is.na(df$achieved) & df$achieved)
        failures <- add_failure(
            failures, "error", paste0(prefix, ":achieved_has_n"),
            safe_sum(reached & !is.finite(df$n_found)),
            "achieved search rows must have a finite sample size"
        )
        failures <- add_failure(
            failures, "error", paste0(prefix, ":achieved_reaches_target"),
            safe_sum(reached & is.finite(df$prob_found) &
                         df$prob_found + 1e-12 < df$target_prob),
            "achieved search rows must meet the target probability"
        )
    }
    failures
}

reference_contract <- function(fixture_set_id, role, validation_dir) {
    ref_files <- list.files(
        validation_dir,
        pattern = paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1",
                                  fixture_set_id),
                         ".*reference.*summary\\.csv$"),
        full.names = TRUE
    )

    rows_checked <- 0L
    unsupported_rows <- NA_integer_
    failures <- 0L
    diagnostics <- 0L

    for (path in ref_files) {
        ref <- read_csv_if_exists(path)
        if (is.null(ref) || nrow(ref) == 0) {
            next
        }
        if ("rows_checked" %in% names(ref)) {
            rows_checked <- rows_checked + safe_sum(ref$rows_checked)
        }
        if ("unsupported_rows" %in% names(ref)) {
            unsupported_rows <- if (is.na(unsupported_rows)) 0L else unsupported_rows
            unsupported_rows <- unsupported_rows + safe_sum(ref$unsupported_rows)
        }

        severe <- if ("severe_interval_failures" %in% names(ref)) {
            safe_sum(ref$severe_interval_failures)
        } else 0L
        explicit <- if ("failures" %in% names(ref)) safe_sum(ref$failures) else 0L
        total <- severe + explicit
        if (identical(role, "diagnostic_approximation")) {
            diagnostics <- diagnostics + total
        } else {
            failures <- failures + total
        }
    }

    if (identical(role, "integrity_only") && length(ref_files) == 0) {
        unsupported_rows <- NA_integer_
    } else if (is.na(unsupported_rows)) {
        unsupported_rows <- 0L
    }

    list(
        rows_checked = rows_checked,
        unsupported_rows = unsupported_rows,
        reference_failures = failures,
        diagnostic_issues = diagnostics
    )
}

validate_fixture <- function(fixture_dir, output_dir) {
    overview <- read_overview(fixture_dir)
    fixture_set_id <- overview$fixture_set_id[[1]]
    family <- overview$family[[1]]
    mode <- overview$mode[[1]]
    role <- fixture_role(fixture_set_id, family, mode)
    failures <- empty_failures()

    required <- c("spec.rds", "manifest.csv", "overview.csv")
    if (identical(mode, "fixed")) {
        required <- c(required, "tail-summary.rds", "decision-summary.rds",
                      "search-summary.rds")
    } else if (identical(mode, "sequential")) {
        required <- c(required, "cumulative-summary.rds", "final-summary.rds",
                      "search-summary.rds")
    } else {
        failures <- add_failure(failures, "error", "fixture:unknown_mode", 1,
                                "fixture mode must be fixed or sequential")
    }

    missing <- required[!file.exists(file.path(fixture_dir, required))]
    failures <- add_failure(
        failures, "error", "fixture:missing_required_file", length(missing),
        if (length(missing)) paste(missing, collapse = ", ") else ""
    )

    if (identical(mode, "fixed") && length(missing) == 0) {
        tail_summary <- readRDS(file.path(fixture_dir, "tail-summary.rds"))
        decision_summary <- readRDS(file.path(fixture_dir, "decision-summary.rds"))
        search_summary <- readRDS(file.path(fixture_dir, "search-summary.rds"))
        failures <- validate_numeric_counts(tail_summary, failures, "tail_summary")
        failures <- validate_numeric_counts(decision_summary, failures,
                                            "decision_summary")
        failures <- validate_search_summary(search_summary, failures,
                                            "search_summary")
    }

    if (identical(mode, "sequential") && length(missing) == 0) {
        cumulative_summary <- readRDS(file.path(fixture_dir,
                                                "cumulative-summary.rds"))
        final_summary <- readRDS(file.path(fixture_dir, "final-summary.rds"))
        search_summary <- readRDS(file.path(fixture_dir, "search-summary.rds"))
        failures <- validate_numeric_counts(cumulative_summary, failures,
                                            "cumulative_summary")
        failures <- validate_numeric_counts(final_summary, failures,
                                            "final_summary")
        failures <- validate_search_summary(search_summary, failures,
                                            "search_summary")
    }

    ref <- reference_contract(fixture_set_id, role, output_dir)
    if (identical(role, "package_reference") && ref$rows_checked == 0) {
        failures <- add_failure(
            failures, "error", "reference:missing_reference_summary", 1,
            "package-reference fixtures must have materialized reference summaries"
        )
    }
    deterministic_failures <- safe_sum(failures$severity == "error")
    passed <- deterministic_failures == 0 && ref$reference_failures == 0

    summary <- data.frame(
        fixture_set_id = fixture_set_id,
        family = family,
        mode = mode,
        validation_role = role,
        deterministic_failures = deterministic_failures,
        reference_failures = ref$reference_failures,
        diagnostic_issues = ref$diagnostic_issues,
        reference_rows_checked = ref$rows_checked,
        reference_rows_unsupported = ref$unsupported_rows,
        passed = passed,
        stringsAsFactors = FALSE
    )

    if (ref$diagnostic_issues > 0) {
        failures <- add_failure(
            failures, "diagnostic",
            "reference:known_approximation_gap",
            ref$diagnostic_issues,
            paste0(
                "known fixed-t exact-t simulation gap; ptbf01() verifies the ",
                "documented normal-effect approximation, so this is not a ",
                "package-reference failure"
            )
        )
    }

    write_csv(summary, file.path(output_dir, paste0(fixture_set_id,
                                                    "-summary.csv")))
    write_csv(failures, file.path(output_dir, paste0(fixture_set_id,
                                                     "-failures.csv")))
    list(summary = summary, failures = failures)
}

main <- function() {
    args <- parse_args()
    corpus_root <- normalizePath(arg_value(args, "corpus-root",
                                           "simulations/corpus/v1"),
                                 winslash = "/", mustWork = TRUE)
    output_dir <- normalizePath(arg_value(args, "output-dir",
                                          file.path(corpus_root,
                                                    "fixture-validation")),
                                winslash = "/", mustWork = FALSE)
    ensure_dir(output_dir)
    fixture_set <- arg_value(args, "fixture-set", NULL)
    write_suite <- !arg_flag(args, "no-suite")

    dirs <- find_fixture_dirs(corpus_root, fixture_set)
    results <- lapply(dirs, validate_fixture, output_dir = output_dir)
    summaries <- do.call(rbind, lapply(results, `[[`, "summary"))
    failures <- do.call(rbind, Map(function(result) {
        f <- result$failures
        if (nrow(f) == 0) {
            return(data.frame())
        }
        cbind(fixture_set_id = result$summary$fixture_set_id[[1]], f,
              stringsAsFactors = FALSE)
    }, results))
    if (is.null(failures) || nrow(failures) == 0) {
        failures <- data.frame(
            fixture_set_id = character(),
            severity = character(),
            rule = character(),
            n_failed = integer(),
            details = character(),
            stringsAsFactors = FALSE
        )
    }

    if (write_suite) {
        write_csv(summaries, file.path(output_dir, "fixture-suite-summary.csv"))
        write_csv(failures, file.path(output_dir, "fixture-suite-failures.csv"))
    }

    n_failed <- sum(!summaries$passed)
    cat("Checked ", nrow(summaries), " simulation output sets: ",
        sum(summaries$passed), " passed, ", n_failed, " failed.\n",
        sep = "")
    if (sum(summaries$reference_failures, na.rm = TRUE) > 0) {
        cat("Package-reference failures: ",
            sum(summaries$reference_failures, na.rm = TRUE), "\n", sep = "")
    }
    if (sum(summaries$diagnostic_issues, na.rm = TRUE) > 0) {
        cat("Documented approximation notes: ",
            sum(summaries$diagnostic_issues, na.rm = TRUE), "\n", sep = "")
    }
    cat("Detailed output: ", output_dir, "\n", sep = "")
    if (any(!summaries$passed)) {
        stop("fixture validation failed for: ",
             paste(summaries$fixture_set_id[!summaries$passed], collapse = ", "),
             call. = FALSE)
    }
    invisible(summaries)
}

main()
