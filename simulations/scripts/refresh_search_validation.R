script_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    if (length(hit) > 0) {
        return(normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                             winslash = "/", mustWork = TRUE))
    }
    normalizePath("simulations/scripts/refresh_search_validation.R",
                  winslash = "/", mustWork = FALSE)
}

source(file.path(dirname(script_path()), "common.R"))

summarize_search <- function(comparison) {
    split_key <- interaction(comparison$family, comparison$mode,
                             comparison$bf_type,
                             comparison$search_profile_id,
                             comparison$target_source,
                             drop = TRUE, sep = "\r")
    rows <- lapply(split(comparison, split_key), function(x) {
        reference_z <- x$package_reference_prob_z
        evaluated <- !is.na(x$package_reached)
        reached <- x$package_reached %in% TRUE
        certified <- reached & x$package_first_crossing_certified %in% TRUE
        package_errors <- evaluated & !is.na(x$package_error) &
            nzchar(x$package_error)
        comparable <- reached
        reference_comparable <- !is.na(reference_z)
        n_comparable <- reached & is.finite(x$package_n) &
            is.finite(x$simulation_n)
        n_different <- n_comparable & x$package_n != x$simulation_n
        data.frame(
            family = x$family[[1]],
            mode = x$mode[[1]],
            bf_type = x$bf_type[[1]],
            search_profile_id = x$search_profile_id[[1]],
            target_source = x$target_source[[1]],
            rows = nrow(x),
            simulation_achieved = safe_sum(x$sim_achieved),
            evaluated_rows = safe_sum(evaluated),
            package_reached = safe_sum(x$package_reached %in% TRUE),
            package_certified = safe_sum(certified),
            package_errors = safe_sum(package_errors),
            package_not_reached = safe_sum(evaluated & !reached),
            package_not_certified = safe_sum(reached & !certified),
            n_different = safe_sum(n_different),
            comparable = safe_sum(comparable),
            n_comparable = safe_sum(n_comparable),
            reference_comparable = safe_sum(reference_comparable),
            median_n_ratio = if (any(comparable, na.rm = TRUE)) {
                stats::median(x$n_ratio[comparable], na.rm = TRUE)
            } else NA_real_,
            p90_abs_fractional_error = if (any(comparable, na.rm = TRUE)) {
                stats::quantile(abs(x$fractional_error[comparable]),
                                probs = 0.9, na.rm = TRUE, names = FALSE)
            } else NA_real_,
            max_abs_fractional_error = if (any(comparable, na.rm = TRUE)) {
                max(abs(x$fractional_error[comparable]), na.rm = TRUE)
            } else NA_real_,
            median_abs_reference_z = if (any(reference_comparable)) {
                stats::median(abs(reference_z[reference_comparable]),
                              na.rm = TRUE)
            } else NA_real_,
            p90_abs_reference_z = if (any(reference_comparable)) {
                stats::quantile(abs(reference_z[reference_comparable]),
                                probs = 0.9, na.rm = TRUE, names = FALSE)
            } else NA_real_,
            max_abs_reference_z = if (any(reference_comparable)) {
                max(abs(reference_z[reference_comparable]), na.rm = TRUE)
            } else NA_real_,
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out[order(out$family, out$mode, out$bf_type, out$search_profile_id), ]
}

classify_search_rows <- function(summary) {
    role <- rep("package_smoke", nrow(summary))
    role[summary$target_source == "simulation_anchor" &
             summary$family == "t" &
             summary$reference_comparable == 0] <- "simulation_only"
    role
}

manifest_rows <- function(comparison) {
    if (!"manifest_role" %in% names(comparison)) {
        return(comparison[FALSE, , drop = FALSE])
    }
    role <- comparison$manifest_role
    comparison[!is.na(role) & nzchar(role), , drop = FALSE]
}

summarize_manifest_coverage <- function(comparison) {
    manifest <- manifest_rows(comparison)
    if (nrow(manifest) == 0) {
        return(data.frame())
    }
    finite_min <- function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0) NA_real_ else min(v)
    }
    finite_max <- function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0) NA_real_ else max(v)
    }
    split_key <- interaction(manifest$expected_package_function,
                             manifest$manifest_role,
                             manifest$simulation_status,
                             drop = TRUE, sep = "\r")
    rows <- lapply(split(manifest, split_key), function(x) {
        reached <- x$package_reached %in% TRUE
        certified <- x$package_first_crossing_certified %in% TRUE
        evaluated <- !is.na(x$package_reached)
        data.frame(
            package_function = x$expected_package_function[[1]],
            manifest_role = x$manifest_role[[1]],
            simulation_status = x$simulation_status[[1]],
            rows = nrow(x),
            simulation_achieved = safe_sum(x$sim_achieved),
            package_reached = safe_sum(reached),
            package_certified = safe_sum(reached & certified),
            package_errors = safe_sum(evaluated &
                                          !is.na(x$package_error) &
                                          nzchar(x$package_error)),
            bf_types = paste(sort(unique(x$bf_type)), collapse = ", "),
            evidence = paste(sort(unique(x$evidence)), collapse = ", "),
            thresholds = paste(sort(unique(x$evidence_threshold)),
                               collapse = ", "),
            schedules = length(unique(x$schedule_id)),
            design_cases = length(unique(x$design_case_id)),
            bf_priors = length(unique(x$bf_prior_id)),
            target_probs = length(unique(x$target_prob)),
            target_prob_values = paste(sort(unique(x$target_prob)),
                                       collapse = ", "),
            look_counts = paste(sort(unique(x$look_count)),
                                collapse = ", "),
            min_n = finite_min(x$simulation_n),
            max_n = finite_max(x$simulation_n),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out[order(out$package_function, out$manifest_role,
              out$simulation_status), , drop = FALSE]
}

validate_search <- function(comparison, summary, metadata, scope) {
    status <- empty_failures()

    package_rows <- !is.na(comparison$package_reached)
    status <- add_failure(
        status, "error", "package_rows:errors",
        safe_sum(package_rows & nzchar(comparison$package_error)),
        "materialized package-search rows must not contain evaluator errors"
    )
    crossing_rows <- comparison$target_source != "simulation_anchor"
    status <- add_failure(
        status, "error", "package_rows:first_crossing_certified",
        safe_sum(package_rows &
                     crossing_rows &
                     comparison$package_reached %in% TRUE &
                     comparison$package_first_crossing_certified %in% FALSE),
        "reached package-search rows must certify the first crossing"
    )
    status <- add_failure(
        status, "error", "package_rows:target_reached",
        safe_sum(package_rows &
                     comparison$package_reached %in% TRUE &
                     is.finite(comparison$package_actual_prob) &
                     comparison$package_actual_prob + 1e-10 <
                         comparison$target_prob),
        "reached package-search rows must meet the requested target"
    )

    reference_rows <- !is.na(comparison$package_reference_prob_z)
    status <- add_failure(
        status, "error", "reference_rows:z_abs_gt_4",
        safe_sum(reference_rows &
                     abs(comparison$package_reference_prob_z) > 4),
        "package probability at exact anchor schedules should stay within 4 MCSE"
    )

    manifest <- manifest_rows(comparison)
    if (nrow(manifest) > 0) {
        expected_functions <- sort(unique(manifest$expected_package_function))
        for (fn in expected_functions) {
            rows <- manifest[manifest$expected_package_function == fn,
                             , drop = FALSE]
            required <- rows$manifest_role == "required_n_search"
            status <- add_failure(
                status, "error",
                paste0("manifest:", fn, ":required_not_evaluated"),
                safe_sum(required & is.na(rows$package_reached)),
                paste0(fn, " required sequential n-finding manifest rows ",
                       "must be evaluated by the package search")
            )
            status <- add_failure(
                status, "error",
                paste0("manifest:", fn, ":required_not_reached"),
                safe_sum(required &
                             !(rows$package_reached %in% TRUE)),
                paste0(fn, " required sequential n-finding cases must ",
                       "be reached by the package search")
            )
            status <- add_failure(
                status, "error",
                paste0("manifest:", fn, ":required_not_certified"),
                safe_sum(required &
                             rows$package_reached %in% TRUE &
                             !(rows$package_first_crossing_certified %in%
                                   TRUE)),
                paste0(fn, " required sequential n-finding cases must ",
                       "certify the first crossing")
            )
        }

        diagnostics <- manifest$manifest_role == "boundary_diagnostic"
        status <- add_failure(
            status, "warning", "manifest:diagnostics_reached",
            safe_sum(diagnostics & manifest$package_reached %in% TRUE),
            paste0("boundary diagnostic rows are non-achieved simulation ",
                   "cases; a package crossing here should be reviewed")
        )
        status <- add_failure(
            status, "warning", "manifest:internal_best_nonachieved",
            safe_sum(diagnostics &
                         manifest$simulation_status ==
                             "not_achieved_internal_best"),
            paste0("some non-achieved simulation diagnostics did not end ",
                   "at the search boundary and should be interpreted case ",
                   "by case")
        )
    }

    if (identical(scope, "smoke")) {
        no_package_profiles <- summary$comparable == 0 &
            summary$target_source != "simulation_anchor"
        status <- add_failure(
            status, "warning", "scope:profiles_without_package_rows",
            safe_sum(no_package_profiles),
            paste0(
                "search validation is a stratified smoke suite; not every ",
                "simulation crossing has a package-search row"
            )
        )
    }

    roles <- classify_search_rows(summary)
    summary$validation_role <- roles
    list(status = status, summary = summary)
}

main <- function() {
    args <- parse_args()
    corpus_root <- normalizePath(arg_value(args, "corpus-root",
                                           "simulations/corpus/v1"),
                                 winslash = "/", mustWork = TRUE)
    validation_dir <- normalizePath(arg_value(args, "validation-dir",
                                              file.path(corpus_root,
                                                        "search-validation")),
                                    winslash = "/", mustWork = TRUE)
    scope <- match.arg(arg_value(args, "scope", "smoke"),
                       choices = c("smoke", "full"))
    comparison_path <- file.path(validation_dir,
                                 "search-validation-comparison.rds")
    comparison_bundle <- read_rds_required(comparison_path)
    comparison <- comparison_bundle$comparison
    metadata <- comparison_bundle$metadata
    summary <- summarize_search(comparison)
    manifest_coverage <- summarize_manifest_coverage(comparison)
    checked <- validate_search(comparison, summary, metadata, scope)
    summary <- checked$summary
    status <- checked$status

    overview <- data.frame(
        generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        validation_scope = scope,
        simulation_rows = nrow(comparison_bundle$simulation),
        package_rows = nrow(comparison_bundle$package),
        comparison_rows = nrow(comparison),
        comparable_rows = safe_sum(comparison$package_reached %in% TRUE),
        fixed_rows = safe_sum(comparison$mode == "fixed"),
        sequential_rows = safe_sum(comparison$mode == "sequential"),
        sequential_manifest_rows = nrow(manifest_rows(comparison)),
        package_recomputed_rows = if (!is.null(metadata$package_recomputed_rows)) {
            metadata$package_recomputed_rows
        } else NA_integer_,
        package_recompute_elapsed_seconds = if (!is.null(
            metadata$package_recompute_elapsed_seconds)) {
            metadata$package_recompute_elapsed_seconds
        } else NA_real_,
        search_diagnostics_evaluated = if (!is.null(
            metadata$sequential_search_diagnostics_evaluated)) {
            metadata$sequential_search_diagnostics_evaluated
        } else NA,
        package_rows_policy = if (!is.null(metadata$package_rows)) {
            metadata$package_rows
        } else NA_character_,
        stringsAsFactors = FALSE
    )

    write_csv(overview, file.path(validation_dir,
                                  "search-validation-overview.csv"))
    write_csv(summary, file.path(validation_dir,
                                 "search-validation-summary.csv"))
    write_csv(manifest_coverage,
              file.path(validation_dir,
                        "sequential-search-validation-coverage.csv"))
    write_csv(status, file.path(validation_dir,
                                "search-validation-status.csv"))

    cat("Compared ", nrow(comparison), " search rows; ",
        safe_sum(comparison$package_reached %in% TRUE),
        " package searches reached a target.\n", sep = "")
    cat("Package search rows recomputed: ",
        overview$package_recomputed_rows[[1]], "\n", sep = "")
    if (is.finite(overview$package_recompute_elapsed_seconds[[1]])) {
        cat("Package search recomputation time: ",
            round(overview$package_recompute_elapsed_seconds[[1]], 1),
            " seconds.\n", sep = "")
    }
    if (nrow(manifest_coverage) > 0) {
        primary <- manifest_coverage[
            manifest_coverage$manifest_role == "required_n_search",
            , drop = FALSE]
        cat("Primary sequential search checks: ",
            sum(primary$package_certified, na.rm = TRUE), " of ",
            sum(primary$rows, na.rm = TRUE),
            " package searches certified.\n", sep = "")
    }
    active_status <- status[status$n_failed > 0, , drop = FALSE]
    if (nrow(active_status) > 0) {
        cat("Sample-size search notes: ",
            paste(paste(active_status$n_failed, active_status$severity,
                        "rows"),
                  collapse = "; "),
            ".\n", sep = "")
    }
    cat("Detailed output: ", validation_dir, "\n", sep = "")
    if (any(status$severity == "error")) {
        stop("search validation failed; see search-validation-status.csv",
             call. = FALSE)
    }
    invisible(summary)
}

main()
