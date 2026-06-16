bfpwr_sim_search_manifest_file <- function(repo_root = ".") {
    file.path(repo_root, "simulations", "registry", "manifests",
              "sequential-search-validation-cases.csv")
}

bfpwr_sim_bind_rows_fill <- function(...) {
    rows <- list(...)
    if (length(rows) == 1L && is.list(rows[[1]]) &&
        !is.data.frame(rows[[1]])) {
        rows <- rows[[1]]
    }
    rows <- rows[vapply(rows, is.data.frame, logical(1))]
    rows <- rows[vapply(rows, nrow, integer(1)) > 0]
    if (length(rows) == 0) {
        return(data.frame())
    }
    cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
    rows <- lapply(rows, function(x) {
        missing <- setdiff(cols, names(x))
        for (col in missing) {
            x[[col]] <- NA
        }
        x[cols]
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

bfpwr_sim_sequential_package_function <- function(family) {
    out <- ifelse(family == "t", "ntbf01seq",
                  ifelse(family == "z", "nbf01seq", NA_character_))
    if (any(is.na(out))) {
        stop("unsupported sequential search validation family: ",
             paste(unique(family[is.na(out)]), collapse = ", "),
             call. = FALSE)
    }
    out
}

bfpwr_sim_relative_corpus_path <- function(path, corpus_root) {
    normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    normalized_root <- normalizePath(corpus_root, winslash = "/",
                                     mustWork = FALSE)
    prefix <- paste0(normalized_root, "/")
    ifelse(startsWith(normalized_path, prefix),
           substring(normalized_path, nchar(prefix) + 1L),
           normalized_path)
}

bfpwr_sim_search_sequential_summary_fast <- function(cumulative_summary,
                                                     search_targets) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        return(bfpwr_sim_search_sequential_summary(
            cumulative_summary = cumulative_summary,
            search_targets = search_targets))
    }

    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "schedule_id",
                    "schedule_family_id", "start_n", "increment", "n_looks",
                    "max_n", "threshold_pair_id", "evidence_threshold",
                    "k1", "k0")
    target_cols <- c("search_target_id", "evidence", "target_prob")
    dt <- data.table::as.data.table(cumulative_summary)
    targets <- data.table::as.data.table(
        search_targets[, c("evidence_threshold", target_cols),
                       drop = FALSE])
    data.table::setorderv(dt, c(group_cols, "look"))
    joined <- dt[targets, on = "evidence_threshold", allow.cartesian = TRUE]
    joined[, prob := ifelse(evidence == "H1", cum_pH1, cum_pH0)]
    joined[, mcse := ifelse(evidence == "H1", mcse_cum_pH1, mcse_cum_pH0)]

    by_cols <- c(group_cols, target_cols)
    out <- joined[, {
        hit <- which(prob >= target_prob[[1]])
        achieved <- length(hit) > 0
        first <- if (achieved) hit[[1]] else NA_integer_
        previous <- if (achieved && first > 1L) first - 1L else NA_integer_
        finite_prob <- is.finite(prob)
        best <- if (any(finite_prob)) which.max(prob) else NA_integer_
        list(
            criterion = "first_look_with_cumulative_prob_ge_target",
            achieved = achieved,
            look_found = if (achieved) look[[first]] else NA_integer_,
            n_found = if (achieved) n[[first]] else NA_integer_,
            prob_found = if (achieved) prob[[first]] else NA_real_,
            mcse_found = if (achieved) mcse[[first]] else NA_real_,
            look_previous = if (is.na(previous)) {
                NA_integer_
            } else {
                look[[previous]]
            },
            n_previous = if (is.na(previous)) NA_integer_ else n[[previous]],
            prob_previous = if (is.na(previous)) {
                NA_real_
            } else {
                prob[[previous]]
            },
            best_look = if (is.na(best)) NA_integer_ else look[[best]],
            best_n = if (is.na(best)) NA_integer_ else n[[best]],
            best_prob = if (is.na(best)) NA_real_ else prob[[best]],
            max_n_searched = max(n, na.rm = TRUE)
        )
    }, by = by_cols]
    as.data.frame(out)
}

bfpwr_sim_read_sequential_search_summaries <- function(
        corpus_root,
        families = c("z", "t"),
        search_targets = NULL) {
    if (is.null(search_targets)) {
        search_targets <- bfpwr_sim_fixture_search_targets_core()
    }
    rows <- list()
    idx <- 0L
    for (family in families) {
        root <- file.path(corpus_root, "fixtures", family, "sequential")
        if (!dir.exists(root)) {
            next
        }
        files <- list.files(root, pattern = "^cumulative-summary\\.rds$",
                            recursive = TRUE, full.names = TRUE)
        for (file in files) {
            cumulative <- readRDS(file)
            if (!is.data.frame(cumulative) || nrow(cumulative) == 0) {
                next
            }
            x <- bfpwr_sim_search_sequential_summary_fast(
                cumulative_summary = cumulative,
                search_targets = search_targets)
            if (!is.data.frame(x) || nrow(x) == 0) next
            idx <- idx + 1L
            ## Sequential fixture specs keep package-search targets disabled.
            ## The validation manifest derives first crossings from the
            ## cumulative summaries so package n-search and simulation rows
            ## use the same current criterion.
            x$source_cumulative_summary <-
                bfpwr_sim_relative_corpus_path(file, corpus_root)
            rows[[idx]] <- x
        }
    }
    out <- bfpwr_sim_bind_rows_fill(rows)
    if (nrow(out) == 0) {
        return(out)
    }
    out$package_function <- bfpwr_sim_sequential_package_function(out$family)
    out$mode <- "sequential"
    out
}

bfpwr_sim_search_validation_key <- function(x) {
    key_cols <- c("package_function", "fixture_set_id", "bf_prior_id",
                  "design_case_id", "schedule_id", "evidence_threshold",
                  "evidence", "target_prob", "start_n", "increment",
                  "n_looks", "max_n")
    paste("seq-search-v1", do.call(paste, c(x[key_cols], sep = "|")),
          sep = "|")
}

bfpwr_sim_stratified_select <- function(x,
                                        n,
                                        strata_cols,
                                        order_cols) {
    if (nrow(x) <= n) {
        return(x)
    }
    keep_order_cols <- order_cols[order_cols %in% names(x)]
    if (length(keep_order_cols) > 0) {
        ord <- do.call(order, c(x[keep_order_cols], list(na.last = TRUE)))
        x <- x[ord, , drop = FALSE]
    }
    strata_cols <- strata_cols[strata_cols %in% names(x)]
    x$.selection_stratum <- if (length(strata_cols) == 0) {
        "all"
    } else {
        bfpwr_sim_row_key(x, strata_cols)
    }
    groups <- split(seq_len(nrow(x)), x$.selection_stratum)
    group_names <- sort(names(groups))
    if (length(group_names) > n) {
        idx <- floor(seq(1, length(group_names) + 1, length.out = n + 1))
        idx <- unique(pmin(idx[seq_len(n)], length(group_names)))
        if (length(idx) < n) {
            idx <- c(idx, setdiff(seq_along(group_names), idx))
            idx <- idx[seq_len(n)]
        }
        group_names <- group_names[idx]
    }
    selected <- integer()
    cursor <- setNames(rep(1L, length(groups)), names(groups))
    while (length(selected) < n) {
        progressed <- FALSE
        for (name in group_names) {
            group <- groups[[name]]
            j <- cursor[[name]]
            if (j <= length(group)) {
                selected <- c(selected, group[[j]])
                cursor[[name]] <- j + 1L
                progressed <- TRUE
                if (length(selected) == n) {
                    break
                }
            }
        }
        if (!progressed) {
            break
        }
    }
    out <- x[selected, , drop = FALSE]
    out$.selection_stratum <- NULL
    rownames(out) <- NULL
    out
}

bfpwr_sim_add_search_selection_fields <- function(x,
                                                  manifest_role) {
    achieved <- x$achieved %in% TRUE
    boundary_limited <- !achieved & is.finite(x$best_n) &
        is.finite(x$max_n_searched) & x$best_n >= x$max_n_searched
    simulation_status <- rep("achieved", nrow(x))
    simulation_status[!achieved & boundary_limited] <-
        "not_achieved_at_search_bound"
    simulation_status[!achieved & !boundary_limited] <-
        "not_achieved_internal_best"

    x$manifest_role <- manifest_role
    x$simulation_status <- simulation_status
    x$expected_package_function <- x$package_function
    x$target_source <- "sequential_schedule_crossing"
    x$search_profile_id <- paste0("seq-", x$schedule_family_id)
    x$search_profile_label <- paste0("look every ", x$increment,
                                     " from N = ", x$start_n)
    x$simulation_schedule <- "incremental_schedule_family"
    x$look_count <- x$n_looks
    x$grid_step <- x$increment
    x$reference_n <- NA_integer_
    x$reference_schedule <- x$schedule_family_id
    x$reference_previous_n <- NA_integer_
    x$reference_previous_prob <- NA_real_
    x$reference_growth <- NA
    x$schedule_start_n <- x$start_n
    x$schedule_increment <- x$increment
    x$sim_achieved <- achieved
    x$simulation_n <- x$n_found
    x$simulation_prob <- x$prob_found
    x$simulation_mcse <- x$mcse_found
    x$simulation_previous_n <- x$n_previous
    x$simulation_previous_prob <- x$prob_previous
    x$simulation_best_n <- x$best_n
    x$simulation_best_prob <- x$best_prob
    x$simulation_max_n_searched <- x$max_n_searched
    x$package_nrange_lower <- x$start_n
    x$package_nrange_upper <- x$max_n_searched
    x$selection_stratum <- paste(x$package_function, x$bf_type, x$evidence,
                                 paste0("bf", x$evidence_threshold),
                                 paste0("p", x$target_prob),
                                 x$schedule_family_id,
                                 paste0("looks", x$n_looks),
                                 sep = ":")
    x$package_key <- bfpwr_sim_search_validation_key(x)
    x
}

bfpwr_sim_build_sequential_search_manifest <- function(
        corpus_root,
        families = c("z", "t"),
        achieved_per_function = 50L,
        diagnostic_per_function = 10L,
        min_achieved_per_function = 50L,
        required_mcse_margin = 2) {
    all_rows <- bfpwr_sim_read_sequential_search_summaries(
        corpus_root = corpus_root, families = families)
    if (nrow(all_rows) == 0) {
        stop("no sequential search summaries found under ", corpus_root,
             call. = FALSE)
    }
    required_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                       "design_case_id", "look_grid_name", "schedule_id",
                       "schedule_family_id", "start_n", "increment",
                       "n_looks", "max_n", "threshold_pair_id",
                       "evidence_threshold", "k1", "k0", "evidence",
                       "target_prob", "achieved", "n_found", "prob_found",
                       "mcse_found", "n_previous", "prob_previous",
                       "best_n", "best_prob", "max_n_searched",
                       "package_function")
    missing <- setdiff(required_cols, names(all_rows))
    if (length(missing) > 0) {
        stop("sequential search summaries are missing columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }

    strata_cols <- c("bf_type", "evidence", "evidence_threshold",
                     "target_prob", "schedule_family_id", "n_looks",
                     "look_grid_name")
    order_cols <- c("bf_type", "evidence", "evidence_threshold",
                    "target_prob", "schedule_family_id", "n_looks",
                    "design_case_id", "bf_prior_id")
    selected <- list()
    idx <- 0L
    functions <- sort(unique(all_rows$package_function))
    for (fn in functions) {
        rows <- all_rows[all_rows$package_function == fn, , drop = FALSE]
        achieved <- rows[rows$achieved %in% TRUE &
                             is.finite(rows$n_found) &
                             is.finite(rows$prob_found) &
                             is.finite(rows$mcse_found) &
                             rows$prob_found - required_mcse_margin *
                                 rows$mcse_found >= rows$target_prob,
                         , drop = FALSE]
        if (nrow(achieved) < min_achieved_per_function) {
            stop("only ", nrow(achieved), " robust achieved sequential search rows ",
                 "available for ", fn, "; need at least ",
                 min_achieved_per_function, call. = FALSE)
        }
        idx <- idx + 1L
        selected[[idx]] <- bfpwr_sim_add_search_selection_fields(
            bfpwr_sim_stratified_select(
                achieved,
                n = as.integer(achieved_per_function),
                strata_cols = strata_cols,
                order_cols = order_cols),
            manifest_role = "required_n_search")

        nonachieved <- rows[!(rows$achieved %in% TRUE), , drop = FALSE]
        if (diagnostic_per_function > 0 && nrow(nonachieved) > 0) {
            idx <- idx + 1L
            selected[[idx]] <- bfpwr_sim_add_search_selection_fields(
                bfpwr_sim_stratified_select(
                    nonachieved,
                    n = as.integer(diagnostic_per_function),
                    strata_cols = c(strata_cols, "best_n"),
                    order_cols = order_cols),
                manifest_role = "boundary_diagnostic")
        }
    }
    manifest <- bfpwr_sim_bind_rows_fill(selected)
    manifest$search_case_id <- sprintf("seq-search-%03d", seq_len(nrow(manifest)))
    manifest$manifest_version <- "sequential-search-v1"
    role_order <- match(manifest$manifest_role,
                        c("required_n_search", "boundary_diagnostic"))
    manifest <- manifest[order(manifest$package_function,
                               role_order,
                               manifest$search_case_id), , drop = FALSE]
    rownames(manifest) <- NULL

    cols <- c("manifest_version", "search_case_id", "manifest_role",
              "simulation_status", "expected_package_function",
              "package_function", "fixture_set_id", "family", "mode",
              "bf_type", "bf_prior_id", "design_case_id", "look_grid_name",
              "schedule_id", "schedule_family_id", "start_n", "increment",
              "n_looks", "max_n", "threshold_pair_id",
              "evidence_threshold", "k1", "k0", "evidence", "target_prob",
              "criterion", "achieved", "look_found", "n_found",
              "prob_found", "mcse_found", "look_previous", "n_previous",
              "prob_previous", "best_look", "best_n", "best_prob",
              "max_n_searched", "target_source", "search_profile_id",
              "search_profile_label", "simulation_schedule", "look_count",
              "grid_step", "reference_n", "reference_schedule",
              "reference_previous_n", "reference_previous_prob",
              "reference_growth", "schedule_start_n", "schedule_increment",
              "sim_achieved", "simulation_n", "simulation_prob",
              "simulation_mcse", "simulation_previous_n",
              "simulation_previous_prob", "simulation_best_n",
              "simulation_best_prob", "simulation_max_n_searched",
              "package_nrange_lower", "package_nrange_upper",
              "selection_stratum", "package_key", "source_cumulative_summary")
    manifest[cols[cols %in% names(manifest)]]
}

bfpwr_sim_manifest_to_search_simulation <- function(manifest) {
    cols <- c("fixture_set_id", "family", "mode", "bf_type", "bf_prior_id",
              "design_case_id", "look_grid_name", "schedule_id",
              "schedule_family_id", "start_n", "increment", "n_looks",
              "max_n", "search_profile_id", "search_profile_label",
              "simulation_schedule", "look_count", "grid_step",
              "target_source", "reference_n", "reference_schedule",
              "reference_previous_n", "reference_previous_prob",
              "reference_growth", "schedule_start_n",
              "schedule_increment", "evidence_threshold", "evidence",
              "target_prob", "sim_achieved", "simulation_n",
              "simulation_prob", "simulation_mcse",
              "simulation_previous_n", "simulation_previous_prob",
              "simulation_best_n", "simulation_best_prob",
              "simulation_max_n_searched", "package_nrange_lower",
              "package_nrange_upper", "package_key", "manifest_version",
              "search_case_id", "achieved", "prob_found", "mcse_found",
              "best_n", "best_prob", "max_n_searched",
              "manifest_role", "simulation_status",
              "expected_package_function", "selection_stratum")
    missing <- setdiff(cols, names(manifest))
    for (col in missing) {
        manifest[[col]] <- NA
    }
    manifest[cols]
}
