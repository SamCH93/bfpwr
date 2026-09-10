bfpwr_sim_package_verification_manifest_file <- function(repo_root = ".") {
    file.path(repo_root, "simulations", "registry", "manifests",
              "package-verification-cases.csv")
}

bfpwr_sim_package_verification_coverage_file <- function(repo_root = ".") {
    file.path(repo_root, "simulations", "registry", "manifests",
              "package-verification-cases-coverage.csv")
}

bfpwr_sim_reference_package_function <- function(family, bf_type, mode) {
    if (identical(mode, "sequential")) {
        if (identical(family, "z")) return("pbf01seq")
        if (identical(family, "t")) return("ptbf01seq")
        return(NA_character_)
    }
    if (identical(family, "z")) {
        if (identical(bf_type, "moment")) return("pnmbf01")
        return("pbf01")
    }
    if (identical(family, "t")) return("ptbf01")
    if (identical(family, "binomial")) return("pbinbf01")
    NA_character_
}

bfpwr_sim_fixture_spec_table <- function(corpus_root) {
    specs <- list.files(file.path(corpus_root, "fixtures"),
                        pattern = "^spec\\.rds$", recursive = TRUE,
                        full.names = TRUE)
    if (length(specs) == 0) {
        stop("no fixture specs found under ", file.path(corpus_root, "fixtures"),
             call. = FALSE)
    }
    rows <- lapply(specs, function(path) {
        spec <- readRDS(path)
        data.frame(
            fixture_set_id = spec$fixture_set_id,
            family = spec$family,
            mode = spec$mode,
            path = dirname(path),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out[order(out$family, out$mode, out$fixture_set_id), , drop = FALSE]
}

bfpwr_sim_package_manifest_case_id <- function(prefix, rows) {
    sprintf("%s-%04d", prefix, seq_len(nrow(rows)))
}

bfpwr_sim_fixed_package_reference_manifest <- function(fixture) {
    ## The package-verification manifest is the runtime contract for package
    ## calls. Keep fixed-n cases as explicit cells so report builds do not
    ## silently expand to the full Monte Carlo corpus.
    rows <- bfpwr_sim_select_fixed_reference_plot_rows(fixture$tail_summary)
    keep_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                   "design_case_id", "look_grid_name", "threshold_id",
                   "evidence_threshold", "tail", "threshold",
                   "log_threshold", "n", "nsim")
    rows <- rows[intersect(keep_cols, names(rows))]
    if (nrow(rows) == 0) {
        return(data.frame())
    }
    rows <- rows[!duplicated(bfpwr_sim_row_key(
        rows,
        c("fixture_set_id", "bf_prior_id", "threshold_id", "tail", "n"))),
        , drop = FALSE]
    rows <- rows[order(rows$family, rows$bf_type, rows$design_case_id,
                       rows$bf_prior_id, rows$threshold_id, rows$tail,
                       rows$n), , drop = FALSE]
    rows$mode <- "fixed"
    rows$validation_role <- "fixed_probability_reference"
    rows$manifest_role <- rows$validation_role
    rows$package_function <- vapply(
        seq_len(nrow(rows)),
        function(i) {
            bfpwr_sim_reference_package_function(rows$family[[i]],
                                                 rows$bf_type[[i]], "fixed")
        },
        character(1))
    rows$reference_scope <- "selected_fixed_tail_cells"
    rows$validation_required <- TRUE
    rows$rationale <- paste(
        "Vectorized fixed-n package probability reference for a deterministic",
        "cell selection spanning each simulated prior, threshold, tail, and",
        "sample-size range.")
    rows
}

bfpwr_sim_sequential_design_class <- function(design_case_id) {
    x <- as.character(design_case_id)
    out <- rep("other", length(x))
    out[grepl("usd50|s20", x)] <- "extreme_scale"
    out[grepl("m0p|dpoint-m", x)] <- "wrong_direction"
    out[grepl("n2x1p1", x)] <- "unequal_two_sample"
    out[grepl("paired", x)] <- "paired"
    out[grepl("one-", x)] <- "one_sample"
    out[grepl("two-dnorm", x)] <- "two_sample_uncertain"
    out[grepl("two-dpoint", x)] <- "two_sample_point"
    out[grepl("dnorm-0-s|dpoint-0-", x)] <- "null"
    out[grepl("0p2|0p35", x) & out == "other"] <- "local_effect"
    out[grepl("0p5|10", x) & out == "other"] <- "moderate_or_large"
    out
}

bfpwr_sim_sequential_prior_class <- function(bf_prior_id) {
    x <- as.character(bf_prior_id)
    out <- rep("other_prior", length(x))
    out[grepl("point", x)] <- "point"
    out[grepl("normal", x)] <- "normal"
    out[grepl("moment", x)] <- "moment"
    out[grepl("directional", x)] <- "directional"
    out[grepl("cauchy", x)] <- "cauchy"
    out[grepl("student-t", x)] <- "student_t"
    out[grepl("scale0p35|psd0p353", x)] <- paste(out[grepl("scale0p35|psd0p353", x)],
                                                  "narrow", sep = "_")
    out[grepl("scale1|psd20", x)] <- paste(out[grepl("scale1|psd20", x)],
                                            "wide", sep = "_")
    out[grepl("df30", x)] <- paste(out[grepl("df30", x)], "high_df",
                                    sep = "_")
    out[grepl("null0p2", x)] <- paste(out[grepl("null0p2", x)],
                                       "shifted_null", sep = "_")
    out
}

bfpwr_sim_sequential_direction_class <- function(bf_prior_id) {
    x <- as.character(bf_prior_id)
    out <- rep("two_sided_or_symmetric", length(x))
    out[grepl("greater", x)] <- "greater"
    out[grepl("less", x)] <- "less"
    out[grepl("directional", x)] <- "directional"
    out
}

bfpwr_sim_add_sequential_selection_fields <- function(rows) {
    if (nrow(rows) == 0) return(rows)
    rows$design_class <- bfpwr_sim_sequential_design_class(
        rows$design_case_id)
    rows$prior_class <- bfpwr_sim_sequential_prior_class(rows$bf_prior_id)
    rows$direction_class <- bfpwr_sim_sequential_direction_class(
        rows$bf_prior_id)
    rows$schedule_class <- ifelse(
        rows$n_looks <= 10, "short_exact",
        ifelse(rows$n_looks <= 20, "medium_search", "stress_search"))
    rows
}

bfpwr_sim_select_diverse_rows <- function(rows,
                                          n,
                                          strata_cols = c("design_class",
                                                          "prior_class",
                                                          "direction_class",
                                                          "look_grid_name"),
                                          order_cols = c("design_class",
                                                         "prior_class",
                                                         "direction_class",
                                                         "design_case_id",
                                                         "bf_prior_id"),
                                          offset = 0L) {
    if (nrow(rows) <= n) {
        return(rows)
    }
    rows <- bfpwr_sim_add_sequential_selection_fields(rows)
    rows$.selection_id <- seq_len(nrow(rows))
    strata_cols <- strata_cols[strata_cols %in% names(rows)]
    rows$.selection_stratum <- bfpwr_sim_row_key(rows, strata_cols)
    groups <- split(seq_len(nrow(rows)), rows$.selection_stratum)
    group_names <- sort(names(groups))
    if (length(group_names) == 0) {
        return(utils::head(rows, n))
    }
    if (offset > 0L && length(group_names) > 1L) {
        shift <- offset %% length(group_names)
        if (shift > 0L) {
            group_names <- c(group_names[(shift + 1L):length(group_names)],
                             group_names[seq_len(shift)])
        }
    }
    order_cols <- order_cols[order_cols %in% names(rows)]
    selected <- integer()
    cursor <- setNames(rep(1L, length(groups)), names(groups))
    while (length(selected) < n) {
        progressed <- FALSE
        for (name in group_names) {
            group <- rows[groups[[name]], , drop = FALSE]
            if (length(order_cols) > 0) {
                ord <- do.call(order,
                                c(group[order_cols], list(na.last = TRUE)))
                group <- group[ord, , drop = FALSE]
            }
            j <- cursor[[name]]
            if (j <= nrow(group)) {
                selected <- c(selected, group$.selection_id[[j]])
                cursor[[name]] <- j + 1L
                progressed <- TRUE
                if (length(selected) == n) break
            }
        }
        if (!progressed) break
    }
    out <- rows[selected, , drop = FALSE]
    out$.selection_stratum <- NULL
    out$.selection_id <- NULL
    rownames(out) <- NULL
    out
}

bfpwr_sim_select_one_diverse_row <- function(rows, offset = 0L) {
    rows <- bfpwr_sim_add_sequential_selection_fields(rows)
    strata_cols <- c("design_class", "prior_class", "direction_class",
                     "look_grid_name")
    strata_cols <- strata_cols[strata_cols %in% names(rows)]
    rows$.selection_stratum <- bfpwr_sim_row_key(rows, strata_cols)
    strata <- sort(unique(rows$.selection_stratum))
    if (length(strata) == 0) {
        return(rows[0, , drop = FALSE])
    }
    selected_stratum <- strata[(offset %% length(strata)) + 1L]
    group <- rows[rows$.selection_stratum == selected_stratum, , drop = FALSE]
    order_cols <- c("design_case_id", "bf_prior_id", "schedule_id",
                    "evidence_threshold", "evidence", "target_prob")
    order_cols <- order_cols[order_cols %in% names(group)]
    if (length(order_cols) > 0) {
        ord <- do.call(order, c(group[order_cols], list(na.last = TRUE)))
        group <- group[ord, , drop = FALSE]
    }
    out <- group[((offset %/% length(strata)) %% nrow(group)) + 1L,
                 , drop = FALSE]
    out$.selection_stratum <- NULL
    rownames(out) <- NULL
    out
}

bfpwr_sim_sequential_reference_manifest <- function(fixture) {
    if (!identical(fixture$spec$mode, "sequential") ||
        !fixture$spec$family %in% c("z", "t")) {
        return(data.frame())
    }
    short_schedules <- sprintf("start20-by10-looks%02d", c(2, 5, 10))
    rows <- fixture$final_summary[
        fixture$final_summary$schedule_id %in% short_schedules &
            fixture$final_summary$evidence_threshold %in% c(3, 10, 30),
        , drop = FALSE]
    if (nrow(rows) == 0) {
        return(data.frame())
    }
    key_cols <- c("fixture_set_id", "bf_prior_id", "schedule_id",
                  "evidence_threshold")
    rows <- rows[!duplicated(bfpwr_sim_row_key(rows, key_cols)),
                 , drop = FALSE]
    split_key <- bfpwr_sim_row_key(rows, c("schedule_id",
                                           "evidence_threshold"))
    groups <- split(rows, split_key)
    per_cell <- if (identical(fixture$spec$family, "z")) 8L else 16L
    selected <- vector("list", length(groups))
    group_names <- sort(names(groups))
    for (i in seq_along(group_names)) {
        group <- groups[[group_names[[i]]]]
        selected[[i]] <- bfpwr_sim_select_diverse_rows(
            group, n = per_cell, offset = i - 1L)
    }
    rows <- bfpwr_sim_bind_rows_fill(selected)
    rows <- bfpwr_sim_add_sequential_selection_fields(rows)
    keep_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                   "design_case_id", "look_grid_name", "schedule_id",
                   "schedule_family_id", "start_n", "increment", "n_looks",
                   "max_n", "threshold_pair_id", "evidence_threshold",
                   "k1", "k0", "nsim", "design_class", "prior_class",
                   "direction_class", "schedule_class")
    rows <- rows[intersect(keep_cols, names(rows))]
    rows$mode <- "sequential"
    rows$validation_role <- if (identical(fixture$spec$family, "t")) {
        "sequential_expected_sample_size_reference"
    } else {
        "sequential_package_reference"
    }
    rows$manifest_role <- rows$validation_role
    rows$package_function <- if (identical(fixture$spec$family, "t")) {
        "ptbf01seq"
    } else {
        "pbf01seq"
    }
    rows$reference_scope <- if (identical(fixture$spec$family, "t")) {
        "short_strict_expected_sample_size_only"
    } else {
        "short_strict_all_tails_all_looks_and_expected_sample_size"
    }
    rows$validation_required <- TRUE
    rows$package_key <- NA_character_
    rows$selection_stratum <- paste(rows$schedule_class, rows$design_class,
                                    rows$prior_class, rows$direction_class,
                                    sep = ":")
    rows$rationale <- paste(
        "Explicit short-schedule sequential package reference selected from",
        "balanced design/prior strata; long schedules are covered by",
        "sequential n-search validation.")
    rows
}

bfpwr_sim_sequential_search_explicit_grid <- function(package_function) {
    ## Package verification should exercise a fixed policy surface, not a
    ## replenished quota. Cells that do not robustly cross remain diagnostics.
    primary <- expand.grid(
        schedule_id = sprintf("start20-by10-looks%02d", c(2, 5, 10, 20)),
        evidence_threshold = c(3, 10, 30),
        evidence = c("H1", "H0"),
        target_prob = c(0.10, 0.30, 0.80, 0.95),
        grid_role = "primary_grid",
        stringsAsFactors = FALSE)
    stress <- expand.grid(
        schedule_id = c("start10-by10-looks50",
                        "start10-by10-looks100",
                        "start10-by10-looks200"),
        evidence_threshold = c(3, 10, 30),
        evidence = c("H1", "H0"),
        target_prob = c(0.30, 0.80),
        grid_role = "stress_grid",
        stringsAsFactors = FALSE)
    out <- rbind(primary, stress)
    out$package_function <- package_function
    out
}

bfpwr_sim_search_robust_achieved <- function(rows,
                                             required_mcse_margin = 2) {
    rows$achieved %in% TRUE &
        is.finite(rows$n_found) &
        is.finite(rows$prob_found) &
        is.finite(rows$mcse_found) &
        rows$prob_found - required_mcse_margin * rows$mcse_found >=
            rows$target_prob
}

bfpwr_sim_add_package_search_rows <- function(rows,
                                              manifest_role,
                                              validation_required,
                                              selection_origin) {
    if (nrow(rows) == 0) return(rows)
    rows <- bfpwr_sim_add_search_selection_fields(rows, manifest_role)
    rows$validation_required <- validation_required
    rows$validation_role <- "sequential_n_search"
    rows$reference_scope <- "explicit_first_crossing_sample_size_search"
    if (!"selection_origin" %in% names(rows)) {
        rows$selection_origin <- selection_origin
    }
    rows$selection_policy <- selection_origin
    rows$rationale <- paste(
        "Explicit sequential sample-size search case selected from schedule,",
        "threshold, evidence-direction, target-probability, and design/prior",
        "strata.")
    rows
}

bfpwr_sim_package_sequential_search_manifest <- function(
        corpus_root,
        families = c("z", "t"),
        required_mcse_margin = 2,
        diagnostics_per_function = 12L) {
    all_rows <- bfpwr_sim_read_sequential_search_summaries(
        corpus_root = corpus_root, families = families)
    if (nrow(all_rows) == 0) {
        return(data.frame())
    }
    all_rows <- bfpwr_sim_add_sequential_selection_fields(all_rows)
    selected <- list()
    idx <- 0L
    functions <- sort(unique(all_rows$package_function))
    for (fn in functions) {
        rows <- all_rows[all_rows$package_function == fn, , drop = FALSE]
        grid <- bfpwr_sim_sequential_search_explicit_grid(fn)
        grid_rows <- list()
        grid_idx <- 0L
        for (i in seq_len(nrow(grid))) {
            g <- grid[i, , drop = FALSE]
            candidates <- rows[
                rows$schedule_id == g$schedule_id[[1]] &
                    rows$evidence_threshold == g$evidence_threshold[[1]] &
                    rows$evidence == g$evidence[[1]] &
                    rows$target_prob == g$target_prob[[1]],
                , drop = FALSE]
            if (nrow(candidates) == 0) {
                next
            }
            grid_idx <- grid_idx + 1L
            one <- bfpwr_sim_select_one_diverse_row(
                candidates, offset = i + grid_idx)
            one$selection_origin <- g$grid_role[[1]]
            grid_rows[[grid_idx]] <- one
        }
        grid_rows <- bfpwr_sim_bind_rows_fill(grid_rows)
        if (nrow(grid_rows) > 0) {
            robust <- bfpwr_sim_search_robust_achieved(
                grid_rows, required_mcse_margin = required_mcse_margin)
            if (any(robust)) {
                idx <- idx + 1L
                selected[[idx]] <- bfpwr_sim_add_package_search_rows(
                    grid_rows[robust, , drop = FALSE],
                    manifest_role = "required_n_search",
                    validation_required = TRUE,
                    selection_origin = "explicit_required_grid")
            }
            if (any(!robust)) {
                idx <- idx + 1L
                selected[[idx]] <- bfpwr_sim_add_package_search_rows(
                    grid_rows[!robust, , drop = FALSE],
                    manifest_role = "boundary_diagnostic",
                    validation_required = FALSE,
                    selection_origin = "explicit_nonrobust_grid")
            }
        }

        selected_key <- if (nrow(grid_rows) > 0) {
            bfpwr_sim_search_validation_key(grid_rows)
        } else character()
        diagnostics <- rows[!(rows$achieved %in% TRUE), , drop = FALSE]
        if (nrow(diagnostics) > 0 && diagnostics_per_function > 0) {
            diagnostics <- diagnostics[
                !bfpwr_sim_search_validation_key(diagnostics) %in%
                    selected_key,
                , drop = FALSE]
            if (nrow(diagnostics) > 0) {
                idx <- idx + 1L
                selected[[idx]] <- bfpwr_sim_add_package_search_rows(
                    bfpwr_sim_select_diverse_rows(
                        diagnostics,
                        n = as.integer(diagnostics_per_function),
                        strata_cols = c("schedule_class", "evidence",
                                        "evidence_threshold", "target_prob",
                                        "design_class", "prior_class"),
                        offset = length(selected)),
                    manifest_role = "boundary_diagnostic",
                    validation_required = FALSE,
                    selection_origin = "explicit_nonachieved_diagnostic")
            }
        }
    }
    rows <- bfpwr_sim_bind_rows_fill(selected)
    if (nrow(rows) == 0) {
        return(rows)
    }
    rows$search_case_id <- sprintf("pkgseq-search-%04d", seq_len(nrow(rows)))
    rows$package_function <- rows$expected_package_function
    rows
}

bfpwr_sim_sequential_reference_manifest_from_search <- function(search_manifest) {
    if (is.null(search_manifest) || nrow(search_manifest) == 0) {
        return(data.frame())
    }
    ## Retained only for manually supplied legacy manifests. The default
    ## package-verification manifest now derives sequential references from
    ## final_summary and n-search rows from cumulative_summary.
    keep <- search_manifest$manifest_role == "required_n_search" &
        search_manifest$family %in% c("z", "t") &
        search_manifest$mode == "sequential"
    rows <- search_manifest[keep, , drop = FALSE]
    if (nrow(rows) == 0) {
        return(data.frame())
    }
    if ("n_looks" %in% names(rows)) {
        rows <- rows[rows$n_looks <= 10, , drop = FALSE]
    }
    if (nrow(rows) == 0) {
        return(data.frame())
    }
    rows$validation_role <- ifelse(
        rows$family == "t",
        "sequential_expected_sample_size_reference",
        "sequential_package_reference")
    rows$manifest_role <- rows$validation_role
    rows$package_function <- ifelse(rows$family == "t", "ptbf01seq",
                                    "pbf01seq")
    rows$reference_scope <- ifelse(
        rows$family == "t",
        "legacy_short_expected_sample_size_only",
        "legacy_short_all_tails_all_looks_and_expected_sample_size")
    rows$validation_required <- TRUE
    rows$source_search_case_id <- if ("search_case_id" %in% names(rows)) {
        rows$search_case_id
    } else NA_character_
    rows$package_key <- NA_character_
    rows$rationale <- paste(
        "Legacy sequential package reference derived from a supplied",
        "sequential search manifest.")
    rows
}

bfpwr_sim_sequential_integrity_manifest <- function(fixture) {
    if (!identical(fixture$spec$mode, "sequential") ||
        !identical(fixture$spec$family, "binomial")) {
        return(data.frame())
    }
    data.frame(
        fixture_set_id = fixture$spec$fixture_set_id,
        family = fixture$spec$family,
        mode = fixture$spec$mode,
        bf_type = fixture$spec$bf_type,
        bf_prior_id = NA_character_,
        design_case_id = NA_character_,
        look_grid_name = NA_character_,
        validation_role = "fixture_integrity_only",
        manifest_role = "fixture_integrity_only",
        package_function = NA_character_,
        reference_scope = "no_sequential_binomial_package_api",
        validation_required = TRUE,
        rationale = paste(
            "Sequential binomial fixtures are retained for manifest and",
            "deterministic integrity checks; the package has no sequential",
            "binomial probability API."),
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_read_package_verification_manifest <- function(path,
                                                        required = FALSE) {
    if (is.null(path) || !nzchar(path)) {
        if (isTRUE(required)) {
            stop("package verification manifest path is empty", call. = FALSE)
        }
        return(NULL)
    }
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!file.exists(path)) {
        if (isTRUE(required)) {
            stop("package verification manifest is missing: ", path,
                 call. = FALSE)
        }
        return(NULL)
    }
    manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
    required_cols <- c("package_verification_case_id", "validation_role",
                       "fixture_set_id", "family", "mode")
    missing <- setdiff(required_cols, names(manifest))
    if (length(missing) > 0) {
        stop("package verification manifest is missing columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    manifest
}

bfpwr_sim_build_package_verification_manifest <- function(
        corpus_root,
        search_manifest = NULL) {
    fixture_specs <- bfpwr_sim_fixture_spec_table(corpus_root)
    rows <- list()
    idx <- 0L
    for (i in seq_len(nrow(fixture_specs))) {
        spec <- fixture_specs[i, , drop = FALSE]
        fixture <- bfpwr_sim_read_fixture_summary(
            corpus_root = corpus_root,
            fixture_set_id = spec$fixture_set_id[[1]],
            family = spec$family[[1]],
            mode = spec$mode[[1]])
        if (identical(spec$mode[[1]], "fixed")) {
            idx <- idx + 1L
            rows[[idx]] <- bfpwr_sim_fixed_package_reference_manifest(fixture)
        } else if (identical(spec$family[[1]], "binomial")) {
            idx <- idx + 1L
            rows[[idx]] <- bfpwr_sim_sequential_integrity_manifest(fixture)
        } else if (spec$family[[1]] %in% c("z", "t")) {
            idx <- idx + 1L
            rows[[idx]] <- bfpwr_sim_sequential_reference_manifest(fixture)
        }
    }

    search_rows <- bfpwr_sim_package_sequential_search_manifest(
        corpus_root = corpus_root)
    if (nrow(search_rows) > 0) {
        idx <- idx + 1L
        rows[[idx]] <- search_rows
    }

    if (!is.null(search_manifest) && nrow(search_manifest) > 0) {
        ## Manual legacy append only. The default package-verification
        ## contract is now rebuilt directly from sequential fixture summaries
        ## so it can use explicit schedule/threshold/target strata.
        seq_ref <- bfpwr_sim_sequential_reference_manifest_from_search(
            search_manifest)
        if (nrow(seq_ref) > 0) {
            idx <- idx + 1L
            rows[[idx]] <- seq_ref
        }
    }

    manifest <- bfpwr_sim_bind_rows_fill(rows)
    if (nrow(manifest) == 0) {
        return(manifest)
    }
    manifest$manifest_version <- "package-verification-v2"
    if (!"package_verification_case_id" %in% names(manifest)) {
        manifest$package_verification_case_id <- NA_character_
    }
    missing_case_id <- is.na(manifest$package_verification_case_id) |
        !nzchar(manifest$package_verification_case_id)
    if (any(missing_case_id)) {
        manifest$package_verification_case_id[missing_case_id] <-
            bfpwr_sim_package_manifest_case_id(
                "pkgver", manifest[missing_case_id, , drop = FALSE])
    }
    if (!"package_key" %in% names(manifest)) {
        manifest$package_key <- NA_character_
    }
    no_package_key <- is.na(manifest$package_key) | !nzchar(manifest$package_key)
    if (any(no_package_key)) {
        key_cols <- c("manifest_version", "validation_role", "fixture_set_id",
                      "family", "mode", "bf_prior_id", "design_case_id",
                      "threshold_id", "tail", "n", "schedule_id",
                      "evidence_threshold",
                      "package_verification_case_id")
        key_cols <- key_cols[key_cols %in% names(manifest)]
        key_source <- manifest[no_package_key, key_cols, drop = FALSE]
        manifest$package_key[no_package_key] <- apply(
            key_source, 1L, function(x) {
                paste(ifelse(is.na(x), "", as.character(x)), collapse = "|")
            })
    }

    preferred <- c("manifest_version", "package_verification_case_id",
                   "validation_role", "manifest_role", "validation_required",
                   "reference_scope", "package_function", "package_key",
                   "fixture_set_id", "family", "mode", "bf_type",
                   "bf_prior_id", "design_case_id", "look_grid_name",
                   "design_class", "prior_class", "direction_class",
                   "threshold_id", "tail", "n", "threshold",
                   "log_threshold",
                   "schedule_id", "schedule_family_id", "start_n",
                   "increment", "n_looks", "max_n", "threshold_pair_id",
                   "evidence_threshold", "k1", "k0", "evidence",
                   "target_prob", "simulation_status", "source_search_case_id",
                   "search_case_id", "schedule_class", "selection_origin",
                   "selection_policy", "selection_stratum", "rationale")
    cols <- c(preferred[preferred %in% names(manifest)],
              setdiff(names(manifest), preferred))
    manifest <- manifest[cols]
    manifest <- manifest[order(manifest$validation_role,
                               manifest$family, manifest$mode,
                               manifest$fixture_set_id,
                               manifest$package_verification_case_id),
                         , drop = FALSE]
    rownames(manifest) <- NULL
    manifest
}

bfpwr_sim_package_verification_coverage <- function(manifest) {
    if (is.null(manifest) || nrow(manifest) == 0) {
        return(data.frame())
    }
    package_function <- manifest$package_function
    package_function[is.na(package_function) | !nzchar(package_function)] <-
        "(none)"
    split_key <- interaction(manifest$validation_role, manifest$family,
                             manifest$mode, package_function,
                             drop = TRUE)
    rows <- lapply(split(manifest, split_key), function(x) {
        fn <- x$package_function[[1]]
        if (is.na(fn) || !nzchar(fn)) fn <- "(none)"
        data.frame(
            validation_role = x$validation_role[[1]],
            family = x$family[[1]],
            mode = x$mode[[1]],
            package_function = fn,
            rows = nrow(x),
            required_rows = sum(x$validation_required %in% TRUE,
                                na.rm = TRUE),
            fixtures = length(unique(x$fixture_set_id)),
            bf_priors = if ("bf_prior_id" %in% names(x)) {
                keep <- !is.na(x$bf_prior_id) & nzchar(x$bf_prior_id)
                length(unique(x$bf_prior_id[keep]))
            } else 0L,
            designs = if ("design_case_id" %in% names(x)) {
                keep <- !is.na(x$design_case_id) &
                    nzchar(x$design_case_id)
                length(unique(x$design_case_id[keep]))
            } else 0L,
            schedules = if ("schedule_id" %in% names(x)) {
                keep <- !is.na(x$schedule_id) & nzchar(x$schedule_id)
                length(unique(x$schedule_id[keep]))
            } else 0L,
            thresholds = if ("evidence_threshold" %in% names(x)) {
                paste(sort(unique(x$evidence_threshold[
                    is.finite(x$evidence_threshold)])), collapse = ", ")
            } else "",
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    out <- out[order(out$validation_role, out$family, out$mode,
                     out$package_function), , drop = FALSE]
    rownames(out) <- NULL
    out
}
