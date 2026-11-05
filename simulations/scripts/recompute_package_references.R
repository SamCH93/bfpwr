support_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    script <- if (length(hit) > 0) {
        normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                      winslash = "/", mustWork = TRUE)
    } else {
        normalizePath("simulations/scripts/recompute_package_references.R",
                      winslash = "/", mustWork = FALSE)
    }
    file.path(dirname(script), "simulation_support.R")
}

source(support_path())
source_package_checkout()
source_simulation_library()

write_csv <- function(x, path) {
    ensure_dir(dirname(path))
    utils::write.csv(x, path, row.names = FALSE, na = "")
    invisible(path)
}

fixture_specs <- function(corpus_root, fixture_set = NULL) {
    specs <- list.files(file.path(corpus_root, "fixtures"),
                        pattern = "^spec\\.rds$", recursive = TRUE,
                        full.names = TRUE)
    if (length(specs) == 0) {
        stop("no fixture specs found under ", file.path(corpus_root, "fixtures"),
             call. = FALSE)
    }
    out <- lapply(specs, function(path) {
        spec <- readRDS(path)
        data.frame(
            fixture_set_id = spec$fixture_set_id,
            family = spec$family,
            mode = spec$mode,
            path = dirname(path),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, out)
    if (!is.null(fixture_set) && nzchar(fixture_set)) {
        requested <- trimws(strsplit(fixture_set, ",", fixed = TRUE)[[1]])
        requested <- requested[nzchar(requested)]
        out <- out[out$fixture_set_id %in% requested, , drop = FALSE]
        missing <- setdiff(requested, out$fixture_set_id)
        if (length(missing) > 0) {
            stop("requested fixture set(s) not found: ",
                 paste(missing, collapse = ", "), call. = FALSE)
        }
    }
    out[order(out$family, out$mode, out$fixture_set_id), , drop = FALSE]
}

write_reference_result <- function(reference, prefix) {
    write_csv(reference$summary, paste0(prefix, "-reference-summary.csv"))
    write_csv(reference$group_summary, paste0(prefix, "-reference-groups.csv"))
    write_csv(reference$outliers, paste0(prefix, "-reference-outliers.csv"))
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
        write_csv(reference$diagnostic_rows,
                  paste0(prefix, "-reference-diagnostic-rows.csv"))
    }
    if (!is.null(reference$en_summary)) {
        write_csv(reference$en_summary,
                  paste0(prefix, "-reference-en-summary.csv"))
    }
    if (!is.null(reference$en_rows)) {
        write_csv(reference$en_rows, paste0(prefix, "-reference-en-rows.csv"))
    }
    if (!is.null(reference$timings)) {
        write_csv(reference$timings, paste0(prefix, "-reference-timings.csv"))
    }
}

recompute_fixture_references <- function(corpus_root,
                                          output_dir,
                                          fixture_set = NULL,
                                          case_set = "production",
                                          package_manifest = NULL) {
    specs <- fixture_specs(corpus_root, fixture_set)
    designs <- bfpwr_sim_design_case_set(case_set)
    bf_priors <- bfpwr_sim_bf_prior_case_set(case_set)
    status <- vector("list", nrow(specs))
    status_path <- file.path(output_dir, "package-reference-recompute-status.csv")

    for (i in seq_len(nrow(specs))) {
        spec <- specs[i, , drop = FALSE]
        prefix <- file.path(output_dir, spec$fixture_set_id)
        role <- if (identical(spec$family[[1]], "binomial") &&
                    identical(spec$mode[[1]], "sequential")) {
            "integrity_only"
        } else {
            "package_reference"
        }

        if (identical(role, "integrity_only")) {
            message("skipping package reference for unsupported fixture ",
                    spec$fixture_set_id)
            status[[i]] <- data.frame(
                fixture_set_id = spec$fixture_set_id,
                family = spec$family,
                mode = spec$mode,
                action = "skipped",
                reason = "sequential binomial has no package probability API",
                stringsAsFactors = FALSE
            )
            write_csv(do.call(rbind, status[seq_len(i)]), status_path)
            next
        }

        message("recomputing fixture package reference ",
                spec$fixture_set_id, " (", i, " of ", nrow(specs), ")")
        t0 <- proc.time()[["elapsed"]]
        fixture <- bfpwr_sim_read_fixture_summary(
            corpus_root = corpus_root,
            fixture_set_id = spec$fixture_set_id,
            family = spec$family,
            mode = spec$mode
        )
        manifest_cases <- if (!is.null(package_manifest) &&
                              nrow(package_manifest) > 0) {
            package_manifest[
                package_manifest$fixture_set_id == spec$fixture_set_id &
                    package_manifest$family == spec$family &
                    package_manifest$mode == spec$mode,
                , drop = FALSE]
        } else {
            NULL
        }
        reference <- if (identical(spec$mode[[1]], "fixed")) {
            bfpwr_sim_validate_fixed_references(
                fixture, designs = designs, bf_priors = bf_priors,
                manifest_cases = manifest_cases)
        } else if (spec$family[[1]] %in% c("z", "t")) {
            bfpwr_sim_validate_sequential_references(
                fixture, designs = designs, bf_priors = bf_priors,
                profile = "curated", strict = TRUE,
                manifest_cases = manifest_cases)
        } else {
            NULL
        }

        if (is.null(reference)) {
            status[[i]] <- data.frame(
                fixture_set_id = spec$fixture_set_id,
                family = spec$family,
                mode = spec$mode,
                action = "skipped",
                reason = "no package-reference adapter",
                stringsAsFactors = FALSE
            )
        } else {
            write_reference_result(reference, prefix)
            status[[i]] <- data.frame(
                fixture_set_id = spec$fixture_set_id,
                family = spec$family,
                mode = spec$mode,
                action = "recomputed",
                reason = "",
                stringsAsFactors = FALSE
            )
        }
        message("finished ", spec$fixture_set_id, " in ",
                round(proc.time()[["elapsed"]] - t0, 1), " seconds")
        write_csv(do.call(rbind, status[seq_len(i)]), status_path)
    }

    status <- do.call(rbind, status)
    write_csv(status, status_path)
    invisible(status)
}

row_value <- function(row, name, default = NA) {
    if (name %in% names(row)) row[[name]][[1]] else default
}

capture_package_eval <- function(expr) {
    warnings <- character()
    t0 <- proc.time()[["elapsed"]]
    value <- tryCatch(
        withCallingHandlers(
            expr,
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )
    list(value = value,
         warnings = unique(warnings),
         elapsed_seconds = proc.time()[["elapsed"]] - t0)
}

empty_package_row <- function(package_key) {
    n <- length(package_key)
    data.frame(
        package_key = package_key,
        package_function = rep(NA_character_, n),
        package_search = rep(NA_character_, n),
        package_reached = rep(NA, n),
        package_n = rep(NA_integer_, n),
        package_maximum_n = rep(NA_integer_, n),
        package_actual_prob = rep(NA_real_, n),
        package_reference_prob = rep(NA_real_, n),
        package_evaluations = rep(NA_integer_, n),
        package_first_crossing_certified = rep(NA, n),
        package_elapsed_seconds = rep(NA_real_, n),
        package_error = rep(NA_character_, n),
        package_warnings = rep(NA_character_, n),
        stringsAsFactors = FALSE
    )
}

finite_scalar <- function(x) {
    length(x) == 1L && is.finite(x)
}

target_tail <- function(row) {
    if (identical(row$evidence[[1]], "H1")) "H1" else "H0"
}

threshold_for_row <- function(row) {
    if (identical(row$evidence[[1]], "H1")) {
        1 / row$evidence_threshold[[1]]
    } else {
        row$evidence_threshold[[1]]
    }
}

lower_tail_for_row <- function(row) {
    identical(row$evidence[[1]], "H1")
}

search_schedule_args <- function(row) {
    if (identical(row$mode[[1]], "fixed")) {
        return(list(looks = 1L))
    }
    start_n <- suppressWarnings(as.integer(row_value(row, "schedule_start_n")))
    increment <- suppressWarnings(as.integer(row_value(row,
                                                       "schedule_increment")))
    if (is.finite(start_n) && is.finite(increment)) {
        return(list(minN = start_n, by = increment))
    }
    list(looks = as.integer(row$look_count[[1]]))
}

search_nrange <- function(row) {
    c(as.integer(row$package_nrange_lower[[1]]),
      as.integer(row$package_nrange_upper[[1]]))
}

solver_to_package_row <- function(package_key,
                                  package_function,
                                  package_search,
                                  solver,
                                  warnings,
                                  elapsed_seconds = NA_real_) {
    n <- if (is.null(solver$n)) NA_real_ else solver$n
    maximum_n <- if (is.null(solver$maximumN)) n else solver$maximumN
    reached <- isTRUE(solver$reached)
    data.frame(
        package_key = package_key,
        package_function = package_function,
        package_search = package_search,
        package_reached = reached,
        package_n = if (finite_scalar(n)) as.integer(n) else NA_integer_,
        package_maximum_n = if (finite_scalar(maximum_n)) {
            as.integer(maximum_n)
        } else NA_integer_,
        package_actual_prob = if (finite_scalar(solver$actualPower)) {
            solver$actualPower
        } else NA_real_,
        package_reference_prob = NA_real_,
        package_evaluations = if (is.null(solver$evaluations)) {
            NA_integer_
        } else as.integer(solver$evaluations),
        package_first_crossing_certified =
            if (is.null(solver$firstCrossingCertified)) {
                NA
            } else solver$firstCrossingCertified,
        package_elapsed_seconds = elapsed_seconds,
        package_error = if (is.null(solver$error)) "" else solver$error,
        package_warnings = paste(warnings, collapse = " | "),
        stringsAsFactors = FALSE
    )
}

fixed_t_actual <- function(n, row, prior, design, dp) {
    if (!is.finite(n)) {
        return(NA_real_)
    }
    n1 <- as.integer(n)
    n2 <- if (identical(design$generation$type, "two.sample")) {
        pmax(2L, as.integer(round(n1 * design$generation$n2_multiplier)))
    } else {
        n1
    }
    ptbf01(
        k = threshold_for_row(row),
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
        lower.tail = lower_tail_for_row(row)
    )
}

eval_fixed_t_search <- function(row, bf_prior, design, metadata) {
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    package_key <- row$package_key[[1]]
    evaluated <- capture_package_eval({
        ntbf01(
            k = threshold_for_row(row),
            power = row$target_prob[[1]],
            null = prior$null,
            plocation = prior$plocation,
            pscale = prior$pscale,
            pdf = prior$pdf,
            type = prior$type,
            alternative = prior$alternative,
            dpm = dp$dpm,
            dpsd = dp$dpsd,
            lower.tail = lower_tail_for_row(row),
            integer = TRUE,
            nrange = search_nrange(row)
        )
    })
    if (inherits(evaluated$value, "condition")) {
        out <- empty_package_row(package_key)
        out$package_function <- "ntbf01"
        out$package_search <- "fixed_wrapper"
        out$package_elapsed_seconds <- evaluated$elapsed_seconds
        out$package_error <- conditionMessage(evaluated$value)
        out$package_warnings <- paste(evaluated$warnings, collapse = " | ")
        return(out)
    }
    n <- as.numeric(evaluated$value)
    actual <- fixed_t_actual(n, row, prior, design, dp)
    data.frame(
        package_key = package_key,
        package_function = "ntbf01",
        package_search = "fixed_wrapper",
        package_reached = is.finite(n) && is.finite(actual) &&
            actual + 1e-10 >= row$target_prob[[1]],
        package_n = if (is.finite(n)) as.integer(n) else NA_integer_,
        package_maximum_n = if (is.finite(n)) as.integer(n) else NA_integer_,
        package_actual_prob = actual,
        package_reference_prob = NA_real_,
        package_evaluations = NA_integer_,
        package_first_crossing_certified = NA,
        package_elapsed_seconds = evaluated$elapsed_seconds,
        package_error = "",
        package_warnings = paste(evaluated$warnings, collapse = " | "),
        stringsAsFactors = FALSE
    )
}

eval_z_search <- function(row, bf_prior, design, metadata) {
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    schedule <- search_schedule_args(row)
    args <- c(
        list(
            k1 = 1 / row$evidence_threshold[[1]],
            k0 = row$evidence_threshold[[1]],
            power = row$target_prob[[1]],
            usd = design$generation$usd,
            null = prior$null,
            psd = prior$psd,
            dpm = dp$dpm,
            dpsd = dp$dpsd,
            type = bf_prior$bf_type,
            target = target_tail(row),
            nrange = search_nrange(row),
            strict = TRUE,
            integer = TRUE,
            search = "adaptive",
            details = TRUE
        ),
        schedule
    )
    if (!identical(bf_prior$bf_type, "moment")) {
        args$pm <- prior$pm
    }
    evaluated <- capture_package_eval(do.call(nbf01seq, args))
    if (inherits(evaluated$value, "condition")) {
        out <- empty_package_row(row$package_key[[1]])
        out$package_function <- "nbf01seq"
        out$package_search <- "adaptive"
        out$package_elapsed_seconds <- evaluated$elapsed_seconds
        out$package_error <- conditionMessage(evaluated$value)
        out$package_warnings <- paste(evaluated$warnings, collapse = " | ")
        return(out)
    }
    solver_to_package_row(row$package_key[[1]], "nbf01seq", "adaptive",
                          evaluated$value, evaluated$warnings,
                          evaluated$elapsed_seconds)
}

eval_t_seq_search <- function(row, bf_prior, design, metadata) {
    prior <- bf_prior$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    schedule <- search_schedule_args(row)
    ratio <- if (identical(prior$type, "two.sample")) {
        design$generation$n2_multiplier
    } else {
        1
    }
    args <- c(
        list(
            k1 = 1 / row$evidence_threshold[[1]],
            k0 = row$evidence_threshold[[1]],
            power = row$target_prob[[1]],
            plocation = prior$plocation,
            pscale = prior$pscale,
            pdf = prior$pdf,
            dpm = dp$dpm,
            dpsd = dp$dpsd,
            type = prior$type,
            alternative = prior$alternative,
            target = target_tail(row),
            nrange = search_nrange(row),
            ratio = ratio,
            strict = TRUE,
            integer = TRUE,
            search = "adaptive",
            details = TRUE
        ),
        schedule
    )
    evaluated <- capture_package_eval(do.call(ntbf01seq, args))
    if (inherits(evaluated$value, "condition")) {
        out <- empty_package_row(row$package_key[[1]])
        out$package_function <- "ntbf01seq"
        out$package_search <- "adaptive"
        out$package_elapsed_seconds <- evaluated$elapsed_seconds
        out$package_error <- conditionMessage(evaluated$value)
        out$package_warnings <- paste(evaluated$warnings, collapse = " | ")
        return(out)
    }
    solver_to_package_row(row$package_key[[1]], "ntbf01seq", "adaptive",
                          evaluated$value, evaluated$warnings,
                          evaluated$elapsed_seconds)
}

eval_package_search_row <- function(row, designs, bf_priors, metadata) {
    bf_prior <- bfpwr_sim_find_bf_prior_case(row$bf_prior_id[[1]], bf_priors)
    design <- bfpwr_sim_find_design_case(row$design_case_id[[1]], designs)
    if (identical(row$family[[1]], "z")) {
        return(eval_z_search(row, bf_prior, design, metadata))
    }
    if (identical(row$family[[1]], "t") &&
        identical(row$mode[[1]], "fixed")) {
        return(eval_fixed_t_search(row, bf_prior, design, metadata))
    }
    if (identical(row$family[[1]], "t") &&
        identical(row$mode[[1]], "sequential")) {
        return(eval_t_seq_search(row, bf_prior, design, metadata))
    }
    out <- empty_package_row(row$package_key[[1]])
    out$package_error <- "no package search adapter for this row"
    out
}

comparison_with_package <- function(simulation, package) {
    package_cols <- setdiff(names(package), "package_key")
    idx <- match(simulation$package_key, package$package_key)
    package_side <- package[idx, package_cols, drop = FALSE]
    out <- cbind(simulation, package_side)
    out$n_ratio <- out$package_n / out$simulation_n
    out$fractional_error <- (out$package_n - out$simulation_n) /
        out$simulation_n
    out$simulation_target_margin <- out$simulation_prob - out$target_prob
    out$package_target_margin <- out$package_actual_prob - out$target_prob
    out$package_reference_prob_gap <- out$package_actual_prob -
        out$package_reference_prob
    out$package_reference_prob_z <- out$package_reference_prob_gap /
        out$simulation_mcse
    out
}

search_manifest_path <- function(path) {
    if (is.null(path) || !nzchar(path)) {
        return(NULL)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

read_search_manifest <- function(path) {
    path <- search_manifest_path(path)
    if (is.null(path) || !file.exists(path)) {
        return(NULL)
    }
    manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
    if ("validation_role" %in% names(manifest)) {
        manifest <- manifest[
            manifest$validation_role == "sequential_n_search",
            , drop = FALSE]
    }
    required <- c("package_key", "manifest_role", "package_function",
                  "fixture_set_id", "family", "mode", "bf_prior_id",
                  "design_case_id", "schedule_id", "evidence_threshold",
                  "evidence", "target_prob")
    missing <- setdiff(required, names(manifest))
    if (length(missing) > 0) {
        stop("search manifest is missing columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    manifest
}

existing_package_keys <- function(bundle, simulation) {
    keys <- if (!is.null(bundle$package) && nrow(bundle$package) > 0) {
        bundle$package$package_key
    } else if (is.data.frame(bundle$comparison) &&
               all(c("package_key", "package_reached") %in%
                   names(bundle$comparison))) {
        comparison <- bundle$comparison
        comparison$package_key[!is.na(comparison$package_reached)]
    } else {
        character()
    }
    keys <- unique(keys)
    keys[keys %in% simulation$package_key]
}

limit_search_rows <- function(selected, manifest_rows, max_rows) {
    if (!is.finite(max_rows)) {
        return(selected)
    }
    max_rows <- as.integer(max_rows)
    if (max_rows < 0L) {
        stop("--max-search-rows must be non-negative", call. = FALSE)
    }
    if (nrow(manifest_rows) == 0) {
        return(utils::head(selected, max_rows))
    }

    manifest_keys <- unique(manifest_rows$package_key)
    explicit <- selected[selected$package_key %in% manifest_keys,
                         , drop = FALSE]
    legacy <- selected[!selected$package_key %in% manifest_keys,
                       , drop = FALSE]
    legacy_limit <- max_rows - nrow(explicit)
    if (legacy_limit < 0L) {
        message("--max-search-rows is below the explicit manifest size; ",
                "retaining all manifest rows and dropping legacy rows")
        legacy_limit <- 0L
    }
    bfpwr_sim_bind_rows_fill(explicit, utils::head(legacy, legacy_limit))
}

recompute_search_validation <- function(corpus_root,
                                        validation_dir,
                                        max_rows = Inf,
                                        require_search = FALSE,
                                        manifest_file = NULL,
                                        evaluate_diagnostics = FALSE) {
    t0 <- proc.time()[["elapsed"]]
    path <- file.path(validation_dir, "search-validation-comparison.rds")
    manifest <- read_search_manifest(manifest_file)
    manifest_rows <- data.frame()
    manifest_is_package_verification <- FALSE
    manifest_only_bundle <- FALSE
    if (!is.null(manifest)) {
        manifest_is_package_verification <- "validation_role" %in%
            names(manifest)
        manifest_rows <- bfpwr_sim_manifest_to_search_simulation(manifest)
    }

    if (file.exists(path)) {
        bundle <- readRDS(path)
        simulation <- bundle$simulation
        if (is.null(simulation) || nrow(simulation) == 0) {
            stop("search comparison bundle has no simulation manifest",
                 call. = FALSE)
        }
    } else {
        if (nrow(manifest_rows) == 0) {
            if (!isTRUE(require_search)) {
                cat("skipping package search recomputation; missing ",
                    path, " and no sequential n-search manifest rows\n",
                    sep = "")
                return(invisible(NULL))
            }
            stop("search comparison bundle is missing and the manifest has ",
                 "no sequential n-search rows: ", path, call. = FALSE)
        }
        ensure_dir(validation_dir)
        cat("creating search validation bundle from package manifest; ",
            "missing ", path, "\n", sep = "")
        bundle <- list(
            metadata = list(
                source = "package_verification_manifest",
                created_by = "recompute_package_references.R",
                generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
            ),
            simulation = manifest_rows[FALSE, , drop = FALSE],
            package = empty_package_row(character()),
            comparison = data.frame(),
            summary = data.frame()
        )
        simulation <- bundle$simulation
        manifest_only_bundle <- TRUE
    }

    if (!is.null(manifest)) {
        ## The explicit manifest is now the source of truth for sequential
        ## schedule-crossing searches. Keep fixed and anchor rows from the
        ## legacy bundle, but drop the older opaque sequential crossing subset.
        old_sequential_crossing <- if (nrow(simulation) > 0 &&
                                       all(c("mode", "target_source") %in%
                                           names(simulation))) {
            simulation$mode == "sequential" &
                simulation$target_source == "sequential_schedule_crossing"
        } else {
            rep(FALSE, nrow(simulation))
        }
        simulation <- simulation[!old_sequential_crossing, , drop = FALSE]
        simulation <- bfpwr_sim_bind_rows_fill(
            simulation[!simulation$package_key %in% manifest_rows$package_key,
                       , drop = FALSE],
            manifest_rows)
    }

    keys <- existing_package_keys(bundle, simulation)
    if (nrow(manifest_rows) > 0) {
        manifest_eval <- manifest_rows$manifest_role == "required_n_search"
        if (isTRUE(evaluate_diagnostics)) {
            manifest_eval <- manifest_eval |
                manifest_rows$manifest_role == "boundary_diagnostic"
        }
        keys <- unique(c(keys, manifest_rows$package_key[manifest_eval]))
    }
    selected <- simulation[simulation$package_key %in% keys, , drop = FALSE]
    selected <- limit_search_rows(selected, manifest_rows, max_rows)
    designs <- bfpwr_sim_design_case_set("production")
    bf_priors <- bfpwr_sim_bf_prior_case_set("production")

    rows <- vector("list", nrow(selected))
    for (i in seq_len(nrow(selected))) {
        if (i == 1 || i %% 25 == 0 || i == nrow(selected)) {
            cat("recomputing package search", i, "of", nrow(selected), "\n")
            flush.console()
        }
        rows[[i]] <- eval_package_search_row(selected[i, , drop = FALSE],
                                             designs, bf_priors,
                                             bundle$metadata)
    }
    package <- if (length(rows) == 0) {
        empty_package_row(character())
    } else {
        do.call(rbind, rows)
    }
    comparison <- comparison_with_package(simulation, package)
    metadata <- bundle$metadata
    metadata$generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    metadata$corpus_root <- normalizePath(corpus_root, winslash = "/",
                                          mustWork = FALSE)
    metadata$package_rows <- if (nrow(manifest_rows) > 0) {
        if (isTRUE(manifest_is_package_verification)) {
            if (isTRUE(manifest_only_bundle)) {
                "package_verification_manifest"
            } else {
                "fixed_anchor_plus_package_verification_manifest"
            }
        } else {
            if (isTRUE(manifest_only_bundle)) {
                "sequential_manifest"
            } else {
                "fixed_anchor_plus_sequential_manifest"
            }
        }
    } else {
        "simulation_crossings_only"
    }
    metadata$package_recomputed_by <- "recompute_package_references.R"
    metadata$package_recomputed_rows <- nrow(package)
    metadata$package_recompute_elapsed_seconds <-
        round(proc.time()[["elapsed"]] - t0, 3)
    metadata$sequential_search_manifest_rows <- nrow(manifest_rows)
    metadata$sequential_search_manifest <- if (is.null(manifest_file)) {
        NA_character_
    } else {
        normalizePath(manifest_file, winslash = "/", mustWork = FALSE)
    }
    metadata$sequential_search_diagnostics_evaluated <-
        isTRUE(evaluate_diagnostics)

    provenance <- bfpwr_sim_package_provenance()
    for (name in names(provenance)) {
        metadata[[name]] <- provenance[[name]][[1]]
    }

    out <- list(
        metadata = metadata,
        simulation = simulation,
        package = package,
        comparison = comparison,
        summary = data.frame()
    )
    saveRDS(out, path, compress = "xz")
    invisible(out)
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
    validation_dir <- normalizePath(arg_value(args, "search-validation-dir",
                                              file.path(corpus_root,
                                                        "search-validation")),
                                     winslash = "/", mustWork = FALSE)
    ensure_dir(output_dir)
    root <- repo_root()
    provenance <- bfpwr_sim_package_provenance()
    write_csv(provenance,
              file.path(output_dir, "package-reference-provenance.csv"))
    package_manifest_file <- arg_value(args, "package-manifest", NULL)
    if (is.null(package_manifest_file)) {
        default_package_manifest <-
            bfpwr_sim_package_verification_manifest_file(root)
        if (file.exists(default_package_manifest)) {
            package_manifest_file <- default_package_manifest
        }
    }
    package_manifest <- bfpwr_sim_read_package_verification_manifest(
        package_manifest_file, required = FALSE)

    if (!arg_flag(args, "skip-fixtures")) {
        recompute_fixture_references(
            corpus_root = corpus_root,
            output_dir = output_dir,
            fixture_set = arg_value(args, "fixture-set", NULL),
            case_set = arg_value(args, "case-set", "production"),
            package_manifest = package_manifest
        )
    }
    if (!arg_flag(args, "skip-search")) {
        manifest_file <- if (!is.null(package_manifest_file)) {
            package_manifest_file
        } else {
            arg_value(args, "search-manifest",
                      bfpwr_sim_search_manifest_file(root))
        }
        recompute_search_validation(
            corpus_root = corpus_root,
            validation_dir = validation_dir,
            max_rows = as.numeric(arg_value(args, "max-search-rows", Inf)),
            require_search = arg_flag(args, "require-search"),
            manifest_file = manifest_file,
            evaluate_diagnostics =
                arg_flag(args, "evaluate-search-diagnostics")
        )
    }
}

if (sys.nframe() == 0L) main()
