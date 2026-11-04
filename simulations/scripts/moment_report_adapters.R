## Read-only adapters: retain each cohort's saved values and provenance.
moment_root <- Sys.getenv("BFPWR_SIM_MOMENT_RESULTS", unset = "")
if (!nzchar(moment_root) || !file.exists(file.path(moment_root, "verification.rds"))) {
    stop("Missing moment report input: set BFPWR_SIM_MOMENT_RESULTS to the saved result root")
}
moment_root <- normalizePath(moment_root, winslash = "/", mustWork = TRUE)
moment_bundle <- readRDS(file.path(moment_root, "verification.rds"))
moment_diagnostic <- moment_bundle$numerical_diagnostics
moment_fixtures <- fixtures[vapply(fixtures, function(spec) {
    spec$family == "z" && identical(spec$bf_types, "moment")
}, logical(1))]
moment_case_ids <- vapply(moment_bundle$cases, `[[`, "", "bf_prior_id")
core_case_ids <- vapply(bf_priors, `[[`, "", "bf_prior_id")
design_ids <- vapply(family_design_cases("z"), `[[`, "", "design_case_id")
saved_design_ids <- vapply(moment_bundle$designs, `[[`, "", "design_case_id")
stopifnot(setequal(design_ids, saved_design_ids))
for (design in family_design_cases("z")) {
    saved <- moment_bundle$designs[[match(design$design_case_id, saved_design_ids)]]
    stopifnot(identical(design$generation, saved$generation),
              identical(design$look_grid, saved$look_grid), identical(design$nsim, saved$nsim))
}
stopifnot(length(moment_case_ids) == 70L, !anyDuplicated(moment_case_ids),
          !any(moment_case_ids %in% core_case_ids),
          nrow(moment_bundle$fixed) == 146400L,
          nrow(moment_bundle$sequential) == 222570L,
          nrow(moment_bundle$moments) == 4332L,
          nrow(moment_bundle$searches) == 4200L,
          nrow(moment_bundle$timings) == 5692L,
          identical(moment_bundle$timing_schema_version, 2L),
          isTRUE(moment_diagnostic$summary$passed))
if (!is.null(one_sided)) stopifnot(!any(moment_case_ids %in% one_sided$cases$case_id))
check_columns <- intersect(c("case_id", "schedule_id", "evidence_threshold", "EN", "VarN",
    "reference_EN", "reference_VarN", "mcse_EN", "mcse_VarN", "tolerance_EN", "tolerance_VarN"),
    names(moment_bundle$moments))
stopifnot(identical(moment_diagnostic$moment_rows[check_columns],
                    moment_bundle$moments[check_columns]))
for (name in c("fixed", "sequential", "moments", "searches")) {
    x <- moment_bundle[[name]]
    key <- switch(name,
        fixed = c("case_id", "n", "threshold_id"),
        sequential = c("case_id", "schedule_id", "threshold_pair_id", "look", "event"),
        moments = c("case_id", "schedule_id", "threshold_pair_id"),
        searches = c("case_id", "mode", "threshold", "event", "target"))
    stopifnot(!anyDuplicated(x[key]), setequal(unique(x$case_id), moment_case_ids))
}

moment_reference_rows <- function(mode) {
    rows <- moment_bundle[[if (mode == "fixed") "fixed" else "sequential"]]
    rows$error <- rows$reference_prob - rows$prob
    rows$abs_error <- abs(rows$error)
    rows$tail_probability <- pmin(rows$prob, 1 - rows$prob)
    scale <- bfpwr_sim_mcse(rows$prob, rows$nsim)
    rows$error_mcse <- ifelse(scale > 0, rows$error/scale, NA_real_)
    if (mode == "sequential") {
        rows$empirical <- rows$prob
        rows$reference <- rows$reference_prob
        rows$observed_mcse <- scale
        rows$n_looks <- rows$look
    }
    rows
}

moment_en_rows <- function() {
    rows <- moment_bundle$moments
    rows$empirical_EN <- rows$EN
    rows
}

pool_moment_rows <- function(core, added, new_only = FALSE) {
    if (new_only) return(added)
    stopifnot(nrow(core) > 0L,
              !any(unique(core$bf_prior_id) %in% unique(added$bf_prior_id)))
    out <- bfpwr_sim_bind_rows_fill(core, added)
    stopifnot(nrow(out) == nrow(core) + nrow(added))
    for (name in intersect(names(core), names(added))) {
        stopifnot(isTRUE(all.equal(out[[name]][seq_len(nrow(core))], core[[name]],
                                  tolerance = 0, check.attributes = FALSE)),
                  isTRUE(all.equal(tail(out[[name]], nrow(added)), added[[name]],
                                   tolerance = 0, check.attributes = FALSE)))
    }
    out
}

moment_search_rows <- function() {
    x <- moment_bundle$searches
    probability <- x$actual_power
    fixed <- which(x$mode == "fixed" & is.finite(x$grid_n))
    f <- moment_bundle$fixed
    keys <- with(f, paste(case_id, n, tail, evidence_threshold, sep = "\r"))
    stopifnot(!anyDuplicated(keys))
    hits <- match(with(x[fixed, ], paste(case_id, grid_n, event, threshold, sep = "\r")), keys)
    if (anyNA(hits)) stop("Missing saved moment probability at a fixed search grid solution")
    probability[fixed] <- f$reference_prob[hits]
    probability[x$mode == "fixed" & !is.finite(x$grid_n)] <- NA_real_
    data.frame(family = "z", mode = x$mode, source = "one-sided-moment",
        function_name = ifelse(x$mode == "fixed", "nnmbf01", "nbf01seq"),
        model = "moment", case_id = x$case_id, target = x$target,
        outcome = x$event, threshold = x$threshold,
        simulation_n = x$simulation_n, package_n = x$n_found,
        comparison_n = x$grid_n, package_prob = probability,
        simulation_prob = x$simulation_power, simulation_mcse = x$simulation_mcse,
        stringsAsFactors = FALSE)
}

one_sided_timing_rows <- function(incremental = FALSE) {
    ## Schema 2 already records each Moment+/- probability and search call.
    rows <- data.frame(Function = moment_bundle$timings$operation,
                       seconds = moment_bundle$timings$seconds)
    if (!incremental) {
        normal_searches <- one_sided$searches
        rows <- rbind(
            data.frame(Function = ifelse(normal_searches$mode == "fixed", "nbf01", "nbf01seq"),
                       seconds = normal_searches$seconds),
            rows,
            data.frame(Function = ifelse(one_sided$timings$operation == "fixed", "pbf01", "pbf01seq"),
                       seconds = one_sided$timings$seconds))
    }
    stopifnot(all(is.finite(rows$seconds)))
    rows
}

latex_one_sided_timing_table <- function(incremental = FALSE) {
    timings <- one_sided_timing_rows(incremental)
    summary <- do.call(rbind, lapply(split(timings, timings$Function), function(x) {
        data.frame(Function = x$Function[1], "Timed runs" = nrow(x),
            Median = fmt_seconds(median(x$seconds)),
            P90 = fmt_seconds(as.numeric(quantile(x$seconds, 0.9))),
            Max = fmt_seconds(max(x$seconds)),
            Total = fmt_seconds(sum(x$seconds)), check.names = FALSE)
    }))
    caption <- "One-sided z-test run times."
    if (!incremental) caption <- paste(caption,
        "Each fixed pbf01 run times the H1 and H0 vector calls together;",
        "each other run times one package call.")
    latex_table(summary, caption = caption, keep_together = TRUE, max_width = incremental)
}

report_input_preflight <- function() {
    required <- c(file.path(results_root, "fixture-validation/package-reference-provenance.csv"),
        file.path(results_root, "one-sided-z/verification.rds"),
        file.path(results_root, "search-validation/search-validation-comparison.rds"))
    for (spec in bfpwr_sim_fixture_specs()) {
        if (spec$family == "binomial" && spec$mode == "sequential") next
        summaries <- if (spec$mode == "fixed") c("tail-summary.rds", "decision-summary.rds") else
            c("cumulative-summary.rds", "final-summary.rds")
        required <- c(required,
            file.path(bfpwr_sim_fixture_dir(corpus_root, spec),
                c("spec.rds", "manifest.csv", "search-summary.rds", summaries)),
            file.path(results_root, "fixture-validation", paste0(spec$fixture_set_id,
                "-reference-diagnostic-rows.csv")),
            file.path(results_root, "fixture-validation", paste0(spec$fixture_set_id,
                "-reference-timings.csv")))
        if (spec$mode == "sequential") required <- c(required,
            file.path(results_root, "fixture-validation", paste0(spec$fixture_set_id,
                "-reference-en-rows.csv")))
        if (spec$mode == "sequential" && spec$family == "t") required <- c(required,
            file.path(results_root, "fixture-validation", paste0(spec$fixture_set_id,
                "-reference-summary.csv")))
    }
    missing <- required[!file.exists(required)]
    if (length(missing)) stop("Incomplete baseline report inputs:\n", paste(missing, collapse = "\n"))
    if (is.null(one_sided) || !nrow(reference_provenance)) stop("Missing baseline normal/prior provenance")
    stopifnot(nrow(one_sided$cases) == 74L, !anyDuplicated(one_sided$cases$case_id),
              setequal(unique(one_sided$probabilities$mode), c("fixed", "sequential")),
              setequal(unique(one_sided$searches$mode), c("fixed", "sequential")))
    for (mode in c("fixed", "sequential")) {
        expected_ids <- if (mode == "fixed") {
            one_sided$cases$case_id[one_sided$cases$group == "core"]
        } else one_sided$cases$case_id
        stopifnot(setequal(unique(one_sided$probabilities$case_id[
            one_sided$probabilities$mode == mode]), expected_ids))
    }
    manifest <- bfpwr_sim_read_package_verification_manifest(
        bfpwr_sim_package_verification_manifest_file(repo_root_path), required = TRUE)
    require_cases <- function(rows, expected, label) {
        if (!"package_verification_case_id" %in% names(rows) ||
            !setequal(unique(rows$package_verification_case_id), expected$package_verification_case_id)) {
            stop("Incomplete saved reference-case coverage for ", label)
        }
    }
    for (spec in bfpwr_sim_fixture_specs()) {
        if (spec$family == "binomial" && spec$mode == "sequential") next
        rows <- read.csv(file.path(results_root, "fixture-validation",
            paste0(spec$fixture_set_id, "-reference-diagnostic-rows.csv")))
        if (!nrow(rows) || !any(is.finite(rows$reference_prob))) {
            stop("Missing prior/family reference results for ", spec$fixture_set_id)
        }
        timings <- read.csv(file.path(results_root, "fixture-validation",
            paste0(spec$fixture_set_id, "-reference-timings.csv")))
        if (!any(is.finite(timings$elapsed_seconds))) {
            stop("Missing saved timing records for ", spec$fixture_set_id)
        }
        role <- if (spec$mode == "fixed") "fixed_probability_reference" else
            "sequential_package_reference"
        expected <- manifest[manifest$fixture_set_id == spec$fixture_set_id &
            manifest$validation_role == role, ]
        require_cases(rows, expected, spec$fixture_set_id)
        expected_rows <- if (spec$mode == "fixed") rep(1L, nrow(expected)) else 3L*expected$n_looks
        observed_rows <- table(rows$package_verification_case_id)
        key <- if (spec$mode == "fixed") "package_verification_case_id" else
            c("package_verification_case_id", "look", "tail")
        if (anyDuplicated(rows[key]) || !all(
            as.numeric(observed_rows[expected$package_verification_case_id]) == expected_rows)) {
            stop("Incomplete saved comparison rows for ", spec$fixture_set_id)
        }
        if (spec$mode == "sequential") {
            en <- read.csv(file.path(results_root, "fixture-validation",
                paste0(spec$fixture_set_id, "-reference-en-rows.csv")))
            expected_en <- manifest[manifest$fixture_set_id == spec$fixture_set_id &
                manifest$validation_role %in% c("sequential_package_reference",
                    "sequential_expected_sample_size_reference"), ]
            require_cases(en, expected_en, paste(spec$fixture_set_id, "EN"))
            if (nrow(en) != nrow(expected_en)) stop("Incomplete saved EN rows for ", spec$fixture_set_id)
        }
    }
    search <- readRDS(file.path(results_root, "search-validation/search-validation-comparison.rds"))
    expected_search <- manifest[manifest$validation_role == "sequential_n_search", ]
    required_search <- expected_search$package_key[expected_search$manifest_role == "required_n_search"]
    if (anyDuplicated(search$package$package_key) ||
        anyDuplicated(search$comparison$package_key) ||
        !all(required_search %in% search$package$package_key) ||
        !all(expected_search$package_key %in% search$comparison$package_key)) {
        stop("Incomplete saved baseline search-case coverage")
    }
    moment_search_keys <- expected_search$package_key[expected_search$manifest_role == "required_n_search" &
        expected_search$family == "z" & expected_search$bf_type == "moment"]
    moment_search_elapsed <- search$comparison$package_elapsed_seconds[
        match(moment_search_keys, search$comparison$package_key)]
    if (length(moment_search_elapsed) != length(moment_search_keys) ||
        !all(is.finite(moment_search_elapsed))) {
        stop("Missing saved baseline moment search timing records")
    }
    required
}
baseline_inputs <- report_input_preflight()
