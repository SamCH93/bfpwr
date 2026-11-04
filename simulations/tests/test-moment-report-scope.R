## Run from the repository root with the complete saved baseline and moment inputs.
source("simulations/scripts/simulation_support.R")
args <- parse_args()
baseline <- arg_value(args, "baseline-root")
moment <- arg_value(args, "moment-root")
if (is.null(baseline) || is.null(moment)) {
    stop("Supply --baseline-root and --moment-root with the saved report inputs", call. = FALSE)
}
Sys.setenv(BFPWR_SIM_REPO = getwd(),
    BFPWR_SIM_CORPUS = normalizePath(baseline, winslash = "/", mustWork = TRUE),
    BFPWR_SIM_RESULTS = normalizePath(baseline, winslash = "/", mustWork = TRUE),
    BFPWR_SIM_MOMENT_RESULTS = normalizePath(moment, winslash = "/", mustWork = TRUE))

run_chunk <- function(lines, start) {
    end <- which(seq_along(lines) > start & lines == "@")[1L]
    eval(parse(text = lines[seq.int(start + 1L, end - 1L)]), envir = .GlobalEnv)
}
template <- readLines("simulations/simulation-verification.Rnw", warn = FALSE)
run_chunk(template, grep("^<<fixture-setup,", template))

# Record actual arguments and results while the unmodified incremental chunks run.
captured <- new.env()
capture_result <- function(name) {
    original <- get(name, envir = .GlobalEnv)
    assign(name, function(...) {
        result <- original(...)
        captured[[name]] <- c(captured[[name]], list(result))
        result
    }, envir = .GlobalEnv)
}
for (name in c("design_table", "analysis_prior_table", "fixture_setup_table",
               "pool_moment_rows", "sample_size_comparisons", "one_sided_timing_rows",
               "package_timing_rows", "package_timing_summary_table")) {
    capture_result(name)
}
points <- function(x, y, ...) {
    if (chunk %in% c("moment-fixed-search", "moment-sequential-search")) {
        captured[[chunk]] <- c(captured[[chunk]], list(data.frame(x, y)))
    }
    graphics::points(x, y, ...)
}
body <- readLines("simulations/moment-extension-body.Rnw", warn = FALSE)
grDevices::pdf(NULL, width = 6.8, height = 6.8)
invisible(capture.output(for (start in grep("^<<", body)) {
    chunk <- sub("^<<([^,>]+).*", "\\1", body[start])
    run_chunk(body, start)
}))
invisible(grDevices::dev.off())

designs_shown <- captured$design_table[[1L]]
priors_shown <- captured$analysis_prior_table[[1L]]
settings <- captured$fixture_setup_table
pooled <- captured$pool_moment_rows
search <- captured$sample_size_comparisons[[1L]]
timings <- captured$one_sided_timing_rows[[1L]]
fixed_points <- captured[["moment-fixed-search"]]
sequential_points <- captured[["moment-sequential-search"]]
cat("Design and analysis rows:", nrow(designs_shown), nrow(priors_shown), "\n")
cat("Fixed and sequential setting pairs:",
    unlist(lapply(settings, `[[`, "Prior/design pairs")), "\n")
plot_rows <- c(2L, 4L, 5L)
cat("Fixed, sequential and EN plot rows:", vapply(pooled[plot_rows], nrow, 0L), "\n")
cat("Search rows and plotted pairs:", nrow(search),
    nrow(fixed_points[[1L]]), nrow(sequential_points[[1L]]), "\n")
cat("One-sided timing rows:", nrow(timings), "| baseline timing calls:",
    length(captured$package_timing_rows), "\n")
if (length(captured$package_timing_rows)) {
    legacy <- captured$package_timing_rows[[1L]]
    cat("Baseline moment timing rows:", sum(legacy$family == "z" & legacy$bf_type == "moment"), "\n")
}

cases <- moment_bundle$cases
case_designs <- vapply(cases, `[[`, "", "design_case_id")
alternatives <- vapply(cases, function(case) case$analysis_prior$alternative, "")
case_priors <- do.call(rbind, lapply(cases, function(case) {
    data.frame(Function = format_function_name(case),
        H1 = format_h1_prior(case), H0 = format_h0_prior(case),
        design_case_id = case$design_case_id)
}))
expected_priors <- unique(case_priors[1:3])
expected_priors$Designs <- vapply(seq_len(nrow(expected_priors)), function(i) {
    matches <- with(case_priors, Function == expected_priors$Function[i] &
        H1 == expected_priors$H1[i] & H0 == expected_priors$H0[i])
    format_design_refs(case_priors$design_case_id[matches], "z")
}, "")
expected_priors <- expected_priors[order(expected_priors$Function,
    expected_priors$H1, expected_priors$H0), ]
rownames(expected_priors) <- rownames(priors_shown) <- NULL

fixed <- moment_bundle$fixed
sequential <- moment_bundle$sequential
en <- moment_bundle$moments
raw_search <- moment_bundle$searches
raw_timing <- moment_bundle$timings
raw_plot_rows <- list(fixed, fixed, sequential, sequential, en)
stopifnot(length(moment_case_ids) == 70L,
    identical(as.integer(table(alternatives)), c(35L, 35L)),
    !any(moment_case_ids %in% c(core_case_ids, one_sided$cases$case_id)),
    length(design_ids) == 15L, setequal(case_designs, design_ids),
    identical(designs_shown$Design, unname(design_short_id_map("z")[design_ids])),
    nrow(priors_shown) == 8L, identical(priors_shown, expected_priors),
    identical(settings[[1L]][["Prior/design pairs"]], c(22L, 48L)),
    identical(settings[[2L]][["Prior/design pairs"]], 70L),
    length(unique(sequential$schedule_id)) == 22L,
    length(pooled) == 5L,
    all(vapply(seq_along(pooled), function(i) {
        raw <- raw_plot_rows[[i]]
        # The plotting adapter uses n_looks for the current look index.
        if ("look" %in% names(raw)) raw$n_looks <- raw$look
        identical(pooled[[i]][names(raw)], raw) && setequal(raw$case_id, moment_case_ids)
    }, FALSE)),
    identical(vapply(pooled[plot_rows], nrow, 0L), c(146400L, 222570L, 4332L)))

fixed_grid <- which(raw_search$mode == "fixed" & is.finite(raw_search$grid_n))
fixed_key <- with(fixed, paste(case_id, n, tail, evidence_threshold, sep = "\r"))
search_key <- with(raw_search[fixed_grid, ], paste(case_id, grid_n, event, threshold, sep = "\r"))
expected_prob <- raw_search$actual_power
expected_prob[raw_search$mode == "fixed"] <- NA_real_
expected_prob[fixed_grid] <- fixed$reference_prob[match(search_key, fixed_key)]
stopifnot(nrow(search) == 4200L, identical(as.integer(table(search$mode)), c(2100L, 2100L)),
    all(search$source == "one-sided-moment"), setequal(search$case_id, moment_case_ids),
    identical(search$case_id, raw_search$case_id), identical(search$target, raw_search$target),
    identical(search$mode, raw_search$mode), identical(search$outcome, raw_search$event),
    identical(search$threshold, raw_search$threshold),
    identical(search$simulation_n, raw_search$simulation_n),
    identical(search$package_n, raw_search$n_found),
    identical(search$comparison_n, raw_search$grid_n),
    identical(search$package_prob, expected_prob),
    identical(search$simulation_prob, raw_search$simulation_power),
    identical(search$simulation_mcse, raw_search$simulation_mcse))
for (mode in c("fixed", "sequential")) {
    selected <- search[search$mode == mode & search$comparable, ]
    actual <- if (mode == "fixed") fixed_points else sequential_points
    stopifnot(length(actual) == 2L, nrow(selected) == if (mode == "fixed") 582L else 890L,
        identical(actual[[1L]], data.frame(x = selected$simulation_n, y = selected$package_n)),
        identical(actual[[2L]], data.frame(x = selected$comparison_n, y = selected$error_mcse)))
}
timing_counts <- table(timings$Function)
stopifnot(length(captured$package_timing_rows) == 0L,
    length(captured$package_timing_summary_table) == 0L, nrow(timings) == 5692L,
    setequal(raw_timing$case_id, moment_case_ids),
    identical(timings$Function, raw_timing$operation), identical(timings$seconds, raw_timing$seconds),
    all(is.finite(timings$seconds)),
    identical(names(timing_counts), c("nbf01seq", "nnmbf01", "pbf01seq", "pnmbf01")),
    identical(as.integer(timing_counts), c(2100L, 2100L, 420L, 1072L)))
cat("Incremental tables, plot inputs, searches, and timings use only saved Moment+/- results.\n")
