## Optional diagnostic export in the existing fixture schema.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source_simulation_library()
source_package_checkout()
source(file.path(repo_root(), "simulations/R/corpus-io.R"))
source(file.path(repo_root(), "simulations/R/moment-extension.R"))
source(file.path(repo_root(), "simulations/R/moment-reference-check.R"))
args <- parse_args()
results <- arg_value(args, "results", file.path(repo_root(), "simulations/corpus/moment-v2"))
results <- normalizePath(results, winslash = "/", mustWork = TRUE)
output <- arg_value(args, "output", paste0(results, "-report"))
paths <- c(results, normalizePath(output, winslash = "/", mustWork = FALSE))
if (.Platform$OS.type == "windows") paths <- tolower(paths)
if (paths[1L] == paths[2L]) {
    stop("--output must be separate from the --results directory.", call. = FALSE)
}
bundle <- readRDS(file.path(results, "verification.rds"))
stopifnot(length(bundle$cases) == 70, all(bundle$fixed$nsim == 10000),
          all(bundle$sequential$nsim == 10000))
ensure_dir(file.path(output, "fixture-validation"))
diagnostic <- bundle$numerical_diagnostics
if (is.null(diagnostic)) diagnostic <- bfpwr_sim_moment_reference_check(bundle)
stopifnot(nrow(diagnostic$moment_rows) == nrow(bundle$moments))
check_columns <- intersect(c("case_id", "schedule_id", "evidence_threshold", "EN", "VarN",
    "reference_EN", "reference_VarN", "mcse_EN", "mcse_VarN", "tolerance_EN", "tolerance_VarN"),
    names(bundle$moments))
if (!identical(diagnostic$moment_rows[check_columns], bundle$moments[check_columns])) {
    stop("Numerical diagnostics do not match the selected result bundle")
}
saveRDS(diagnostic, file.path(output, "numerical-diagnostics.rds"))
write.csv(diagnostic$moment_rows, file.path(output, "numerical-diagnostics.csv"), row.names = FALSE)
if (!isTRUE(diagnostic$summary$passed)) stop("Reference convergence diagnostics failed; inspect numerical-diagnostics.csv")
fixtures <- list()
for (alternative in c("less", "greater")) for (mode in c("fixed", "sequential")) {
    id <- paste0("z-nmbf01-", alternative, "-", mode, "-core-v1")
    spec <- if (mode == "fixed") bfpwr_sim_fixed_fixture_spec(id, "z", "moment",
        evidence_thresholds = c(3,10,30), search_targets = bfpwr_sim_fixture_search_targets_core())
    else bfpwr_sim_sequential_fixture_spec(id, "z", "moment", evidence_thresholds = c(3,10,30),
        schedules = bfpwr_sim_z_fixture_schedules_core(),
        search_targets = bfpwr_sim_search_targets(target_prob = numeric(0)))
    spec$bf_prior_ids <- vapply(bundle$cases[vapply(bundle$cases,
        function(x) x$analysis_prior$alternative == alternative, logical(1))], `[[`, "", "bf_prior_id")
    fixtures[[length(fixtures)+1L]] <- spec
    path <- bfpwr_sim_fixture_dir(output, spec)
    ensure_dir(path)
    select <- function(name) {
        x <- bundle[[name]]
        x <- x[x$alternative == alternative, , drop = FALSE]
        x$fixture_set_id <- id
        x
    }
    if (mode == "fixed") {
        tails <- select("fixed")
        parts <- list("spec" = spec, "tail-summary" = tails,
            "decision-summary" = select("decisions"),
            "search-summary" = bfpwr_sim_search_fixed_summary(tails, spec$search_targets))
        diagnostics <- bfpwr_sim_mc_reference_diagnostics(tails, label = id)
    } else {
        cumulative <- select("cumulative")
        parts <- list("spec" = spec, "cumulative-summary" = cumulative,
            "final-summary" = select("final"), "search-summary" =
                bfpwr_sim_search_sequential_summary(cumulative, spec$search_targets))
        diagnostics <- bfpwr_sim_mc_reference_diagnostics(select("sequential"), label = id,
            group_cols = c("alternative", "event"))
        en <- select("moments")
        en$empirical_EN <- en$EN
        write.csv(en, file.path(output, "fixture-validation", paste0(id, "-reference-en-rows.csv")), row.names = FALSE)
    }
    files <- vapply(names(parts), function(name) {
        bfpwr_sim_write_rds(parts[[name]], file.path(path, paste0(name, ".rds")))
    }, character(1))
    overview <- data.frame(fixture_set_id = id, family = "z", mode = mode,
        bf_priors = 35L, search_rows = nrow(parts[["search-summary"]]))
    if (mode == "fixed") {
        overview$tail_rows <- nrow(parts[["tail-summary"]])
        overview$decision_rows <- nrow(parts[["decision-summary"]])
    } else {
        overview$cumulative_rows <- nrow(parts[["cumulative-summary"]])
        overview$final_rows <- nrow(parts[["final-summary"]])
    }
    overview_file <- file.path(path,"overview.csv")
    write.csv(overview, overview_file, row.names = FALSE)
    bfpwr_sim_write_fixture_manifest(path, c(files, overview_file))
    integrity <- bfpwr_sim_validate_fixture_summary(output, id, family="z", mode=mode,
        references=FALSE, output_dir=file.path(output,"fixture-validation"))
    if (!integrity$summary$passed) {
        print(integrity$failures)
        stop("Extension fixture integrity check failed: ",id)
    }
    rows <- diagnostics$diagnostic_rows
    rows$error <- rows$reference_prob - rows$prob
    rows$tail_probability <- pmin(rows$prob, 1-rows$prob)
    scale <- bfpwr_sim_mcse(rows$prob, rows$nsim)
    rows$error_mcse <- ifelse(scale > 0, rows$error/scale, NA_real_)
    if (mode == "sequential") rows$n_looks <- rows$look
    write.csv(rows, file.path(output, "fixture-validation", paste0(id, "-reference-diagnostic-rows.csv")), row.names = FALSE)
    write.csv(diagnostics$summary, file.path(output, "fixture-validation", paste0(id, "-reference-summary.csv")), row.names = FALSE)
}
saveRDS(fixtures, file.path(output, "report-fixtures.rds"))
cat("Prepared four report fixtures with existing schemas.\n")
