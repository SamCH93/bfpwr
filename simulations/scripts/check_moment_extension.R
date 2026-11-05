script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source_simulation_library()
source_package_checkout()
source(file.path(repo_root(), "simulations/R/bayes-factors.R"))
source(file.path(repo_root(), "simulations/R/moment-extension.R"))
source(file.path(repo_root(), "simulations/R/moment-identity.R"))
cat(R.version.string, "\n", R.home(), "\n")
print(.libPaths())
stopifnot(requireNamespace("mvtnorm"), requireNamespace("lamW"), requireNamespace("knitr"))
cases <- bfpwr_sim_moment_extension_cases()
plan <- bfpwr_sim_moment_extension_plan(cases)
stopifnot(length(cases) == 70, nrow(plan) == 700,
          sum(vapply(cases, function(x) x$analysis_prior$alternative == "less", logical(1))) == 35)
designs <- lapply(cases, function(x) bfpwr_sim_find_design_case(x$design_case_id))
schedules <- bfpwr_sim_z_fixture_schedules_core()
compatible <- lapply(designs, function(d) schedules[vapply(schedules,
    bfpwr_sim_schedule_is_compatible, logical(1), design = d)])
coverage <- c(fixed = sum(vapply(designs, function(d) 6*length(d$look_grid), numeric(1))),
    sequential = sum(vapply(compatible, function(s) 9*sum(vapply(s,
        function(x) length(x$n), integer(1))), numeric(1))),
    moments = 3*sum(lengths(compatible)), searches = length(cases)*60)
stopifnot(identical(unname(coverage), c(146400, 222570, 4332, 4200)))
print(coverage)
for (i in seq(1, 70, by = 2)) {
    a <- cases[[i]]; b <- cases[[i+1]]
    prior <- a$analysis_prior
    estimate <- prior$null + c(-.3, 0, .3)
    data <- data.frame(estimate = estimate, se = .2)
    actual <- bfpwr_sim_log_bf(a, data)
    expected <- nmbf01(estimate, se = .2, null = prior$null, psd = prior$psd,
                       alternative = "less", log = TRUE)
    stopifnot(identical(actual, expected))
    data$estimate <- 2*prior$null - estimate
    stopifnot(isTRUE(all.equal(actual, bfpwr_sim_log_bf(b, data), tolerance = 1e-12)))
}
## Direct first-crossing checks distinguish per-look from stopped trajectories.
x <- matrix(c(log(.05), log(20), log(2), log(20), log(.05), log(.05)), nrow = 2, byrow = TRUE)
counts <- bfpwr_sim_sequential_counts(x, data.frame(k1 = .1, k0 = 10))
stopifnot(identical(counts$stop_h1, c(1,0,0)), identical(counts$stop_h0, c(1,0,0)))
source_files <- bfpwr_sim_moment_source_files(repo_root(), "reference")
for (f in file.path(repo_root(), source_files[grepl("[.]R$", source_files)])) parse(f)
cat("70 cases, 700 chunks; BF forwarding, reflection, stopping, and R syntax checks passed.\n")
source(file.path(repo_root(), "simulations/R/corpus-io.R"))
original_find <- bfpwr_sim_find_design_case
design <- original_find(cases[[1]]$design_case_id)
design$look_grid <- seq(10,50,10)
design$nsim <- 20L
design$design_prior$mean <- -.6
design$chunks <- bfpwr_sim_make_chunks(20L, 10L, design$master_seed)
bfpwr_sim_find_design_case <- function(...) design
input <- file.path(repo_root(), "simulations/logs/smoke-input")
output <- file.path(repo_root(), "simulations/logs/smoke-output")
bfpwr_sim_write_rds(design, file.path(bfpwr_sim_design_dir(input, design), "design.rds"))
set.seed(728)
for (j in 1:2) {
    rows <- do.call(rbind, lapply(seq(design$chunks$replicate_start[j], design$chunks$replicate_end[j]), function(id) {
        data.frame(design_case_id = design$design_case_id, replicate_id = id,
                   chunk_id = j, look = 1:5, n = design$look_grid,
                   estimate = cumsum(rnorm(5, mean = -6, sd = sqrt(10)))/design$look_grid,
                   se = 1/sqrt(design$look_grid))
    }))
    bfpwr_sim_write_rds(rows, bfpwr_sim_design_chunk_files(input, design)[j])
    index <- data.frame(chunk_id = 1:2, replicate_start = design$chunks$replicate_start,
        replicate_end = design$chunks$replicate_end, n_rows = 50L, bytes = file.info(bfpwr_sim_design_chunk_files(input, design))$size,
        sha256 = unname(tools::sha256sum(bfpwr_sim_design_chunk_files(input, design))),
        file = file.path("chunks", basename(bfpwr_sim_design_chunk_files(input, design))))
    bfpwr_sim_write_rds(index, file.path(bfpwr_sim_design_dir(input, design), "chunk-index.rds"))
}
for (j in 1:2) bfpwr_sim_moment_extension_chunk(input, output, cases[[1]], j)
result <- bfpwr_sim_moment_extension_case(output, cases[[1]], max_looks = 2)
stopifnot(all(result$fixed$nsim == 20), all(result$sequential$nsim == 20),
          nrow(result$searches) == 30, all(is.finite(result$fixed$reference_prob)),
          all(is.finite(result$sequential$reference_prob)),
          any(is.finite(result$searches$simulation_n[result$searches$event == "H1"])))
cat("Synthetic two-chunk fixture and package-reference smoke check passed.\n")
