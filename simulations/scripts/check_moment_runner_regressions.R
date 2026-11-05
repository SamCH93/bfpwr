## Lightweight engineering regression: real corpus summaries, deterministic package-call probes.
source("simulations/scripts/simulation_support.R")
repo_root <- function() normalizePath(getwd(), winslash = "/")
source_simulation_library()
source_package_checkout()
source("simulations/R/corpus-io.R")
e <- new.env(parent = .GlobalEnv)
source("simulations/R/moment-extension.R", local = e)
case <- e$bfpwr_sim_moment_extension_cases()[[1]]
design <- bfpwr_sim_find_design_case(case$design_case_id)
design$nsim <- 20L
design$chunks <- bfpwr_sim_make_chunks(20L, 10L, design$master_seed)
e$bfpwr_sim_find_design_case <- function(...) design
output <- tempfile("call-probe-", tmpdir = "simulations/logs/runner-fixes")
set.seed(928)
for (j in 1:2) {
    ids <- seq(design$chunks$replicate_start[j], design$chunks$replicate_end[j])
    rows <- expand.grid(n = design$look_grid, replicate_id = ids)
    rows$design_case_id <- design$design_case_id
    rows$chunk_id <- j
    rows$look <- match(rows$n, design$look_grid)
    rows$bf_prior_id <- case$bf_prior_id
    rows$log_bf01 <- rnorm(nrow(rows))
    bfpwr_sim_write_rds(rows, bfpwr_sim_bf_chunk_files(output, case, design)[j])
}
calls <- c(pbf01seq = 0L, nnmbf01 = 0L, nbf01seq = 0L)
e$pbf01seq <- function(n, ...) {
    calls["pbf01seq"] <<- calls["pbf01seq"] + 1L
    h1 <- cumsum(rep(.001, length(n)))
    h0 <- cumsum(rep(.0005, length(n)))
    m <- .bfseq_sample_size_moments(diff(c(0, h1)), diff(c(0, h0)), n)
    c(list(n = n, cumpH1 = h1, cumpH0 = h0, cumpInc = 1 - h1 - h0), m)
}
e$pnmbf01 <- function(n, ...) rep(.01, length(n))
e$nnmbf01 <- function(...) { calls["nnmbf01"] <<- calls["nnmbf01"] + 1L; NA_real_ }
e$nbf01seq <- function(...) { calls["nbf01seq"] <<- calls["nbf01seq"] + 1L; NA_real_ }
result <- e$bfpwr_sim_moment_extension_case(output, case)
print(calls)
stopifnot(nrow(result$fixed) == 300L, nrow(result$sequential) == 2331L,
          nrow(result$moments) == 60L, nrow(result$searches) == 60L,
          calls["nnmbf01"] == 30L, calls["nbf01seq"] == 30L)
cat("Unchanged fixture/search coverage and actual solver-call coverage passed.\n")
stopifnot(calls["pbf01seq"] == 6L)
stopifnot(sum(result$timings$operation == "pbf01seq") == 6L,
          sum(result$timings$operation == "nnmbf01") == 30L,
          sum(result$timings$operation == "nbf01seq") == 30L,
          sum(result$reuse$rows_covered[result$reuse$phase == "search_reference_scan"]) == 300L,
          all(result$moments$max_observed_n <= result$moments$max_n),
          all(result$moments$tolerance_EN == pmax(5*result$moments$mcse_EN,
                                                .002*result$moments$max_observed_n)))
short <- result$moments[result$moments$schedule_id == "start20-by10-looks02", ]
expected <- .bfseq_sample_size_moments(c(.001, .001), c(.0005, .0005), c(20, 30))
stopifnot(all(abs(short$reference_EN - expected$EN) < 1e-12),
          all(abs(short$reference_VarN - expected$VarN) < 1e-12))
cat("Six unique predictions, search reuse, prefix moments and timing schema passed.\n")

## A cancellation-sized negative empirical variance retains each consumer's sign.
n <- c(20, 30)
cumulative <- data.frame(n = n, stop_H1_at_look = 0L, stop_H0_at_look = 0L)
final <- data.frame(nInc = 20L, EN = 30, VarN = -1e-10)
moment <- e$bfpwr_sim_moment_schedule_row(final, cumulative, c(.001, .002), c(.0005, .001), n)
diagnostic <- e$bfpwr_sim_moment_mc_fields(moment, cumulative)
stopifnot(identical(moment$relative_error_VarN, moment$error_VarN/abs(final$VarN)),
          identical(diagnostic$relative_error_VarN, moment$error_VarN/final$VarN),
          identical(moment$z_EN, NA_real_), identical(moment$z_VarN, NA_real_))
final$VarN <- 0
moment <- e$bfpwr_sim_moment_schedule_row(final, cumulative, c(.001, .002), c(.0005, .001), n)
stopifnot(identical(moment$relative_error_VarN, NA_real_))
cat("Negative-variance relative-error signs and zero-denominator behavior passed.\n")
