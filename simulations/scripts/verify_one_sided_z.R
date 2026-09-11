## Reuse the released data; no random observations or extra replicates are drawn.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source_simulation_library()
source(file.path(repo_root(), "simulations", "R", "one-sided-z.R"))

args <- parse_args()
root <- repo_root()
corpus_root <- normalizePath(arg_value(args, "corpus-root",
    file.path(root, "simulations", "corpus", "v1")), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(arg_value(args, "package-root", root),
                              winslash = "/", mustWork = TRUE)
source_package_checkout(package_root)
output <- file.path(corpus_root, "one-sided-z")
ensure_dir(file.path(output, "bayes-factors"))
ensure_dir(file.path(output, "case-results"))
cases <- bfpwr_sim_one_sided_z_cases()
provenance <- bfpwr_sim_package_provenance(package_root)
provenance$package_git_executable <- unname(Sys.which("git"))
provenance$replicates_per_condition <- 10000L
provenance$source_archive <- "bfpwr-sim-corpus-v1-designs.tar.zst"
checksums <- readLines(file.path(root, "simulations", "downloads", "sim-corpus-v1",
                                 "SHA256SUMS.txt"))
provenance$source_archive_sha256 <- strsplit(
    grep(provenance$source_archive, checksums, fixed = TRUE, value = TRUE),
    "[[:space:]]+")[[1]][1]
utils::write.csv(cases, file.path(output, "cases.csv"), row.names = FALSE)

bf_checks <- list()
for (design_id in unique(cases$design_id)) {
    subset <- cases[cases$design_id == design_id, ]
    cat("Reading", design_id, "\n")
    data <- bfpwr_sim_one_sided_z_read(corpus_root, subset[1, ])
    for (j in seq_len(nrow(subset))) {
        case <- subset[j, ]
        cat("  Calculating", case$prior_id, case$alternative, "\n")
        se <- rep(1/sqrt(data$n), each = 10000)
        logbf <- bfpwr_sim_one_sided_z_logbf(
            as.vector(data$estimate), se, case$pm, case$psd, case$alternative)
        stopifnot(all(is.finite(logbf)))
        selected <- unique(round(seq(1, length(logbf), length.out = 250)))
        package_bf <- bf01(as.vector(data$estimate)[selected], se[selected],
            pm = case$pm, psd = case$psd, alternative = case$alternative, log = TRUE)
        error <- max(abs(logbf[selected] - package_bf))
        stopifnot(error < 1e-8)
        bf_checks[[length(bf_checks) + 1L]] <- data.frame(
            case_id = case$case_id, checked = length(selected), max_logbf_error = error)
        logbf <- matrix(logbf, nrow = 10000)
        saveRDS(list(case = case, n = data$n, replicate_id = data$replicate_id,
                     true_effect = data$true_effect, log_bf01 = logbf),
                file.path(output, "bayes-factors", paste0(case$case_id, ".rds")))
        result <- bfpwr_sim_one_sided_z_compare(case, logbf, data$n)
        saveRDS(result, file.path(output, "case-results", paste0(case$case_id, ".rds")))
    }
}
results <- lapply(cases$case_id, function(id) {
    readRDS(file.path(output, "case-results", paste0(id, ".rds")))
})
combine <- function(name) do.call(rbind, lapply(results, `[[`, name))
probabilities <- combine("probabilities")
moments <- combine("moments")
cat("Checking sample-size searches\n")
searches <- bfpwr_sim_one_sided_z_search(cases, probabilities)
diagnostics <- bfpwr_sim_one_sided_z_diagnostics(probabilities, moments)
cat("Checking integration convergence on selected cases\n")
convergence <- bfpwr_sim_one_sided_z_convergence(cases, probabilities)
bundle <- c(list(cases = cases, provenance = provenance,
                 probabilities = probabilities, searches = searches,
                 timings = combine("timings"), bf_checks = do.call(rbind, bf_checks),
                 convergence = convergence),
            diagnostics)
saveRDS(bundle, file.path(output, "verification.rds"))
for (name in c("provenance", "probabilities", "searches", "moments",
               "timings", "bf_checks", "probability_summary", "convergence")) {
    utils::write.csv(bundle[[name]], file.path(output, paste0(name, ".csv")),
                     row.names = FALSE)
}
diagnostic_rows <- bundle$mc_diagnostics$diagnostic_rows
utils::write.csv(diagnostic_rows[diagnostic_rows$holm_p < 0.01, ],
                 file.path(output, "probability_failures.csv"), row.names = FALSE)
print(bundle$probability_summary)
print(bundle$mc_diagnostics$summary)
cat("Searches:", nrow(searches), "; invalid:", sum(!searches$search_valid), "\n")
cat("Moment checks outside tolerance:",
    sum(!bundle$moments$EN_ok), "EN;", sum(!bundle$moments$VarN_ok), "variance\n")

figure_dir <- file.path(root, "simulations", "figure")
ensure_dir(figure_dir)
for (kind in c("curves", "agreement")) {
    grDevices::png(file.path(figure_dir, paste0("one-sided-z-", kind, ".png")),
                   width = 1800, height = if (kind == "curves") 1800 else 1350,
                   res = 180)
    bfpwr_sim_one_sided_z_plot(bundle, kind)
    grDevices::dev.off()
}
if (!arg_flag(args, "skip-report")) {
    run_rscript(file.path(root, "simulations", "scripts",
                          "render_simulation_verification_background.R"),
                c("--corpus-root", corpus_root))
}
