## Recompute package predictions against the existing simulation corpus.
## Separate output directories retain the original simulation and prediction data.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source(file.path(dirname(script), "recompute_package_references.R"))
source(file.path(dirname(script), "..", "R", "one-sided-z.R"))

args <- parse_args()
root <- repo_root()
corpus_root <- normalizePath(arg_value(args, "corpus-root",
    file.path(root, "simulations", "corpus", "v1")), winslash = "/", mustWork = TRUE)
if (!is.null(args$ngrid)) {
    stop("verification uses package defaults; --ngrid is no longer supported")
}
ngrid <- as.numeric(.bfseq_integration_settings(list())$ngrid)
numerical_defaults <- bfpwrOptions()
workers <- as.numeric(arg_value(args, "workers", 6))
stopifnot(length(ngrid) == 1L, is.finite(ngrid), ngrid >= 1, ngrid == floor(ngrid),
          length(workers) == 1L, is.finite(workers), workers >= 1,
          workers == floor(workers))
results_root <- normalizePath(arg_value(args, "results-root",
    file.path(corpus_root, paste0("predictions-default-ngrid-", ngrid))),
    winslash = "/", mustWork = FALSE)
output <- file.path(results_root, "one-sided-z")
ensure_dir(file.path(output, "case-results"))
assemble_only <- arg_flag(args, "assemble-only")

for (directory in c("fixture-validation", "search-validation", "reference-results")) {
    ensure_dir(file.path(results_root, directory))
}
original <- readRDS(file.path(corpus_root, "one-sided-z", "verification.rds"))
cases <- original$cases
manifest_file <- bfpwr_sim_package_verification_manifest_file(root)
manifest <- bfpwr_sim_read_package_verification_manifest(manifest_file, required = TRUE)
specs <- fixture_specs(corpus_root)
reference_cases <- manifest[manifest$validation_role == "sequential_package_reference", ]
provenance_file <- file.path(results_root, "fixture-validation",
                             "package-reference-provenance.csv")
if (assemble_only) {
    provenance <- utils::read.csv(provenance_file, stringsAsFactors = FALSE)
    bfpwr_sim_check_cache(bfpwr_sim_cache_identity(provenance),
        bfpwr_sim_cache_identity(), "parent provenance")
} else {
    provenance <- bfpwr_sim_package_provenance()
    provenance$integration_rng_kind <- "Mersenne-Twister/Inversion/Rejection"
    provenance$prediction_workers <- workers
    write_csv(provenance, provenance_file)
}

identity <- bfpwr_sim_cache_identity(provenance)

if (!assemble_only) {
    cluster <- parallel::makePSOCKcluster(workers,
        rscript = file.path(R.home("bin"),
            if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
        outfile = file.path(results_root, "workers.log"))
    parallel::clusterExport(cluster, c("root", "corpus_root", "results_root",
                                        "output", "ngrid", "numerical_defaults",
                                        "manifest", "manifest_file",
                                        "reference_cases"))
    invisible(parallel::clusterEvalQ(cluster, {
        source(file.path(root, "simulations", "scripts", "recompute_package_references.R"))
        source(file.path(root, "simulations", "R", "one-sided-z.R"))
        bfpwrOptions(numerical_defaults)
        stopifnot(identical(bfpwrOptions(), numerical_defaults))
        ## The package's fixed integration seed also depends on RNGkind().
        RNGkind("Mersenne-Twister", "Inversion", "Rejection")
        designs <- bfpwr_sim_design_case_set("production")
        bf_priors <- bfpwr_sim_bf_prior_case_set("production")
        fixture_cache <- new.env(parent = emptyenv())
        read_fixture <- function(id, family, mode) {
            if (!exists(id, fixture_cache, inherits = FALSE)) {
                assign(id, bfpwr_sim_read_fixture_summary(
                    corpus_root, id, family, mode), fixture_cache)
            }
            get(id, fixture_cache, inherits = FALSE)
        }
        NULL
    }))

    ## Dispatch sequential references by case so expensive strict integrations
    ## do not leave one worker processing an entire fixture on its own.
    fixed <- specs[specs$mode == "fixed", ]
    jobs <- c(lapply(seq_len(nrow(fixed)), function(i) {
        list(kind = "fixed", id = fixed$fixture_set_id[i], family = fixed$family[i])
    }), lapply(seq_len(nrow(reference_cases)), function(i) {
        list(kind = "sequential", id = reference_cases$package_verification_case_id[i])
    }), list(list(kind = "search", id = "legacy-searches")),
    lapply(cases$case_id, function(id) list(kind = "one-sided", id = id)))
    cat("Recomputing", length(jobs), "tasks with ngrid =", ngrid, "\n")
    completed <- tryCatch(parallel::parLapplyLB(cluster, jobs, function(job) {
        cat(format(Sys.time()), "Starting", job$id, "\n")
        if (job$kind == "fixed") {
            fixture <- read_fixture(job$id, job$family, "fixed")
            selected <- manifest[manifest$fixture_set_id == job$id, ]
            reference <- bfpwr_sim_validate_fixed_references(fixture,
                designs, bf_priors, manifest_cases = selected)
            write_reference_result(reference,
                file.path(results_root, "fixture-validation", job$id))
            saveRDS(bfpwr_sim_cache_identity(), file.path(results_root,
                "fixture-validation", job$id, "cache-identity.rds"))
        } else if (job$kind == "sequential") {
            selected <- reference_cases[
                reference_cases$package_verification_case_id == job$id, ]
            fixture <- read_fixture(selected$fixture_set_id, selected$family, "sequential")
            evaluate <- if (selected$family == "t") {
                bfpwr_sim_t_sequential_reference_table
            } else bfpwr_sim_z_sequential_reference_table
            reference <- evaluate(fixture, designs, bf_priors, strict = TRUE,
                                  cases = selected)
            stopifnot(nrow(reference$reference_rows) == 3L*selected$n_looks,
                      nrow(reference$en_rows) == 1L,
                      all(is.finite(reference$reference_rows$reference_prob)))
            reference$cache_identity <- bfpwr_sim_cache_identity()
            saveRDS(reference, file.path(results_root, "reference-results",
                                        paste0(job$id, ".rds")))
        } else if (job$kind == "search") {
            recompute_search_validation(corpus_root,
                file.path(results_root, "search-validation"), require_search = TRUE,
                manifest_file = manifest_file)
            path <- file.path(results_root, "search-validation",
                               "search-validation-comparison.rds")
            bundle <- readRDS(path)
            bundle$metadata$cache_identity <- bfpwr_sim_cache_identity()
            bundle$metadata$integration_grid <- ngrid
            bundle$metadata$integration_rng_kind <- paste(RNGkind(), collapse = "/")
            saveRDS(bundle, path, compress = "xz")
        } else {
            cached <- readRDS(file.path(corpus_root, "one-sided-z", "bayes-factors",
                                        paste0(job$id, ".rds")))
            evaluated <- bfpwr_sim_one_sided_z_evaluate(corpus_root,
                cached$case, cached = cached, search = TRUE)
            result <- evaluated$result
            saveRDS(result, file.path(output, "case-results", paste0(job$id, ".rds")))
        }
        cat(format(Sys.time()), "Finished", job$id, "\n")
        job$id
    }, chunk.size = 1L), finally = parallel::stopCluster(cluster))
    stopifnot(identical(unlist(completed), vapply(jobs, `[[`, character(1), "id")))
}

## Diagnose the combined fixture rows once, preserving the same Holm correction
## and Monte Carlo tolerances as a serial reference recomputation.
for (id in unique(reference_cases$fixture_set_id)) {
    selected <- reference_cases[reference_cases$fixture_set_id == id, ]
    references <- lapply(selected$package_verification_case_id, function(key) {
        reference <- readRDS(file.path(results_root, "reference-results", paste0(key, ".rds")))
        bfpwr_sim_check_cache(reference$cache_identity, identity, key)
        reference
    })
    combined <- lapply(c("reference_rows", "en_rows", "timings"), function(name) {
        bfpwr_sim_bind_rows_fill(lapply(references, `[[`, name))
    })
    names(combined) <- c("reference_rows", "en_rows", "timings")
    stopifnot(nrow(combined$reference_rows) == 3L*sum(selected$n_looks),
              nrow(combined$en_rows) == nrow(selected))
    fixture <- bfpwr_sim_read_fixture_summary(corpus_root, id, selected$family[1],
                                             "sequential")
    reference <- bfpwr_sim_validate_sequential_references(fixture,
        designs = NULL, bf_priors = NULL, reference = combined)
    write_reference_result(reference, file.path(results_root, "fixture-validation", id))
}
status <- specs[c("fixture_set_id", "family", "mode")]
unsupported <- status$family == "binomial" & status$mode == "sequential"
status$action <- ifelse(unsupported, "skipped", "recomputed")
status$reason <- ifelse(unsupported, "sequential binomial has no package probability API", "")
write_csv(status, file.path(results_root, "fixture-validation",
                            "package-reference-recompute-status.csv"))

for (id in specs$fixture_set_id[specs$mode == "fixed"]) {
    bfpwr_sim_check_cache(readRDS(file.path(results_root, "fixture-validation",
        id, "cache-identity.rds")), identity, id)
}
results <- bfpwr_sim_one_sided_z_results(output, cases, identity)
search_metadata <- readRDS(file.path(results_root, "search-validation",
                                    "search-validation-comparison.rds"))$metadata
status <- utils::read.csv(file.path(results_root, "fixture-validation",
                                    "package-reference-recompute-status.csv"))
bfpwr_sim_check_cache(search_metadata$cache_identity, identity, "sample-size searches")
stopifnot(identical(as.numeric(search_metadata$integration_grid), ngrid),
          setequal(status$fixture_set_id, fixture_specs(corpus_root)$fixture_set_id))
combine <- function(name) do.call(rbind, lapply(results, `[[`, name))
probabilities <- combine("probabilities")
## Exact equality of the empirical columns verifies that no trajectories or
## stopping rules changed during the refresh. Saved diagnostics may be sorted
## differently, so align by condition before comparing.
probability_key <- function(x) {
    do.call(paste, c(x[c("case_id", "mode", "schedule", "pair", "n", "outcome")],
                    sep = "\r"))
}
old_key <- probability_key(original$probabilities)
new_key <- probability_key(probabilities)
stopifnot(!anyDuplicated(old_key), !anyDuplicated(new_key),
          setequal(old_key, new_key))
probabilities <- probabilities[match(old_key, new_key), ]
rownames(probabilities) <- NULL
original_probabilities <- original$probabilities
rownames(original_probabilities) <- NULL
empirical <- setdiff(names(probabilities), c("reference_prob", "error"))
stopifnot(identical(probabilities[empirical], original_probabilities[empirical]))
provenance$replicates_per_condition <- original$provenance$replicates_per_condition
provenance$source_archive <- original$provenance$source_archive
provenance$source_archive_sha256 <- original$provenance$source_archive_sha256
bundle <- c(list(cases = cases, provenance = provenance,
                 probabilities = probabilities, searches = combine("searches"),
                 timings = combine("timings"), bf_checks = combine("bf_checks")),
            bfpwr_sim_one_sided_z_diagnostics(probabilities, combine("moments")))
saveRDS(bundle, file.path(output, "verification.rds"))
for (name in c("cases", "provenance", "probabilities", "searches", "moments",
               "timings", "bf_checks", "probability_summary")) {
    write_csv(bundle[[name]], file.path(output, paste0(name, ".csv")))
}
failures <- bundle$mc_diagnostics$diagnostic_rows
write_csv(failures[failures$holm_p < 0.01, ], file.path(output, "probability_failures.csv"))
print(bundle$probability_summary)
print(bundle$mc_diagnostics$summary)
cat("Moment failures:", sum(!bundle$moments$EN_ok), "means;",
    sum(!bundle$moments$VarN_ok), "variances\n")

## Finish both sets of diagnostics and render their results even when a check
## fails. Preserve the failing exit status after producing the report.
validation_errors <- bfpwr_sim_one_sided_z_failures(bundle)
validate <- function(script, arguments) {
    tryCatch(run_rscript(file.path(root, "simulations", "scripts", script), arguments),
        error = function(e) {
            validation_errors <<- c(validation_errors, conditionMessage(e))
        })
}
validate("validate_fixture_suite.R",
         c("--corpus-root", corpus_root, "--output-dir",
           file.path(results_root, "fixture-validation")))
validate("refresh_search_validation.R",
         c("--corpus-root", corpus_root, "--validation-dir",
           file.path(results_root, "search-validation")))
if (!arg_flag(args, "skip-report")) {
    run_rscript(file.path(root, "simulations", "scripts",
                          "render_simulation_verification_background.R"),
                c("--corpus-root", corpus_root, "--results-root", results_root))
}
if (length(validation_errors) > 0L) {
    stop(paste(validation_errors, collapse = "\n"), call. = FALSE)
}
