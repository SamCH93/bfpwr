## Run only the additional one-sided moment-prior cases on an existing corpus.
## Sourcing this file only defines the entry point.
bfpwr_sim_run_moment_extension <- function(args, root) {
    source(file.path(root, "simulations/scripts/simulation_support.R"), local = .GlobalEnv)
    source_simulation_library(root)
    source_package_checkout(root)
    for (name in c("corpus-io", "corpus-validation", "bayes-factors", "moment-extension", "moment-identity")) {
        source(file.path(root, "simulations", "R", paste0(name, ".R")), local = .GlobalEnv)
    }
    args <- parse_args(args)
    phase <- arg_value(args, "phase", "plan")
    output <- arg_value(args, "output", file.path(root, "simulations/corpus/moment-v2"))
    input <- arg_value(args, "input-corpus")
    provenance_path <- file.path(output, "provenance.rds")
    if (is.null(input) && file.exists(provenance_path)) input <- readRDS(provenance_path)$source_corpus
    if (is.null(input)) stop("--input-corpus is required for a fresh plan")
    worker <- as.integer(arg_value(args, "worker", Sys.getenv("PBS_ARRAY_INDEX", "1")))
    workers <- as.integer(arg_value(args, "workers", "1"))
    stopifnot(worker >= 1L, workers >= worker)
    ensure_dir(output)
    cases <- bfpwr_sim_moment_extension_cases()
    plan <- bfpwr_sim_moment_extension_plan(cases)
    design_ids <- unique(vapply(cases, `[[`, "", "design_case_id"))
    designs <- setNames(lapply(design_ids, bfpwr_sim_find_design_case), design_ids)
    inputs <- bfpwr_sim_moment_input_indexes(input, designs)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    bf_identity <- bfpwr_sim_moment_identity(root, "bf", cases, designs, inputs)
    identity <- bfpwr_sim_moment_identity(root, "reference", cases, designs, inputs)
    identity_path <- file.path(output, "calculation-identity.rds")
    if (file.exists(identity_path)) {
        bfpwr_sim_moment_check_identity(readRDS(file.path(output, "bf-dependency-identity.rds")),
                                       bf_identity, "BF dependency identity")
        bfpwr_sim_moment_verify_sources(root, output, "bf")
        if (!phase %in% c("chunks", "finalize")) {
            bfpwr_sim_moment_check_identity(readRDS(identity_path), identity, "reference identity")
            bfpwr_sim_moment_verify_sources(root, output, "reference")
        }
    } else if (phase == "plan") {
        saveRDS(identity, identity_path)
        saveRDS(identity, file.path(output, "reference-identity.rds"))
        saveRDS(bf_identity, file.path(output, "bf-dependency-identity.rds"))
        saveRDS(bf_identity, file.path(output, "bf-evaluation-identity.rds"))
        bf_manifest <- bfpwr_sim_moment_snapshot(root, output, "bf")
        manifest <- bfpwr_sim_moment_snapshot(root, output, "reference")
        saveRDS(bfpwr_sim_moment_provenance(root, input, manifest, cases, designs, inputs), provenance_path)
        saveRDS(bfpwr_sim_moment_provenance(root, input, bf_manifest, cases, designs, inputs, phase = "bf"),
                file.path(output, "bf-evaluation-provenance.rds"))
    } else stop("Run --phase plan before workers")
    if (phase == "plan") {
        saveRDS(cases, file.path(output, "cases.rds"))
        saveRDS(designs, file.path(output, "designs.rds"))
        saveRDS(inputs, file.path(output, "input-indexes.rds"))
        write.csv(plan, file.path(output, "chunk-plan.csv"), row.names = FALSE)
        for (case in cases) bfpwr_sim_write_rds(case,
            file.path(bfpwr_sim_bf_prior_dir(output, case), "bf-prior.rds"))
        cat("Cases:", length(cases), "chunks:", nrow(plan), "\n")
    } else if (phase == "import") {
        source <- arg_value(args, "reuse-from")
        bf_source <- arg_value(args, "bf-source-root")
        if (is.null(source) || is.null(bf_source)) stop("import requires --reuse-from and --bf-source-root")
        if (normalizePath(source, winslash = "/") == normalizePath(output, winslash = "/")) {
            stop("import destination must be a NEW moment-v2 directory")
        }
        bfpwr_sim_moment_import_v1(source, output, input, root, bf_source,
                                   official_metadata = arg_value(args, "official-metadata"))
        cat("Verified historical BF chunks imported; all reference cases remain to be recomputed.\n")
    } else if (phase %in% c("chunks", "summaries")) {
        if (phase == "chunks" && file.exists(file.path(output, "bf-import.rds"))) {
            stop("Imported BF chunks retain their historical provenance; do not evaluate them again here")
        }
        count <- if (phase == "chunks") nrow(plan) else length(cases)
        directory <- file.path(output, if (phase == "chunks") "chunk-status" else "case-results")
        ensure_dir(directory)
        started <- proc.time()[["elapsed"]]
        assigned <- seq(worker, count, by = workers)
        order <- c(assigned, setdiff(seq_len(count), assigned))
        for (i in order) {
            target <- file.path(directory, sprintf("%04d.rds", i))
            lock <- paste0(target, ".lock")
            expected <- if (phase == "chunks") bf_identity else identity
            if (file.exists(target)) {
                bfpwr_sim_moment_check_identity(readRDS(target)$identity, expected, target)
                if (phase == "chunks") {
                    row <- plan[i, ]
                    design <- designs[[cases[[row$case_index]]$design_case_id]]
                    bfpwr_sim_moment_verify_input(input, design, row$chunk_id)
                    path <- bfpwr_sim_bf_chunk_files(output, cases[[row$case_index]], design)[row$chunk_id]
                    stopifnot(identical(unname(tools::sha256sum(path)), readRDS(target)$output_sha256))
                }
                next
            }
            if (!dir.create(lock, showWarnings = FALSE)) next
            cat(format(Sys.time()), phase, i, "of", count, "\n")
            task_started <- Sys.time()
            task_start <- proc.time()[["elapsed"]]
            result <- tryCatch({
                sources <- bfpwr_sim_moment_verify_sources(root, output,
                    if (phase == "chunks") "bf" else "reference")
                if (phase == "chunks") {
                    row <- plan[i, ]
                    bfpwr_sim_moment_extension_chunk(input, output, cases[[row$case_index]], row$chunk_id)
                } else {
                    case <- cases[[i]]
                    design <- designs[[case$design_case_id]]
                    status_ids <- which(plan$case_index == i)
                    for (j in status_ids) {
                        bfpwr_sim_moment_verify_input(input, design, plan$chunk_id[j])
                        status <- readRDS(file.path(output, "chunk-status", sprintf("%04d.rds", j)))
                        path <- bfpwr_sim_bf_chunk_files(output, case, design)[plan$chunk_id[j]]
                        stopifnot(identical(unname(tools::sha256sum(path)), status$output_sha256))
                    }
                    bfpwr_sim_moment_extension_case(output, case)
                }
            }, error = function(e) {
                writeLines(conditionMessage(e), paste0(target, ".error"))
                stop(e)
            })
            if (phase == "summaries") {
                total <- proc.time()[["elapsed"]] - task_start
                result$case_timings$preflight_seconds <- total - result$case_timings$elapsed_seconds
                result$case_timings$started <- task_started
                result$case_timings$completed <- Sys.time()
                result$case_timings$elapsed_seconds <- total
            }
            result$identity <- expected
            result$actual_source_sha256 <- setNames(sources$sha256, sources$file)
            result$completed <- Sys.time()
            saveRDS(result, paste0(target, ".tmp"), compress = "xz")
            stopifnot(file.rename(paste0(target, ".tmp"), target))
            unlink(lock, recursive = TRUE)
            ## Leave time for a graceful exit before the six-hour PBS limit.
            if (proc.time()[["elapsed"]] - started > 5*3600) break
        }
    } else if (phase == "finalize") {
        ensure_dir(file.path(output, "index-status"))
        for (i in seq(worker, length(cases), by = workers)) {
            case <- cases[[i]]
            design <- designs[[case$design_case_id]]
            target <- file.path(output, "index-status", sprintf("%04d.rds", i))
            index <- bfpwr_sim_write_bf_prior_chunk_index(output, case, design, input_root = input)
            stopifnot(nrow(index) == nrow(design$chunks),
                sum(index$n_rows) == design$nsim*length(design$look_grid),
                all(index$n_na_log_bf01 == 0), all(index$n_nan_log_bf01 == 0),
                all(!is.na(index$source_design_chunk_sha256)))
            saveRDS(index, target)
            cat("Indexed and validated case", i, "\n")
        }
    } else if (phase == "assemble") {
        chunk_files <- file.path(output, "chunk-status", sprintf("%04d.rds", seq_len(nrow(plan))))
        stopifnot(all(file.exists(chunk_files)))
        chunks <- lapply(chunk_files, readRDS)
        chunk_identity <- readRDS(file.path(output, "bf-evaluation-identity.rds"))
        stopifnot(all(vapply(chunks, function(x) identical(x$identity, chunk_identity), logical(1))))
        for (i in seq_len(nrow(plan))) {
            case <- cases[[plan$case_index[i]]]
            design <- designs[[case$design_case_id]]
            path <- bfpwr_sim_bf_chunk_files(output, case, design)[plan$chunk_id[i]]
            stopifnot(identical(unname(tools::sha256sum(path)), chunks[[i]]$output_sha256))
        }
        if (file.exists(file.path(output, "bf-import.rds"))) {
            bfpwr_sim_moment_check_identity(readRDS(file.path(output, "bf-import.rds"))$compatible_identity,
                bf_identity, "import proof")
        }
        files <- file.path(output, "case-results", sprintf("%04d.rds", seq_along(cases)))
        stopifnot(all(file.exists(files)))
        results <- lapply(files, readRDS)
        stopifnot(all(vapply(results, function(x) identical(x$identity, identity), logical(1))))
        combine <- function(name) do.call(rbind, lapply(seq_along(results), function(i) {
            x <- results[[i]][[name]]
            if (is.null(x) || !nrow(x)) return(NULL)
            x$case_id <- cases[[i]]$bf_prior_id
            x$alternative <- cases[[i]]$analysis_prior$alternative
            x
        }))
        names <- c("fixed", "decisions", "cumulative", "final", "sequential",
                   "moments", "searches", "timings", "reuse", "case_timings")
        bundle <- setNames(lapply(names, combine), names)
        bundle$timing_schema_version <- 2L
        bundle$cases <- cases
        bundle$designs <- unname(designs)
        bundle$provenance <- readRDS(provenance_path)
        bundle$provenance$bf_evaluation <- readRDS(file.path(output, "bf-evaluation-provenance.rds"))
        bundle$provenance$bf_rows <- sum(vapply(chunks, `[[`, numeric(1), "rows"))
        empirical <- file.path(output, "imported-empirical.rds")
        if (file.exists(empirical)) bfpwr_sim_moment_check_empirical(bundle, readRDS(empirical))
        bundle$fixed_diagnostics <- bfpwr_sim_mc_reference_diagnostics(bundle$fixed,
            label = "one-sided moment fixed", group_cols = c("alternative", "tail"))
        bundle$sequential_diagnostics <- bfpwr_sim_mc_reference_diagnostics(bundle$sequential,
            label = "one-sided moment sequential", group_cols = c("alternative", "event"))
        source(file.path(root, "simulations/R/moment-reference-check.R"), local = .GlobalEnv)
        bundle$numerical_diagnostics <- bfpwr_sim_moment_reference_check(bundle)
        saveRDS(bundle$numerical_diagnostics, file.path(output, "numerical-diagnostics.rds"), compress = "xz")
        write.csv(bundle$numerical_diagnostics$moment_rows,
                  file.path(output, "numerical-diagnostics.csv"), row.names = FALSE)
        saveRDS(bundle, file.path(output, "verification.rds"), compress = "xz")
        for (name in names) write.csv(bundle[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
        print(bundle$fixed_diagnostics$summary)
        print(bundle$sequential_diagnostics$summary)
        print(bundle$numerical_diagnostics$summary)
        cat("Stopping moment failures:", sum(!bundle$moments$EN_ok), sum(!bundle$moments$VarN_ok), "\n")
        cat("Search minimum mismatches:", sum(!bundle$searches$minimum_correct), "\n")
    } else stop("Unknown phase: ", phase)
    invisible(output)
}

if (sys.nframe() == 0L) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
    bfpwr_sim_run_moment_extension(commandArgs(TRUE), root)
}
