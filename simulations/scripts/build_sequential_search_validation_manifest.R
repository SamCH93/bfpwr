support_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    script <- if (length(hit) > 0) {
        normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                      winslash = "/", mustWork = TRUE)
    } else {
        normalizePath("simulations/scripts/build_sequential_search_validation_manifest.R",
                      winslash = "/", mustWork = FALSE)
    }
    file.path(dirname(script), "simulation_support.R")
}

source(support_path())
source_simulation_library()

write_csv <- function(x, path) {
    ensure_dir(dirname(path))
    utils::write.csv(x, path, row.names = FALSE, na = "")
    invisible(path)
}

coverage_summary <- function(manifest) {
    split_key <- interaction(manifest$package_function,
                             manifest$manifest_role,
                             manifest$bf_type,
                             manifest$evidence,
                             drop = TRUE, sep = "\r")
    rows <- lapply(split(manifest, split_key), function(x) {
        data.frame(
            package_function = x$package_function[[1]],
            manifest_role = x$manifest_role[[1]],
            bf_type = x$bf_type[[1]],
            evidence = x$evidence[[1]],
            rows = nrow(x),
            achieved_rows = sum(x$sim_achieved %in% TRUE),
            unique_designs = length(unique(x$design_case_id)),
            unique_priors = length(unique(x$bf_prior_id)),
            unique_schedules = length(unique(x$schedule_id)),
            unique_thresholds = paste(sort(unique(x$evidence_threshold)),
                                      collapse = ","),
            unique_target_probs = paste(sort(unique(x$target_prob)),
                                        collapse = ","),
            unique_n_looks = paste(sort(unique(x$n_looks)), collapse = ","),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out[order(out$package_function, out$manifest_role, out$bf_type,
              out$evidence), , drop = FALSE]
}

main <- function() {
    args <- parse_args()
    root <- repo_root()
    corpus_root <- normalizePath(arg_value(args, "corpus-root",
                                           file.path(root, "simulations",
                                                     "corpus", "v1")),
                                 winslash = "/", mustWork = TRUE)
    output_file <- normalizePath(
        arg_value(args, "output-file",
                  bfpwr_sim_search_manifest_file(root)),
        winslash = "/", mustWork = FALSE)
    achieved_per_function <- as.integer(arg_value(args,
                                                  "achieved-per-function",
                                                  50L))
    diagnostic_per_function <- as.integer(arg_value(args,
                                                    "diagnostic-per-function",
                                                    10L))
    min_achieved_per_function <- as.integer(arg_value(
        args, "min-achieved-per-function", achieved_per_function))
    required_mcse_margin <- as.numeric(arg_value(
        args, "required-mcse-margin", 2))
    families <- trimws(strsplit(arg_value(args, "families", "z,t"),
                                ",", fixed = TRUE)[[1]])
    families <- families[nzchar(families)]

    manifest <- bfpwr_sim_build_sequential_search_manifest(
        corpus_root = corpus_root,
        families = families,
        achieved_per_function = achieved_per_function,
        diagnostic_per_function = diagnostic_per_function,
        min_achieved_per_function = min_achieved_per_function,
        required_mcse_margin = required_mcse_margin)
    write_csv(manifest, output_file)
    summary_file <- sub("[.]csv$", "-coverage.csv", output_file)
    write_csv(coverage_summary(manifest), summary_file)

    cat("wrote sequential search validation manifest:", output_file, "\n")
    cat("rows:", nrow(manifest), "\n")
    print(coverage_summary(manifest))
}

main()
