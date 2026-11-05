script_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    if (length(hit) > 0) {
        return(normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                             winslash = "/", mustWork = TRUE))
    }
    normalizePath("simulations/scripts/run_verification.R",
                  winslash = "/", mustWork = FALSE)
}

source(file.path(dirname(script_path()), "common.R"))

main <- function() {
    args <- parse_args()
    corpus_root <- normalizePath(arg_value(args, "corpus-root",
                                           "simulations/corpus/v1"),
                                 winslash = "/", mustWork = TRUE)
    scripts_dir <- dirname(script_path())
    fixture_script <- file.path(scripts_dir, "validate_fixture_suite.R")
    search_script <- file.path(scripts_dir, "refresh_search_validation.R")

    fixture_cmd <- c(fixture_script, "--corpus-root", corpus_root)
    fixture_set <- arg_value(args, "fixture-set", NULL)
    if (!is.null(fixture_set) && nzchar(fixture_set)) {
        fixture_cmd <- c(fixture_cmd, "--fixture-set", fixture_set)
    }
    search_cmd <- c(search_script, "--corpus-root", corpus_root,
                    "--scope", arg_value(args, "search-scope", "smoke"))
    search_bundle <- file.path(corpus_root, "search-validation",
                               "search-validation-comparison.rds")
    run_search <- !arg_flag(args, "skip-search") &&
        (file.exists(search_bundle) || arg_flag(args, "require-search"))

    cat("== simulation output checks ==\n")
    status <- system2(file.path(R.home("bin"), "Rscript"), fixture_cmd)
    if (!identical(status, 0L)) {
        stop("simulation output checks failed with exit code ", status,
             call. = FALSE)
    }

    if (run_search) {
        cat("\n== sample-size search checks ==\n")
        status <- system2(file.path(R.home("bin"), "Rscript"), search_cmd)
        if (!identical(status, 0L)) {
            stop("sample-size search checks failed with exit code ", status,
                 call. = FALSE)
        }
    } else {
        cat("\n== sample-size search checks skipped ==\n")
        if (arg_flag(args, "skip-search")) {
            cat("sample-size search checks were disabled by --skip-search\n")
        } else {
            cat("missing optional sample-size search bundle:",
                search_bundle, "\n")
        }
    }
}

main()
