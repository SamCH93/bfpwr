## Source this file with local = TRUE so exit_file() resolves to tinytest's
## masked file-stopping helper in each test environment.

bfpwr_run_extended_tests <- function() {
    tolower(Sys.getenv("BFPWR_RUN_EXTENDED_TESTS")) %in% c("true", "1", "yes")
}

bfpwr_exit_if_not_extended <- function(reason = "extended test") {
    if (!bfpwr_run_extended_tests()) {
        ## Keep slow/reference checks out of CRAN by default. Maintainers can
        ## opt in locally with BFPWR_RUN_EXTENDED_TESTS=true.
        ## tinytest masks exit_file() inside test files; keep this unqualified
        ## so the masked version stops the current file.
        exit_file(paste0(
            reason,
            "; set BFPWR_RUN_EXTENDED_TESTS=true to run outside CRAN checks"
        ))
    }
}
