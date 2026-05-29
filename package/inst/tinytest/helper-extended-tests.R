## Source this file with local = TRUE so exit_file() resolves to tinytest's
## masked file-stopping helper in each test environment.

bfpwr_run_extended_tests <- function() {
    tolower(Sys.getenv("BFPWR_RUN_EXTENDED_TESTS", unset = "")) %in%
        c("true", "1", "yes")
}

bfpwr_extended_skip_message <- function(reason = "extended test") {
    paste0(
        reason,
        "; set BFPWR_RUN_EXTENDED_TESTS=true to run outside CRAN checks"
    )
}

bfpwr_timing_multiplier <- function() {
    multiplier <- suppressWarnings(as.numeric(
        Sys.getenv("BFPWR_TIMING_MULTIPLIER", "1")
    ))
    if (!is.finite(multiplier) || multiplier <= 0) {
        return(1)
    }
    multiplier
}

bfpwr_expect_elapsed_under <- function(label, seconds, expr) {
    limit <- seconds * bfpwr_timing_multiplier()
    timing <- system.time(value <- force(expr))
    elapsed <- unname(timing[["elapsed"]])
    expect_true(
        elapsed <= limit,
        info = sprintf(
            paste0(
                "%s took %.2f seconds; expected <= %.2f seconds. ",
                "This suggests a substantial performance regression; ",
                "use BFPWR_TIMING_MULTIPLIER only for known slow machines."
            ),
            label, elapsed, limit
        )
    )
    value
}
