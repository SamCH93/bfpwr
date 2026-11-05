## Run from the repository root with --results pointing to a saved v2 bundle.
source("simulations/scripts/simulation_support.R")

check_export <- function(results) {
    script <- normalizePath("simulations/scripts/prepare_moment_extension_report.R",
        winslash = "/", mustWork = TRUE)
    scratch <- tempfile("moment-export-")
    stopifnot(dir.create(scratch))
    on.exit(unlink(scratch, recursive = TRUE))
    input <- file.path(scratch, "saved results")
    stopifnot(dir.create(input))
    stopifnot(all(file.copy(file.path(results, c("verification.rds", "numerical-diagnostics.rds")), input)))
    original_files <- list.files(input)
    original_hashes <- tools::sha256sum(file.path(input, original_files))

    run <- function(label, args, expect = 0L) {
        log <- file.path(scratch, "export.log")
        status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(script, args)),
            stdout = log, stderr = log)
        messages <- readLines(log)
        cat(label, "exit", status, "expected", expect, "\n")
        cat(messages, sep = "\n")
        stopifnot(status == expect,
            identical(list.files(input, all.files = TRUE, no.. = TRUE), original_files),
            identical(tools::sha256sum(file.path(input, original_files)), original_hashes))
        if (expect != 0L) {
            stopifnot(any(grepl("--output must be separate from the --results directory",
                messages, fixed = TRUE)))
        }
    }

    run("same path", c("--results", input, "--output", input), expect = 1L)
    old_wd <- getwd()
    setwd(scratch)
    on.exit(setwd(old_wd), add = TRUE, after = FALSE)
    run("relative alias", c("--results", input, "--output", "./saved results/"), expect = 1L)
    if (.Platform$OS.type == "windows") {
        run("case alias", c("--results", input, "--output", toupper(input)), expect = 1L)
    }
    run("default sibling", c("--results", paste0(input, "/")))
    explicit <- file.path(scratch, "explicit output")
    run("explicit output", c("--results", input, "--output", explicit))
    for (output in c(paste0(input, "-report"), explicit)) {
        stopifnot(all(file.exists(file.path(output,
            c("numerical-diagnostics.rds", "numerical-diagnostics.csv", "report-fixtures.rds")))),
            dir.exists(file.path(output, "fixture-validation")))
    }
    cat("Export destinations verified and private input bytes and directory entries unchanged.\n")
}

results <- arg_value(parse_args(), "results")
if (is.null(results)) stop("Supply --results with the saved v2 result directory.", call. = FALSE)
check_export(results)
