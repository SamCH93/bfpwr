script_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    if (length(hit) > 0) {
        return(normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                             winslash = "/", mustWork = TRUE))
    }
    normalizePath("simulations/scripts/compare_verification_timings.R",
                  winslash = "/", mustWork = FALSE)
}

source(file.path(dirname(script_path()), "common.R"))

read_fixture_timings <- function(path, label) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    files <- list.files(path, pattern = "-reference-timings[.]csv$",
                        full.names = TRUE)
    if (length(files) == 0) {
        stop(label, " fixture directory has no reference timing files: ",
             path, call. = FALSE)
    }
    rows <- lapply(files, function(file) {
        x <- utils::read.csv(file, stringsAsFactors = FALSE)
        x$timing_file <- basename(file)
        x
    })
    do.call(rbind, rows)
}

timing_key <- function(x) {
    ignored <- c("elapsed_seconds", "timing_status", "timing_error")
    columns <- setdiff(names(x), ignored)
    apply(x[columns], 1L, function(row) {
        paste(ifelse(is.na(row), "", as.character(row)), collapse = "\r")
    })
}

align_timings <- function(baseline, candidate) {
    baseline$key <- timing_key(baseline)
    candidate$key <- timing_key(candidate)
    if (anyDuplicated(baseline$key) || anyDuplicated(candidate$key)) {
        stop("timing rows do not have unique verification keys", call. = FALSE)
    }
    index <- match(baseline$key, candidate$key)
    if (anyNA(index) || nrow(baseline) != nrow(candidate)) {
        stop("baseline and candidate timing rows do not cover the same cases",
             call. = FALSE)
    }
    candidate <- candidate[index, , drop = FALSE]
    data.frame(
        package_function = baseline$package_function,
        baseline_seconds = baseline$elapsed_seconds,
        candidate_seconds = candidate$elapsed_seconds,
        stringsAsFactors = FALSE
    )
}

summarize_timings <- function(x, max_ratio, min_seconds) {
    groups <- split(x, x$package_function)
    out <- do.call(rbind, lapply(groups, function(group) {
        baseline_total <- sum(group$baseline_seconds)
        candidate_total <- sum(group$candidate_seconds)
        ratio <- candidate_total / baseline_total
        data.frame(
            package_function = group$package_function[[1]],
            calls = nrow(group),
            baseline_total_seconds = baseline_total,
            candidate_total_seconds = candidate_total,
            candidate_to_baseline_ratio = ratio,
            median_paired_delta_seconds = median(
                group$candidate_seconds - group$baseline_seconds),
            regression = baseline_total >= min_seconds && ratio > max_ratio,
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out[order(out$package_function), , drop = FALSE]
}

main <- function() {
    args <- parse_args()
    baseline_dir <- arg_value(args, "baseline-fixture-dir", NULL)
    candidate_dir <- arg_value(args, "candidate-fixture-dir", NULL)
    if (is.null(baseline_dir) || is.null(candidate_dir)) {
        stop("--baseline-fixture-dir and --candidate-fixture-dir are required",
             call. = FALSE)
    }
    max_ratio <- as.numeric(arg_value(args, "max-total-ratio", "1.10"))
    min_seconds <- as.numeric(arg_value(args, "min-baseline-seconds", "1"))
    if (!is.finite(max_ratio) || max_ratio <= 0 ||
        !is.finite(min_seconds) || min_seconds < 0) {
        stop("invalid timing comparison thresholds", call. = FALSE)
    }

    baseline <- read_fixture_timings(baseline_dir, "baseline")
    candidate <- read_fixture_timings(candidate_dir, "candidate")
    paired <- align_timings(baseline, candidate)
    summary <- summarize_timings(paired, max_ratio = max_ratio,
                                 min_seconds = min_seconds)
    print(summary, row.names = FALSE)

    output_file <- arg_value(args, "output-file", NULL)
    if (!is.null(output_file)) {
        utils::write.csv(summary, output_file, row.names = FALSE)
    }
    if (any(summary$regression)) {
        stop("verification timing regression for: ",
             paste(summary$package_function[summary$regression],
                   collapse = ", "), call. = FALSE)
    }
}

main()
