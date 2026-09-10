## Shared helpers for the simulation verification scripts.
## These scripts validate the materialized fixture bundle used by the PR review;
## they do not regenerate the multi-GB Monte Carlo corpus from raw simulation
## chunks.

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop("unexpected positional argument: ", key, call. = FALSE)
        }
        key <- sub("^--", "", key)
        if (grepl("=", key, fixed = TRUE)) {
            parts <- strsplit(key, "=", fixed = TRUE)[[1]]
            out[[parts[[1]]]] <- parts[[2]]
            i <- i + 1
        } else if (i == length(args) || startsWith(args[[i + 1]], "--")) {
            out[[key]] <- TRUE
            i <- i + 1
        } else {
            out[[key]] <- args[[i + 1]]
            i <- i + 2
        }
    }
    out
}
arg_value <- function(args, name, default = NULL) {
    if (!is.null(args[[name]])) args[[name]] else default
}

arg_flag <- function(args, name) {
    isTRUE(args[[name]])
}

ensure_dir <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(path)
}

empty_failures <- function() {
    data.frame(
        severity = character(),
        rule = character(),
        n_failed = integer(),
        details = character(),
        stringsAsFactors = FALSE
    )
}

add_failure <- function(failures, severity, rule, n_failed, details) {
    if (is.na(n_failed) || n_failed <= 0) {
        return(failures)
    }
    rbind(
        failures,
        data.frame(
            severity = severity,
            rule = rule,
            n_failed = as.integer(n_failed),
            details = details,
            stringsAsFactors = FALSE
        )
    )
}

write_csv <- function(x, path) {
    ensure_dir(dirname(path))
    utils::write.csv(x, file = path, row.names = FALSE, na = "")
    invisible(path)
}

read_csv_if_exists <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
}

read_rds_required <- function(path) {
    if (!file.exists(path)) {
        stop("required file is missing: ", path, call. = FALSE)
    }
    readRDS(path)
}

safe_sum <- function(x, na.rm = TRUE) {
    if (length(x) == 0) {
        return(0)
    }
    sum(x, na.rm = na.rm)
}

count_bad_prob <- function(x) {
    safe_sum(!is.finite(x) | x < -1e-12 | x > 1 + 1e-12)
}

count_not_close <- function(x, target, tolerance = 1e-10) {
    safe_sum(!is.finite(x) | abs(x - target) > tolerance)
}

find_fixture_dirs <- function(corpus_root, fixture_set = NULL) {
    root <- file.path(corpus_root, "fixtures")
    if (!dir.exists(root)) {
        stop("fixture root is missing: ", root, call. = FALSE)
    }
    specs <- list.files(root, pattern = "^spec\\.rds$", recursive = TRUE,
                        full.names = TRUE)
    dirs <- dirname(specs)
    if (!is.null(fixture_set) && nzchar(fixture_set)) {
        requested <- strsplit(fixture_set, ",", fixed = TRUE)[[1]]
        requested <- trimws(requested[nzchar(trimws(requested))])
        dirs <- dirs[basename(dirs) %in% requested]
        missing <- setdiff(requested, basename(dirs))
        if (length(missing) > 0) {
            stop("requested fixture set(s) not found: ",
                 paste(missing, collapse = ", "), call. = FALSE)
        }
    }
    dirs[order(basename(dirs))]
}

fixture_role <- function(fixture_set_id, family, mode) {
    if (identical(fixture_set_id, "t-tbf01-fixed-core-v1")) {
        return("diagnostic_approximation")
    }
    if (identical(family, "binomial") && identical(mode, "sequential")) {
        return("integrity_only")
    }
    "package_reference"
}

read_overview <- function(fixture_dir) {
    path <- file.path(fixture_dir, "overview.csv")
    if (file.exists(path)) {
        x <- utils::read.csv(path, stringsAsFactors = FALSE)
        if (nrow(x) > 0) {
            return(x[1, , drop = FALSE])
        }
    }

    spec <- readRDS(file.path(fixture_dir, "spec.rds"))
    data.frame(
        fixture_set_id = basename(fixture_dir),
        family = if (!is.null(spec$family)) spec$family else NA_character_,
        mode = if (!is.null(spec$mode)) spec$mode else NA_character_,
        stringsAsFactors = FALSE
    )
}
