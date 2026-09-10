script_path <- function(default = NULL) {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    if (length(hit) > 0) {
        return(normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                             winslash = "/", mustWork = TRUE))
    }
    if (is.null(default)) {
        default <- "simulations/scripts/simulation_support.R"
    }
    normalizePath(default, winslash = "/", mustWork = FALSE)
}

repo_root <- function() {
    normalizePath(file.path(dirname(script_path()), "..", ".."),
                  winslash = "/", mustWork = TRUE)
}

source_package_checkout <- function(root = repo_root()) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    r_dir <- file.path(root, "package", "R")
    files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
    if (length(files) == 0) {
        stop("no package R files found under ", r_dir, call. = FALSE)
    }
    for (file in files) {
        source(file, local = .GlobalEnv)
    }
    options(bfpwr.sim.package_checkout = root)
    invisible(files)
}

bfpwr_sim_package_provenance <- function(
        root = getOption("bfpwr.sim.package_checkout", repo_root())) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    git_output <- function(args) {
        result <- try(system2("git", c("-C", shQuote(root), args),
                              stdout = TRUE, stderr = FALSE), silent = TRUE)
        if (inherits(result, "try-error") ||
            !is.null(attr(result, "status"))) {
            return(character())
        }
        result
    }

    revision <- git_output(c("rev-parse", "HEAD"))
    package_status <- git_output(c("status", "--porcelain",
                                   "--untracked-files=normal", "--",
                                   "package"))
    description <- file.path(root, "package", "DESCRIPTION")
    version <- if (file.exists(description)) {
        read.dcf(description, fields = "Version")[[1]]
    } else {
        NA_character_
    }

    data.frame(
        package_git_revision = if (length(revision)) {
            revision[[1]]
        } else {
            NA_character_
        },
        package_git_dirty = if (length(revision)) {
            length(package_status) > 0
        } else {
            NA
        },
        package_version = version,
        package_source_root = root,
        package_recomputed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        stringsAsFactors = FALSE
    )
}

source_simulation_library <- function(root = repo_root()) {
    files <- c(
        "simulations/R/corpus-release.R",
        "simulations/R/grids.R",
        "simulations/R/cases.R",
        "simulations/registry/factories.R",
        "simulations/R/registry.R",
        "simulations/R/fixture-summaries.R",
        "simulations/R/fixture-validation.R",
        "simulations/R/search-validation.R",
        "simulations/R/package-verification.R",
        "simulations/registry/designs/z.R",
        "simulations/registry/designs/t.R",
        "simulations/registry/designs/binomial.R",
        "simulations/registry/analyses/z.R",
        "simulations/registry/analyses/t.R",
        "simulations/registry/analyses/binomial.R",
        "simulations/registry/bf-priors/z.R",
        "simulations/registry/bf-priors/t.R",
        "simulations/registry/bf-priors/binomial.R",
        "simulations/registry/design-cases.R",
        "simulations/registry/analysis-cases.R",
        "simulations/registry/bf-prior-cases.R",
        "simulations/registry/fixtures/thresholds.R",
        "simulations/registry/fixtures/z.R",
        "simulations/registry/fixtures/t.R",
        "simulations/registry/fixtures/binomial.R",
        "simulations/registry/fixture-cases.R"
    )
    for (file in file.path(root, files)) {
        source(file, local = .GlobalEnv)
    }
    invisible(files)
}

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

run_rscript <- function(script, args = character()) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c(script, args))
    if (!identical(status, 0L)) {
        stop("Rscript failed with exit code ", status, ": ", script,
             call. = FALSE)
    }
    invisible(status)
}
