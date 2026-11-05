## Versioned append-only release packaging; no work runs when sourced.
bfpwr_sim_package_moment_extension <- function(args, root) {
    source(file.path(root, "simulations/scripts/simulation_support.R"), local = .GlobalEnv)
    source(file.path(root, "simulations/R/moment-identity.R"), local = .GlobalEnv)
    args <- parse_args(args)
    output <- normalizePath(arg_value(args, "output"), winslash = "/", mustWork = TRUE)
    release <- arg_value(args, "release")
    version <- arg_value(args, "version", "moment-v2")
    stopifnot(version == "moment-v2", !is.null(release))
    ensure_dir(release)
    release <- normalizePath(release, winslash = "/", mustWork = TRUE)
    prefix <- paste0("bfpwr-sim-corpus-v1-", version)
    required <- c("cases.rds", "designs.rds", "input-indexes.rds", "chunk-plan.csv",
        "calculation-identity.rds", "reference-identity.rds", "bf-dependency-identity.rds",
        "bf-evaluation-identity.rds", "provenance.rds", "bf-evaluation-provenance.rds",
        "chunk-status", "case-results", "verification.rds", "fixed.csv", "sequential.csv",
        "moments.csv", "searches.csv", "timings.csv", "reuse.csv", "case_timings.csv",
        "numerical-diagnostics.rds", "numerical-diagnostics.csv", "validation-checks.csv",
        "numerical-validation-checks.csv")
    stopifnot(all(file.exists(file.path(output, required))),
              dir.exists(file.path(output, "sources/reference")),
              dir.exists(file.path(output, "sources/bf")),
              dir.exists(file.path(output, "index-status")))
    cases <- readRDS(file.path(output, "cases.rds"))
    indexes <- c(file.path("index-status", sprintf("%04d.rds", seq_along(cases))),
        unlist(lapply(cases, function(case) file.path("bf-priors", case$test_family,
            case$bf_prior_id, c("bf-prior.rds", "bf-chunk-index.rds"))), use.names = FALSE))
    if (!all(file.exists(file.path(output, indexes)))) stop("missing per-case BF metadata indexes")
    checks <- read.csv(file.path(output, "validation-checks.csv"))
    if (any(checks$failures > 0)) stop("MC/search validation failures must be resolved before packaging")
    numerical_checks <- read.csv(file.path(output, "numerical-validation-checks.csv"))
    if (any(numerical_checks$failures > 0)) stop("independent numerical reference did not converge")
    for (phase in c("bf", "reference", if (file.exists(file.path(output, "bf-import.rds"))) "bf-historical")) {
        snapshot <- file.path(output, "sources", phase)
        sources <- readRDS(file.path(snapshot, "source-manifest.rds"))
        stopifnot(identical(unname(tools::sha256sum(file.path(snapshot, sources$file))), sources$sha256))
    }
    base <- arg_value(args, "release-base", "https://github.com/FBartos/bfpwr/releases/download/sim-corpus-v1/")
    if (!endsWith(base, "/")) base <- paste0(base, "/")
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(output)
    payload <- sort(list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE), method = "radix")
    payload <- payload[!grepl("(^|/)[.]|[.](lock|tmp|error)(/|$)", payload)]
    identity <- data.frame(file = payload, sha256 = unname(tools::sha256sum(payload)))
    plan <- read.csv("chunk-plan.csv")
    designs <- readRDS("designs.rds")
    inputs <- readRDS("input-indexes.rds")
    expected_plan <- do.call(rbind, lapply(seq_along(cases), function(i) {
        case <- cases[[i]]
        design <- designs[[case$design_case_id]]
        data.frame(case_index = i, bf_prior_id = case$bf_prior_id,
            design_case_id = design$design_case_id, chunk_id = design$chunks$chunk_id)
    }))
    if (!nrow(plan) || !identical(plan, expected_plan) ||
        anyDuplicated(vapply(cases, `[[`, "", "bf_prior_id"))) {
        stop("saved BF chunk plan does not match the complete case and design keys")
    }
    bf_indexes <- lapply(seq_along(cases), function(i) {
        case <- cases[[i]]
        directory <- file.path("bf-priors", case$test_family, case$bf_prior_id)
        index <- readRDS(file.path(directory, "bf-chunk-index.rds"))
        chunk_ids <- designs[[case$design_case_id]]$chunks$chunk_id
        input_index <- inputs[[case$design_case_id]]
        index_checks <- c(
            prior_metadata_matches_case = identical(readRDS(file.path(directory, "bf-prior.rds")), case),
            index_status_matches_chunk_index = identical(
                readRDS(file.path("index-status", sprintf("%04d.rds", i))), index),
            index_chunk_ids_match_design = identical(index$chunk_id, chunk_ids),
            input_chunk_ids_match_design = identical(input_index$chunk_id, chunk_ids),
            index_prior_ids_match_case = identical(index$bf_prior_id, rep(case$bf_prior_id, length(chunk_ids))),
            index_design_ids_match_case = identical(index$design_case_id, rep(case$design_case_id, length(chunk_ids))),
            index_test_families_match_case = identical(index$test_family, rep(case$test_family, length(chunk_ids))),
            index_paths_match_plan = identical(index$file, file.path("bf_chunks", sprintf("chunk-%04d.rds", chunk_ids))),
            input_hashes_present = !anyNA(input_index$sha256),
            index_input_hashes_match = identical(index$source_design_chunk_sha256, input_index$sha256))
        if (!all(index_checks)) {
            stop("Invalid BF metadata indexes for case ", case$bf_prior_id, ": ",
                paste(names(index_checks)[!index_checks], collapse = ", "), call. = FALSE)
        }
        index
    })
    planned_paths <- vapply(seq_len(nrow(plan)), function(i) {
        case <- cases[[plan$case_index[i]]]
        file.path("bf-priors", case$test_family, case$bf_prior_id,
                  "bf_chunks", sprintf("chunk-%04d.rds", plan$chunk_id[i]))
    }, "")
    if (!setequal(planned_paths, payload[grepl("(^|/)bf_chunks/", payload)])) {
        stop("saved BF chunk files do not match the complete plan")
    }
    actual <- do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
        case <- cases[[plan$case_index[i]]]
        index <- bf_indexes[[plan$case_index[i]]]
        chunk <- match(plan$chunk_id[i], index$chunk_id)
        path <- planned_paths[i]
        status <- readRDS(file.path("chunk-status", sprintf("%04d.rds", i)))
        sha256 <- identity$sha256[match(path, identity$file)]
        if (is.na(sha256) || !identical(sha256, status$output_sha256) ||
            !identical(sha256, index$sha256[[chunk]]) ||
            !identical(status$input_sha256, index$source_design_chunk_sha256[[chunk]])) {
            stop("BF chunk SHA256 mismatch: ", path)
        }
        data.frame(case_id = case$bf_prior_id, chunk_id = plan$chunk_id[i],
            input_sha256 = status$input_sha256, output_sha256 = sha256)
    }))
    historical_assets <- NULL
    if (file.exists("bf-import.rds")) {
        proof <- readRDS("bf-import.rds")
        if (!identical(actual, proof$chunks)) stop("current BF payload differs from import proof")
        if (!is.null(proof$official_payload)) {
            official <- proof$official_payload
            catalog <- bfpwr_sim_moment_v1_catalog()
            actual <- actual[order(actual$case_id, actual$chunk_id, method = "radix"), ]
            rownames(actual) <- NULL
            stopifnot(identical(official$metadata_sha256, catalog$metadata_sha256),
                identical(official$assets, catalog$assets), identical(actual, official$chunks))
            historical_assets <- catalog$assets
        }
    }
    identity_file <- file.path(release, paste0(version, "-packaging-identity.rds"))
    if (file.exists(identity_file)) {
        if (!identical(readRDS(identity_file), identity)) stop("release destination already has a different payload")
    } else saveRDS(identity, identity_file)
    assets <- list()
    old_manifest_file <- file.path(release, paste0(version, "-manifest.csv"))
    old_manifest <- if (file.exists(old_manifest_file)) read.csv(old_manifest_file) else NULL
    archive <- function(suffix, files, exclude_bf = FALSE, set = suffix) {
        selected <- payload %in% files
        for (directory in files[dir.exists(files)]) {
            selected <- selected | startsWith(payload, paste0(directory, "/"))
        }
        files <- payload[selected]
        if (exclude_bf) files <- files[!grepl("(^|/)bf_chunks/", files)]
        stopifnot(length(files) > 0L)
        name <- paste0(prefix, "-", suffix, ".tar.zst")
        path <- file.path(release, name)
        if (file.exists(path)) {
            previous <- if (is.null(old_manifest)) NULL else old_manifest[old_manifest$name == name, , drop = FALSE]
            expected <- if (!is.null(previous) && nrow(previous) == 1L) previous$sha256[[1]] else {
                if (file.exists(paste0(path, ".sha256"))) readLines(paste0(path, ".sha256"))[1L] else NA_character_
            }
            if (!identical(unname(tools::sha256sum(path)), expected)) {
                stop("existing release asset has no matching verified manifest: ", name)
            }
            bfpwr_sim_moment_archive_members(path)
        } else {
            list_file <- tempfile()
            temporary <- paste0(path, ".tmp")
            writeLines(files, list_file)
            on.exit(unlink(c(list_file, temporary)), add = TRUE)
            args <- c("--zstd", "-cf", shQuote(temporary), "-T", shQuote(list_file))
            stopifnot(system2("tar", args) == 0L)
            bfpwr_sim_moment_archive_members(temporary)
            stopifnot(file.rename(temporary, path))
            writeLines(unname(tools::sha256sum(path)), paste0(path, ".sha256"))
        }
        assets[[length(assets) + 1L]] <<- data.frame(name = name, asset_set = set,
            sha256 = unname(tools::sha256sum(path)), bytes = file.info(path)$size,
            url = paste0(base, name), historical = FALSE, stringsAsFactors = FALSE)
    }
    optional <- c("bf-import.rds", "imported-empirical.rds", "known-search-limitation.csv")
    archive("results", c(required, optional[file.exists(optional)]))
    archive("metadata", c("bf-priors", "index-status"), exclude_bf = TRUE)
    archive("sources", "sources")
    if (!is.null(historical_assets)) {
        ## Unchanged v1 archives are explicit immutable references, not duplicated uploads.
        assets[[length(assets) + 1L]] <- historical_assets
    } else for (direction in unique(vapply(cases,
        function(x) x$analysis_prior$alternative, ""))) {
        directories <- list.dirs("bf-priors/z", recursive = FALSE)
        directories <- directories[grepl(paste0("-", direction, "-on-"), directories, fixed = TRUE)]
        stopifnot(length(directories) > 0L)
        archive(paste0("bf-chunks-", direction), directories, set = "bf-chunks")
    }
    manifest <- do.call(rbind, assets)
    manifest$version <- version
    manifest_file <- file.path(release, paste0(version, "-manifest.csv"))
    utils::write.csv(manifest, manifest_file, row.names = FALSE)
    paths <- c(file.path(release, manifest$name[!manifest$historical]), manifest_file, identity_file)
    writeLines(paste(unname(tools::sha256sum(paths)), basename(paths)),
               file.path(release, paste0("SHA256SUMS-", version, ".txt")))
    cat("Packaged", sum(!manifest$historical), "new archives and",
        sum(manifest$historical), "historical BF archive references at", release, "\n")
    invisible(manifest)
}

if (sys.nframe() == 0L) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
    bfpwr_sim_package_moment_extension(commandArgs(TRUE), root)
}
