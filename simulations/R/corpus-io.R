bfpwr_sim_design_dir <- function(corpus_root, design_or_id, test_family = NULL) {
    design_case_id <- if (is.list(design_or_id)) design_or_id$design_case_id else design_or_id
    family <- if (is.list(design_or_id)) design_or_id$test_family else test_family
    if (is.null(family)) stop("test_family is required when design_or_id is not a design case")
    file.path(corpus_root, "designs", family, design_case_id)
}

bfpwr_sim_analysis_dir <- function(corpus_root, analysis_or_id, test_family = NULL) {
    analysis_case_id <- if (is.list(analysis_or_id)) analysis_or_id$analysis_case_id else analysis_or_id
    family <- if (is.list(analysis_or_id)) analysis_or_id$test_family else test_family
    if (is.null(family)) stop("test_family is required when analysis_or_id is not an analysis case")
    file.path(corpus_root, "analyses", family, analysis_case_id)
}

bfpwr_sim_bf_prior_dir <- function(corpus_root, bf_prior_or_id,
                                   test_family = NULL) {
    bf_prior_id <- if (is.list(bf_prior_or_id)) {
        bf_prior_or_id$bf_prior_id
    } else {
        bf_prior_or_id
    }
    family <- if (is.list(bf_prior_or_id)) {
        bf_prior_or_id$test_family
    } else {
        test_family
    }
    if (is.null(family)) {
        stop("test_family is required when bf_prior_or_id is not a BF prior case")
    }
    file.path(corpus_root, "bf-priors", family, bf_prior_id)
}

bfpwr_sim_write_rds <- function(object, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(object, path, compress = "xz")
    invisible(path)
}

bfpwr_sim_write_hashes <- function(paths, output_file) {
    hashes <- tools::sha256sum(paths)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste(names(hashes), hashes, sep = "  "), output_file)
    invisible(hashes)
}

bfpwr_sim_design_chunk_files <- function(corpus_root, design) {
    file.path(bfpwr_sim_design_dir(corpus_root, design), "chunks",
              sprintf("chunk-%04d.rds", design$chunks$chunk_id))
}

bfpwr_sim_read_design_chunk <- function(corpus_root, design, chunk_id) {
    file <- file.path(bfpwr_sim_design_dir(corpus_root, design), "chunks",
                      sprintf("chunk-%04d.rds", chunk_id))
    if (!file.exists(file)) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = file,
            asset_set = "chunk-validation",
            label = "design chunk files")
    }
    readRDS(file)
}

bfpwr_sim_bf_chunk_files <- function(corpus_root, bf_prior, design) {
    file.path(bfpwr_sim_bf_prior_dir(corpus_root, bf_prior), "bf_chunks",
              sprintf("chunk-%04d.rds", design$chunks$chunk_id))
}

bfpwr_sim_read_bf_prior_chunk <- function(corpus_root, bf_prior, design,
                                          chunk_id) {
    file <- file.path(bfpwr_sim_bf_prior_dir(corpus_root, bf_prior),
                      "bf_chunks", sprintf("chunk-%04d.rds", chunk_id))
    if (!file.exists(file)) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = file,
            asset_set = "chunk-validation",
            label = "BF chunk files")
    }
    readRDS(file)
}

bfpwr_sim_read_design_trajectories <- function(corpus_root, design) {
    design_dir <- bfpwr_sim_design_dir(corpus_root, design)
    trajectories_file <- file.path(design_dir, "trajectories.rds")
    if (file.exists(trajectories_file)) {
        return(readRDS(trajectories_file))
    }

    chunk_files <- bfpwr_sim_design_chunk_files(corpus_root, design)
    missing <- chunk_files[!file.exists(chunk_files)]
    if (length(missing) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing,
            asset_set = "chunk-validation",
            label = "design chunk files")
    }
    do.call(rbind, lapply(chunk_files, readRDS))
}

bfpwr_sim_write_design_chunk_index <- function(corpus_root, design) {
    design_dir <- bfpwr_sim_design_dir(corpus_root, design)
    chunk_files <- bfpwr_sim_design_chunk_files(corpus_root, design)
    missing <- chunk_files[!file.exists(chunk_files)]
    if (length(missing) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing,
            asset_set = "chunk-validation",
            label = "design chunk files")
    }

    rows <- vector("list", length(chunk_files))
    for (i in seq_along(chunk_files)) {
        chunk_id <- design$chunks$chunk_id[[i]]
        trajectories <- readRDS(chunk_files[[i]])
        bfpwr_sim_validate_design_trajectories(design, trajectories)

        expected <- seq(design$chunks$replicate_start[[i]],
                        design$chunks$replicate_end[[i]])
        observed <- sort(unique(trajectories$replicate_id))
        if (!identical(observed, expected)) {
            stop("chunk ", chunk_id, " does not contain the expected replicate range")
        }
        if (!all(trajectories$chunk_id == chunk_id)) {
            stop("chunk ", chunk_id, " contains unexpected chunk_id values")
        }

        info <- file.info(chunk_files[[i]])
        rows[[i]] <- data.frame(
            design_case_id = design$design_case_id,
            test_family = design$test_family,
            chunk_id = chunk_id,
            replicate_start = design$chunks$replicate_start[[i]],
            replicate_end = design$chunks$replicate_end[[i]],
            n_replicates = length(expected),
            n_looks = length(design$look_grid),
            n_rows = nrow(trajectories),
            bytes = info$size,
            sha256 = unname(tools::sha256sum(chunk_files[[i]])),
            file = file.path("chunks", basename(chunk_files[[i]])),
            stringsAsFactors = FALSE
        )
    }

    index <- do.call(rbind, rows)
    index_file <- file.path(design_dir, "chunk-index.rds")
    index_csv <- file.path(design_dir, "chunk-index.csv")
    bfpwr_sim_write_rds(index, index_file)
    utils::write.csv(index, index_csv, row.names = FALSE)
    bfpwr_sim_write_hashes(c(file.path(design_dir, "design.rds"),
                             index_file, index_csv, chunk_files),
                           file.path(design_dir, "sha256.txt"))
    index
}

bfpwr_sim_write_bf_prior_chunk_index <- function(corpus_root, bf_prior,
                                                 design, input_root = corpus_root) {
    bf_dir <- bfpwr_sim_bf_prior_dir(corpus_root, bf_prior)
    chunk_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior, design)
    missing <- chunk_files[!file.exists(chunk_files)]
    if (length(missing) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing,
            asset_set = "chunk-validation",
            label = "BF chunk files")
    }

    design_index_file <- file.path(bfpwr_sim_design_dir(input_root, design),
                                   "chunk-index.rds")
    design_index <- if (file.exists(design_index_file)) {
        readRDS(design_index_file)
    } else {
        data.frame()
    }

    rows <- vector("list", length(chunk_files))
    for (i in seq_along(chunk_files)) {
        chunk_id <- design$chunks$chunk_id[[i]]
        bf_rows <- readRDS(chunk_files[[i]])
        bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)

        expected <- seq(design$chunks$replicate_start[[i]],
                        design$chunks$replicate_end[[i]])
        observed <- sort(unique(bf_rows$replicate_id))
        if (!identical(observed, expected)) {
            stop("BF chunk ", chunk_id,
                 " does not contain the expected replicate range")
        }
        if (!all(bf_rows$chunk_id == chunk_id)) {
            stop("BF chunk ", chunk_id,
                 " contains unexpected chunk_id values")
        }

        info <- file.info(chunk_files[[i]])
        source <- if (nrow(design_index) > 0) {
            design_index[match(chunk_id, design_index$chunk_id), ,
                         drop = FALSE]
        } else {
            data.frame()
        }
        source_file <- if (nrow(source) == 1 && "file" %in% names(source)) {
            as.character(source$file[[1]])
        } else {
            NA_character_
        }
        source_sha256 <- if (nrow(source) == 1 && "sha256" %in% names(source)) {
            as.character(source$sha256[[1]])
        } else {
            NA_character_
        }
        source_bytes <- if (nrow(source) == 1 && "bytes" %in% names(source)) {
            source$bytes[[1]]
        } else {
            NA_real_
        }
        if (nrow(source) == 1 && !is.na(source_sha256)) {
            source_path <- file.path(bfpwr_sim_design_dir(input_root, design), source_file)
            if (!file.exists(source_path) ||
                !identical(unname(tools::sha256sum(source_path)), source_sha256) ||
                file.info(source_path)$size != source_bytes) {
                stop("source design chunk SHA256 mismatch: ", source_path)
            }
        }
        rows[[i]] <- data.frame(
            bf_prior_id = bf_prior$bf_prior_id,
            design_case_id = design$design_case_id,
            test_family = bf_prior$test_family,
            chunk_id = chunk_id,
            replicate_start = design$chunks$replicate_start[[i]],
            replicate_end = design$chunks$replicate_end[[i]],
            n_replicates = length(expected),
            n_looks = length(design$look_grid),
            n_rows = nrow(bf_rows),
            n_finite_log_bf01 = sum(is.finite(bf_rows$log_bf01)),
            n_infinite_log_bf01 = sum(is.infinite(bf_rows$log_bf01)),
            n_nan_log_bf01 = sum(is.nan(bf_rows$log_bf01)),
            n_na_log_bf01 = sum(is.na(bf_rows$log_bf01) &
                                !is.nan(bf_rows$log_bf01)),
            bytes = info$size,
            sha256 = unname(tools::sha256sum(chunk_files[[i]])),
            file = file.path("bf_chunks", basename(chunk_files[[i]])),
            source_design_chunk_file = source_file,
            source_design_chunk_sha256 = source_sha256,
            source_design_chunk_bytes = source_bytes,
            stringsAsFactors = FALSE
        )
    }

    index <- do.call(rbind, rows)
    index_file <- file.path(bf_dir, "bf-chunk-index.rds")
    index_csv <- file.path(bf_dir, "bf-chunk-index.csv")
    bfpwr_sim_write_rds(index, index_file)
    utils::write.csv(index, index_csv, row.names = FALSE)
    bfpwr_sim_write_hashes(c(file.path(bf_dir, "bf-prior.rds"),
                             index_file, index_csv, chunk_files),
                           file.path(bf_dir, "sha256.txt"))
    index
}

bfpwr_sim_source_package <- function(package_dir = file.path(getwd(), "package")) {
    r_files <- list.files(file.path(package_dir, "R"), pattern = "[.]R$",
                          full.names = TRUE)
    if (length(r_files) == 0) {
        stop("no package R files found under ", package_dir)
    }
    for (file in r_files) {
        source(file, local = .GlobalEnv)
    }
    invisible(r_files)
}

bfpwr_sim_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop("unexpected positional argument: ", key)
        }
        key <- sub("^--", "", key)
        if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
            out[[key]] <- TRUE
            i <- i + 1L
        } else {
            out[[key]] <- args[[i + 1L]]
            i <- i + 2L
        }
    }
    out
}
