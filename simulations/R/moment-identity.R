## Phase identities for the one-sided moment extension. No work runs on source().
bfpwr_sim_moment_text_sha256 <- function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    while (length(lines) && !nzchar(tail(lines, 1L))) lines <- head(lines, -1L)
    file <- tempfile()
    on.exit(unlink(file))
    writeBin(charToRaw(paste(lines, collapse = "\n")), file)
    unname(tools::sha256sum(file))
}

bfpwr_sim_moment_source_manifest <- function(root, files) {
    files <- sort(unique(gsub("\\\\", "/", files)), method = "radix")
    paths <- file.path(root, files)
    if (any(!file.exists(paths))) stop("missing calculation sources: ", paste(files[!file.exists(paths)], collapse = ", "))
    data.frame(file = files, sha256 = unname(tools::sha256sum(paths)),
        normalized_sha256 = vapply(paths, bfpwr_sim_moment_text_sha256, ""),
        stringsAsFactors = FALSE, row.names = NULL)
}

bfpwr_sim_moment_source_files <- function(root, phase) {
    if (phase == "bf") return(c("package/R/nmbf01.R", "simulations/R/bayes-factors.R"))
    c(paste0("package/", c("DESCRIPTION", "NAMESPACE", "LICENSE", "NEWS.md", "README.md", "Makefile")),
      paste0("package/R/", list.files(file.path(root, "package", "R"), "[.]R$")),
      paste0("package/man/", list.files(file.path(root, "package/man"), "[.]Rd$")),
      paste0("package/tests/", list.files(file.path(root, "package/tests"), "[.]R$", recursive = TRUE)),
      paste0("package/vignettes/", list.files(file.path(root, "package/vignettes"), "[.](Rnw|bib)$")),
      paste0("simulations/R/", list.files(file.path(root, "simulations", "R"), "[.]R$")),
      paste0("simulations/registry/", list.files(file.path(root, "simulations", "registry"),
          "[.]R$", recursive = TRUE)),
      "simulations/registry/manifests/package-verification-cases.csv",
      paste0("simulations/", list.files(file.path(root, "simulations"), "[.](Rnw|md)$")),
      paste0("simulations/scripts/", list.files(file.path(root, "simulations/scripts"), "[.](R|py)$")),
      paste0("simulations/jobs/", list.files(file.path(root, "simulations/jobs"), "[.]sh$")))
}

## Presentation and release tools remain in the raw snapshot, not the cache key.
bfpwr_sim_moment_dependency_files <- function(root, phase) {
    if (phase == "bf") return(bfpwr_sim_moment_source_files(root, phase))
    c("package/DESCRIPTION", "package/NAMESPACE",
      paste0("package/R/", list.files(file.path(root, "package/R"), "[.]R$")),
      paste0("simulations/R/", list.files(file.path(root, "simulations/R"), "[.]R$")),
      paste0("simulations/registry/", list.files(file.path(root, "simulations/registry"),
          "[.]R$", recursive = TRUE)),
      paste0("simulations/scripts/", c("run_moment_extension.R", "simulation_support.R",
          "common.R", "finalize_moment_extension_chunks.R", "validate_moment_extension.R")))
}

bfpwr_sim_moment_runtime <- function() {
    list(R = as.character(getRversion()), platform = R.version$platform,
         RNGkind = RNGkind(), packages = setNames(lapply(c("mvtnorm", "lamW", "withr", "qrng"),
             function(x) as.character(utils::packageVersion(x))), c("mvtnorm", "lamW", "withr", "qrng")),
         LAPACK = La_version(), BLAS = basename(extSoftVersion()[["BLAS"]]))
}

bfpwr_sim_moment_bf_runtime <- function() {
    list(R = as.character(getRversion()), platform = R.version$platform,
         stats = as.character(utils::packageVersion("stats")))
}

bfpwr_sim_moment_input_indexes <- function(input, designs) {
    setNames(lapply(designs, function(design) {
        design_path <- file.path(bfpwr_sim_design_dir(input, design), "design.rds")
        if (!file.exists(design_path) || !identical(readRDS(design_path), design)) {
            stop("original input design differs from the frozen design: ", design$design_case_id)
        }
        path <- file.path(bfpwr_sim_design_dir(input, design), "chunk-index.rds")
        if (!file.exists(path)) stop("missing original input index: ", path)
        index <- readRDS(path)
        required <- c("chunk_id", "replicate_start", "replicate_end", "n_rows", "bytes", "sha256", "file")
        stopifnot(all(required %in% names(index)), !anyNA(index$sha256),
                  identical(index$chunk_id, design$chunks$chunk_id))
        index[required]
    }), vapply(designs, `[[`, "", "design_case_id"))
}

bfpwr_sim_moment_identity <- function(root, phase, cases, designs, inputs) {
    sources <- bfpwr_sim_moment_source_manifest(root, bfpwr_sim_moment_dependency_files(root, phase))
    values <- list(cases = cases, designs = designs, inputs = inputs)
    if (phase == "bf") values$runtime <- bfpwr_sim_moment_bf_runtime()
    if (phase == "reference") values <- c(values, list(
        schedules = bfpwr_sim_z_fixture_schedules_core(),
        thresholds = bfpwr_sim_fixture_evidence_thresholds_core(),
        targets = bfpwr_sim_fixture_target_probabilities_core(), defaults = bfpwrOptions(),
        runtime = bfpwr_sim_moment_runtime()))
    list(schema_version = 2L, phase = phase,
         sources = sources[c("file", "normalized_sha256")], values = values)
}

bfpwr_sim_moment_check_identity <- function(saved, current, label) {
    if (!identical(saved, current)) stop("incompatible ", label,
        "; use a new output directory and recompute this phase", call. = FALSE)
    invisible(TRUE)
}

bfpwr_sim_moment_snapshot <- function(root, output, phase) {
    files <- bfpwr_sim_moment_source_files(root, phase)
    manifest <- bfpwr_sim_moment_source_manifest(root, files)
    snapshot <- file.path(output, "sources", phase)
    for (file in files) {
        target <- file.path(snapshot, file)
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        stopifnot(file.copy(file.path(root, file), target, overwrite = FALSE))
    }
    bfpwr_sim_write_rds(manifest, file.path(snapshot, "source-manifest.rds"))
    utils::write.csv(manifest, file.path(snapshot, "source-manifest.csv"), row.names = FALSE)
    manifest
}

bfpwr_sim_moment_verify_sources <- function(root, output, phase) {
    snapshot <- file.path(output, "sources", phase)
    saved <- readRDS(file.path(snapshot, "source-manifest.rds"))
    ## Raw hashes verify the physical source snapshot; normalized hashes decide compatibility.
    actual <- bfpwr_sim_moment_source_manifest(snapshot, saved$file)
    stopifnot(identical(saved$sha256, actual$sha256))
    dependencies <- bfpwr_sim_moment_dependency_files(snapshot, phase)
    saved <- saved[match(sort(dependencies, method = "radix"), saved$file), c("file", "normalized_sha256")]
    rownames(saved) <- NULL
    current <- bfpwr_sim_moment_source_manifest(root, bfpwr_sim_moment_dependency_files(root, phase))
    stopifnot(identical(saved, current[c("file", "normalized_sha256")]))
    current
}

bfpwr_sim_moment_revision <- function(root) {
    result <- suppressWarnings(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
                                       stdout = TRUE, stderr = FALSE))
    if (!length(result) || !is.null(attr(result, "status"))) {
        Sys.getenv("BFPWR_SIM_SOURCE_REVISION", NA_character_)
    } else result[[1]]
}

bfpwr_sim_moment_provenance <- function(root, input, manifest, cases, designs, inputs,
                                        phase = "reference") {
    list(schema_version = 2L, package_revision = bfpwr_sim_moment_revision(root),
         source_root = file.path("sources", phase), source_path_base = "corpus_output",
         source_checksums = setNames(manifest$sha256, manifest$file),
         normalized_source_checksums = setNames(manifest$normalized_sha256, manifest$file),
         defaults = bfpwrOptions(), session = sessionInfo(), runtime = bfpwr_sim_moment_runtime(),
         created = Sys.time(), source_corpus = input, input_path_base = "source_corpus",
         cases = cases, designs = designs, inputs = inputs,
         schedules = bfpwr_sim_z_fixture_schedules_core(),
         thresholds = bfpwr_sim_fixture_evidence_thresholds_core(),
         targets = bfpwr_sim_fixture_target_probabilities_core(),
         replicates = unique(vapply(designs, `[[`, integer(1), "nsim")))
}

bfpwr_sim_moment_verify_input <- function(input, design, chunk_id) {
    record <- bfpwr_sim_moment_input_indexes(input, list(design))[[1]]
    record <- record[record$chunk_id == chunk_id, , drop = FALSE]
    path <- bfpwr_sim_design_chunk_files(input, design)[match(chunk_id, design$chunks$chunk_id)]
    if (nrow(record) != 1L || !file.exists(path) ||
        !identical(unname(tools::sha256sum(path)), record$sha256[[1]]) ||
        file.info(path)$size != record$bytes[[1]]) stop("original input SHA256 mismatch: ", path)
    record
}

bfpwr_sim_moment_v1_catalog <- function() {
    names <- paste0("bfpwr-sim-corpus-v1-moment-bf-chunks-", c("greater", "less"), ".tar.zst")
    list(metadata_sha256 = "2c5906d4dd598074a75f9d641b1e8183a1b8b99418782d58efbdfd3a1557e4a9",
        assets = data.frame(name = names, asset_set = "bf-chunks",
            sha256 = c("4f23a0a6f9436d36c8bb6f6cefa5e500658e140b044c19551f16185277dffc9a",
                       "7c0bdcd2ae6433479716e3f5671f688eecc578cf52c715346bef696fe5434d7c"),
            bytes = c(796447610, 764644799), url = paste0(
                "https://github.com/FBartos/bfpwr/releases/download/sim-corpus-v1/", names),
            historical = TRUE, stringsAsFactors = FALSE))
}

bfpwr_sim_moment_archive_members <- function(archive) {
    members <- system2("tar", c("-tf", shQuote(archive)), stdout = TRUE)
    types <- system2("tar", c("-tvf", shQuote(archive)), stdout = TRUE)
    paths <- gsub("\\\\", "/", members)
    stopifnot(is.null(attr(members, "status")), is.null(attr(types, "status")),
        length(members) > 0L, length(types) == length(members),
        !anyDuplicated(paths),
        !any(grepl("^/|^[A-Za-z]:|(^|/)[.][.](/|$)", paths)))
    if (!all(substr(types, 1L, 1L) %in% c("-", "d"))) {
        stop("archive contains unsupported entry types; only regular files and directories are allowed")
    }
    members
}

## Only the immutable published metadata can bind an import to published BF assets.
bfpwr_sim_moment_official_bf_proof <- function(archive, chunks,
                                             catalog = bfpwr_sim_moment_v1_catalog()) {
    if (!identical(unname(tools::sha256sum(archive)), catalog$metadata_sha256)) {
        stop("official metadata archive SHA256 mismatch")
    }
    bfpwr_sim_moment_archive_members(archive)
    scratch <- tempfile("official-moment-metadata-")
    dir.create(scratch)
    on.exit(unlink(scratch, recursive = TRUE))
    stopifnot(system2("tar", c("-xf", shQuote(archive), "-C", shQuote(scratch))) == 0L)
    indexes <- list.files(file.path(scratch, "bf-priors"), "^bf-chunk-index[.]rds$",
                          recursive = TRUE, full.names = TRUE)
    if (!length(indexes)) stop("official metadata has no BF indexes")
    authoritative <- do.call(rbind, lapply(indexes, function(path) {
        index <- readRDS(path)
        stopifnot(all(c("bf_prior_id", "chunk_id", "source_design_chunk_sha256", "sha256") %in% names(index)))
        data.frame(case_id = index$bf_prior_id, chunk_id = index$chunk_id,
            input_sha256 = index$source_design_chunk_sha256, output_sha256 = index$sha256)
    }))
    order_chunks <- function(x) {
        x <- x[order(x$case_id, x$chunk_id, method = "radix"), c("case_id", "chunk_id", "input_sha256", "output_sha256")]
        rownames(x) <- NULL
        x
    }
    authoritative <- order_chunks(authoritative)
    chunks <- order_chunks(chunks)
    stopifnot(!anyNA(authoritative), !anyDuplicated(authoritative[c("case_id", "chunk_id")]))
    if (!identical(chunks, authoritative)) stop("imported BF payload does not match complete official metadata indexes")
    list(metadata_sha256 = catalog$metadata_sha256, chunks = authoritative, assets = catalog$assets)
}

## Import preserves historical identities; only the optional official proof authorizes archive links.
bfpwr_sim_moment_import_v1 <- function(source, output, input, root, bf_source_root,
                                      official_metadata = NULL) {
    if (file.exists(file.path(output, "bf-import.rds"))) stop("BF import already completed")
    if (length(list.files(file.path(output, "case-results")))) {
        stop("import requires a new destination without reference results")
    }
    cases <- readRDS(file.path(output, "cases.rds"))
    old_cases <- readRDS(file.path(source, "cases.rds"))
    stopifnot(identical(cases, old_cases))
    old_provenance <- readRDS(file.path(source, "bf-evaluation-provenance.rds"))
    old_identity <- readRDS(file.path(source, "bf-evaluation-identity.rds"))
    files <- bfpwr_sim_moment_source_files(root, "bf")
    old_files <- file.path(bf_source_root, files)
    stopifnot(all(file.exists(old_files)), all(files %in% names(old_provenance$source_checksums)),
              identical(unname(tools::md5sum(old_files)),
                        unname(old_provenance$source_checksums[files])))
    old_manifest <- bfpwr_sim_moment_source_manifest(bf_source_root, files)
    current <- bfpwr_sim_moment_verify_sources(root, output, "bf")
    if (!identical(old_manifest$normalized_sha256, current$normalized_sha256)) {
        stop("unknown BF source change; regenerate BF chunks")
    }
    for (file in files) {
        target <- file.path(output, "sources", "bf-historical", file)
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(target)) stopifnot(file.copy(file.path(bf_source_root, file), target))
        stopifnot(identical(unname(tools::sha256sum(target)), old_manifest$sha256[match(file, old_manifest$file)]))
    }
    bfpwr_sim_write_rds(old_manifest, file.path(output, "sources/bf-historical/source-manifest.rds"))
    plan <- bfpwr_sim_moment_extension_plan(cases)
    proof <- vector("list", nrow(plan))
    for (i in seq_len(nrow(plan))) {
        row <- plan[i, ]
        case <- cases[[row$case_index]]
        stopifnot(identical(case, readRDS(file.path(bfpwr_sim_bf_prior_dir(source, case), "bf-prior.rds"))))
        design <- bfpwr_sim_find_design_case(case$design_case_id)
        old_index <- readRDS(file.path(bfpwr_sim_bf_prior_dir(source, case), "bf-chunk-index.rds"))
        index <- old_index[old_index$chunk_id == row$chunk_id, , drop = FALSE]
        input_record <- bfpwr_sim_moment_verify_input(input, design, row$chunk_id)
        stopifnot(nrow(index) == 1L,
            identical(index$source_design_chunk_sha256[[1]], input_record$sha256[[1]]))
        old_path <- bfpwr_sim_bf_chunk_files(source, case, design)[row$chunk_id]
        stopifnot(file.exists(old_path),
                  identical(unname(tools::sha256sum(old_path)), index$sha256[[1]]))
        status <- readRDS(file.path(source, "chunk-status", sprintf("%04d.rds", i)))
        stopifnot(identical(status$identity, old_identity),
                  identical(unname(tools::md5sum(old_path)), status$output_md5),
                  identical(unname(tools::md5sum(bfpwr_sim_design_chunk_files(input, design)[row$chunk_id])),
                            status$input_md5))
        target <- bfpwr_sim_bf_chunk_files(output, case, design)[row$chunk_id]
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        if (!file.exists(target)) stopifnot(file.copy(old_path, target))
        stopifnot(identical(unname(tools::sha256sum(target)), index$sha256[[1]]))
        ## Keep the exact original identity/completion provenance, making paths portable.
        old_input <- status$input
        status$input <- file.path("designs", design$test_family, design$design_case_id, input_record$file)
        status$input_path_base <- "provenance$source_corpus"
        status$historical_input <- old_input
        status$historical_input_path_base <- if (grepl("^/|^[A-Za-z]:|^\\\\", old_input)) {
            "absolute historical path"
        } else "bf-evaluation-provenance$source_corpus"
        status$input_sha256 <- input_record$sha256[[1]]
        status$output_sha256 <- index$sha256[[1]]
        status$reused_from <- "moment-v1"
        bfpwr_sim_write_rds(status, file.path(output, "chunk-status", sprintf("%04d.rds", i)))
        proof[[i]] <- data.frame(case_id = case$bf_prior_id, chunk_id = row$chunk_id,
            input_sha256 = input_record$sha256[[1]], output_sha256 = index$sha256[[1]])
    }
    bfpwr_sim_write_rds(old_identity, file.path(output, "bf-evaluation-identity.rds"))
    bfpwr_sim_write_rds(old_provenance, file.path(output, "bf-evaluation-provenance.rds"))
    proof <- list(source_version = "moment-v1", sources = old_manifest,
        original_provenance = old_provenance, chunks = do.call(rbind, proof),
        compatible_identity = readRDS(file.path(output, "bf-dependency-identity.rds")),
        imported = Sys.time())
    if (!is.null(official_metadata)) {
        proof$official_payload <- bfpwr_sim_moment_official_bf_proof(official_metadata, proof$chunks)
    }
    bfpwr_sim_write_rds(proof, file.path(output, "bf-import.rds"))
    old_bundle <- readRDS(file.path(source, "verification.rds"))
    empirical <- lapply(old_bundle[c("fixed", "decisions", "cumulative", "final", "sequential", "moments")],
        function(x) x[setdiff(names(x), c("reference_prob", "reference_EN", "reference_VarN",
            "mcse_EN", "mcse_VarN", "EN_ok", "VarN_ok"))])
    bfpwr_sim_write_rds(empirical, file.path(output, "imported-empirical.rds"))
    invisible(proof)
}

bfpwr_sim_moment_check_empirical <- function(bundle, expected) {
    for (name in names(expected)) {
        x <- expected[[name]]
        y <- bundle[[name]][names(x)]
        key <- intersect(c("case_id", "schedule_id", "evidence_threshold", "n", "event", "tail"), names(x))
        order_rows <- function(z) do.call(order, c(z[key], list(method = "radix")))
        x <- x[order_rows(x), , drop = FALSE]
        y <- y[order_rows(y), , drop = FALSE]
        rownames(x) <- rownames(y) <- NULL
        if (!identical(x, y)) stop("import changed v1 empirical columns: ", name)
    }
    invisible(TRUE)
}
