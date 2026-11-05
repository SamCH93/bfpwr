bfpwr_sim_corpus_release_repo <- function() {
    "FBartos/bfpwr"
}

bfpwr_sim_corpus_release_tag <- function() {
    "sim-corpus-v1"
}

bfpwr_sim_corpus_release_base_url <- function(repo = bfpwr_sim_corpus_release_repo(),
                                              tag = bfpwr_sim_corpus_release_tag()) {
    paste0("https://github.com/", repo, "/releases/download/", tag)
}

bfpwr_sim_corpus_release_support_files <- function() {
    c("SHA256SUMS.txt", "corpus-manifest.json", "README-sim-corpus-v1.md")
}

bfpwr_sim_corpus_release_assets <- function(asset_set = "all") {
    assets <- c(
        metadata = "bfpwr-sim-corpus-v1-metadata.tar.zst",
        fixtures = "bfpwr-sim-corpus-v1-fixtures.tar.zst",
        fixture_validation = "bfpwr-sim-corpus-v1-fixture-validation.tar.zst",
        designs = "bfpwr-sim-corpus-v1-designs.tar.zst",
        bf_priors_binomial = "bfpwr-sim-corpus-v1-bf-priors-binomial.tar.zst",
        bf_priors_t_part1 = "bfpwr-sim-corpus-v1-bf-priors-t-part1.tar.zst",
        bf_priors_t_part2 = "bfpwr-sim-corpus-v1-bf-priors-t-part2.tar.zst",
        bf_priors_z_part1 = "bfpwr-sim-corpus-v1-bf-priors-z-part1.tar.zst",
        bf_priors_z_part2 = "bfpwr-sim-corpus-v1-bf-priors-z-part2.tar.zst"
    )
    sets <- list(
        all = unname(assets),
        metadata = unname(assets["metadata"]),
        fixture_tests = unname(assets[c("metadata", "fixtures")]),
        fixture_validation = unname(assets[c("metadata", "fixtures",
                                            "fixture_validation")]),
        chunk_validation = unname(assets[c("metadata", "designs",
                                          "bf_priors_binomial",
                                          "bf_priors_t_part1",
                                          "bf_priors_t_part2",
                                          "bf_priors_z_part1",
                                          "bf_priors_z_part2")]),
        bf_priors = unname(assets[c("metadata", "bf_priors_binomial",
                                    "bf_priors_t_part1",
                                    "bf_priors_t_part2",
                                    "bf_priors_z_part1",
                                    "bf_priors_z_part2")]),
        designs = unname(assets[c("metadata", "designs")])
    )
    normalized <- gsub("-", "_", asset_set)
    requested <- unique(trimws(strsplit(normalized, ",", fixed = TRUE)[[1]]))
    if (length(requested) == 0 || any(!nzchar(requested))) {
        requested <- "all"
    }
    unknown <- setdiff(requested, names(sets))
    if (length(unknown) > 0) {
        stop("unknown corpus asset set: ", paste(unknown, collapse = ", "),
             ". Expected one of: ", paste(names(sets), collapse = ", "))
    }
    unique(unlist(sets[requested], use.names = FALSE))
}

bfpwr_sim_default_corpus_root <- function(repo_root = ".") {
    file.path(repo_root, "simulations", "corpus", "v1")
}

bfpwr_sim_corpus_download_instructions <- function(
        corpus_root = bfpwr_sim_default_corpus_root(),
        repo_root = ".",
        asset_set = "fixture-tests",
        reason = NULL) {
    corpus_root <- normalizePath(corpus_root, winslash = "/", mustWork = FALSE)
    repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
    command <- paste(
        "Rscript simulations/scripts/download_corpus_release.R",
        "--corpus-root", shQuote(corpus_root),
        "--asset-set", asset_set
    )
    lines <- c()
    if (!is.null(reason) && nzchar(reason)) {
        lines <- c(lines, reason, "")
    }
    c(
        lines,
        paste0("The bfpwr simulation corpus is required at: ", corpus_root),
        "",
        "Download and extract the needed release assets with:",
        paste0("  ", command),
        "",
        paste0("Release: https://github.com/",
               bfpwr_sim_corpus_release_repo(), "/releases/tag/",
               bfpwr_sim_corpus_release_tag()),
        "",
        "Then run tests or scripts with:",
        paste0("  BFPWR_SIM_CORPUS=", corpus_root),
        paste0("  BFPWR_SIM_REPO=", repo_root)
    )
}

bfpwr_sim_corpus_instruction_text <- function(..., collapse = "\n") {
    paste(bfpwr_sim_corpus_download_instructions(...), collapse = collapse)
}

bfpwr_sim_require_corpus_root <- function(corpus_root,
                                          required_paths = character(),
                                          asset_set = "fixture-tests",
                                          purpose = NULL) {
    if (is.null(corpus_root) || !nzchar(corpus_root)) {
        reason <- if (is.null(purpose)) {
            "No corpus root was supplied."
        } else {
            paste0("No corpus root was supplied for ", purpose, ".")
        }
        stop(bfpwr_sim_corpus_instruction_text(reason = reason,
                                               asset_set = asset_set),
             call. = FALSE)
    }
    normalized <- normalizePath(corpus_root, winslash = "/", mustWork = FALSE)
    if (!dir.exists(normalized)) {
        stop(bfpwr_sim_corpus_instruction_text(
            corpus_root = normalized,
            reason = paste0("Corpus directory does not exist: ", normalized),
            asset_set = asset_set),
            call. = FALSE)
    }
    if (length(required_paths) > 0) {
        full_paths <- file.path(normalized, required_paths)
        missing <- full_paths[!file.exists(full_paths)]
        if (length(missing) > 0) {
            shown <- paste(utils::head(missing, 10), collapse = "\n  ")
            if (length(missing) > 10) {
                shown <- paste0(shown, "\n  ... and ",
                                length(missing) - 10, " more")
            }
            stop(bfpwr_sim_corpus_instruction_text(
                corpus_root = normalized,
                reason = paste0("Corpus is missing required files:\n  ",
                                shown),
                asset_set = asset_set),
                call. = FALSE)
        }
    }
    normalizePath(normalized, winslash = "/", mustWork = TRUE)
}

bfpwr_sim_stop_missing_corpus_files <- function(corpus_root,
                                                missing,
                                                asset_set = "chunk-validation",
                                                label = "corpus files") {
    shown <- paste(utils::head(missing, 10), collapse = "\n  ")
    if (length(missing) > 10) {
        shown <- paste0(shown, "\n  ... and ", length(missing) - 10, " more")
    }
    stop(bfpwr_sim_corpus_instruction_text(
        corpus_root = corpus_root,
        reason = paste0("Missing ", label, ":\n  ", shown),
        asset_set = asset_set),
        call. = FALSE)
}
