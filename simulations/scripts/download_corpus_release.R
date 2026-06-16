support_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    script <- if (length(hit) > 0) {
        normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                      winslash = "/", mustWork = TRUE)
    } else {
        normalizePath("simulations/scripts/download_corpus_release.R",
                      winslash = "/", mustWork = FALSE)
    }
    file.path(dirname(script), "simulation_support.R")
}

source(support_path())
source_simulation_library()

args <- parse_args()
root <- repo_root()

corpus_root <- if (is.null(args[["corpus-root"]])) {
    bfpwr_sim_default_corpus_root(root)
} else {
    args[["corpus-root"]]
}
download_dir <- if (is.null(args[["download-dir"]])) {
    file.path(root, "simulations", "downloads",
              bfpwr_sim_corpus_release_tag())
} else {
    args[["download-dir"]]
}
asset_set <- if (is.null(args[["asset-set"]])) "all" else args[["asset-set"]]
repo <- if (is.null(args[["repo"]])) {
    bfpwr_sim_corpus_release_repo()
} else {
    args[["repo"]]
}
tag <- if (is.null(args[["tag"]])) {
    bfpwr_sim_corpus_release_tag()
} else {
    args[["tag"]]
}
skip_download <- isTRUE(args[["skip-download"]])
no_extract <- isTRUE(args[["no-extract"]])
verify_only <- isTRUE(args[["verify-only"]])

archive_assets <- bfpwr_sim_corpus_release_assets(asset_set)
support_assets <- bfpwr_sim_corpus_release_support_files()
assets <- unique(c(support_assets, archive_assets))

download_dir <- normalizePath(download_dir, winslash = "/", mustWork = FALSE)
corpus_root <- normalizePath(corpus_root, winslash = "/", mustWork = FALSE)
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- bfpwr_sim_corpus_release_base_url(repo = repo, tag = tag)
asset_url <- function(name) paste0(base_url, "/", name)
asset_path <- function(name) file.path(download_dir, name)

download_asset <- function(name, force = FALSE) {
    path <- asset_path(name)
    if (!force && file.exists(path)) {
        return(path)
    }
    if (skip_download || verify_only) {
        stop("missing downloaded release asset: ", path, call. = FALSE)
    }
    cat("downloading", name, "\n")
    status <- utils::download.file(asset_url(name), destfile = path,
                                   mode = "wb", quiet = FALSE)
    if (!identical(status, 0L)) {
        stop("download failed for release asset: ", name, call. = FALSE)
    }
    path
}

sha_file <- download_asset("SHA256SUMS.txt")
sha_lines <- readLines(sha_file, warn = FALSE)
sha_parts <- strsplit(sha_lines, "[[:space:]]+", perl = TRUE)
sha <- stats::setNames(
    vapply(sha_parts, `[`, character(1), 1L),
    vapply(sha_parts, function(x) x[[length(x)]], character(1)))

verify_asset <- function(name) {
    path <- asset_path(name)
    if (!file.exists(path)) {
        return(FALSE)
    }
    if (!name %in% names(sha)) {
        stop("SHA256SUMS.txt does not contain release asset: ", name,
             call. = FALSE)
    }
    actual <- unname(tools::sha256sum(path))
    identical(tolower(actual), tolower(sha[[name]]))
}

for (name in setdiff(assets, "SHA256SUMS.txt")) {
    if (file.exists(asset_path(name)) && verify_asset(name)) {
        cat("verified", name, "\n")
        next
    }
    if (file.exists(asset_path(name))) {
        if (skip_download || verify_only) {
            stop("checksum verification failed for release asset: ", name,
                 call. = FALSE)
        }
        cat("checksum mismatch; re-downloading", name, "\n")
    }
    download_asset(name, force = TRUE)
    if (!verify_asset(name)) {
        stop("checksum verification failed for release asset: ", name,
             call. = FALSE)
    }
    cat("verified", name, "\n")
}

if (!no_extract && !verify_only) {
    tar_bin <- Sys.which("tar")
    if (!nzchar(tar_bin)) {
        stop("could not find 'tar' on PATH; install tar or extract archives manually",
             call. = FALSE)
    }
    dir.create(corpus_root, recursive = TRUE, showWarnings = FALSE)
    for (name in archive_assets) {
        archive <- asset_path(name)
        cat("extracting", name, "to", corpus_root, "\n")
        status <- system2(tar_bin, c("-xf", shQuote(archive),
                                     "-C", shQuote(corpus_root)))
        if (!identical(status, 0L)) {
            stop("tar extraction failed for release asset: ", name,
                 call. = FALSE)
        }
    }
}

cat("\nCorpus release assets are ready.\n")
cat("Download directory:", download_dir, "\n")
cat("Corpus root:", corpus_root, "\n")
cat("\nSet these environment variables when running verification checks:\n")
cat("BFPWR_SIM_CORPUS=", corpus_root, "\n", sep = "")
cat("BFPWR_SIM_REPO=", root, "\n", sep = "")
