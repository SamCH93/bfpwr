support_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    script <- if (length(hit) > 0) {
        normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                      winslash = "/", mustWork = TRUE)
    } else {
        normalizePath("simulations/scripts/build_package_verification_manifest.R",
                      winslash = "/", mustWork = FALSE)
    }
    file.path(dirname(script), "simulation_support.R")
}

source(support_path())
source_simulation_library()

write_csv <- function(x, path) {
    ensure_dir(dirname(path))
    utils::write.csv(x, path, row.names = FALSE, na = "")
    invisible(path)
}

main <- function() {
    args <- parse_args()
    root <- repo_root()
    corpus_root <- normalizePath(
        arg_value(args, "corpus-root",
                  file.path(root, "simulations", "corpus", "v1")),
        winslash = "/", mustWork = TRUE)
    output_file <- normalizePath(
        arg_value(args, "output-file",
                  bfpwr_sim_package_verification_manifest_file(root)),
        winslash = "/", mustWork = FALSE)
    coverage_file <- normalizePath(
        arg_value(args, "coverage-file",
                  bfpwr_sim_package_verification_coverage_file(root)),
        winslash = "/", mustWork = FALSE)
    search_manifest <- NULL
    if (arg_flag(args, "use-legacy-search-manifest")) {
        search_manifest_file <- normalizePath(
            arg_value(args, "search-manifest",
                      bfpwr_sim_search_manifest_file(root)),
            winslash = "/", mustWork = FALSE)
        if (!file.exists(search_manifest_file)) {
            stop("legacy sequential search manifest is missing: ",
                 search_manifest_file, call. = FALSE)
        }
        search_manifest <- utils::read.csv(search_manifest_file,
                                           stringsAsFactors = FALSE)
    }

    manifest <- bfpwr_sim_build_package_verification_manifest(
        corpus_root = corpus_root,
        search_manifest = search_manifest)
    coverage <- bfpwr_sim_package_verification_coverage(manifest)

    write_csv(manifest, output_file)
    write_csv(coverage, coverage_file)

    cat("wrote package verification manifest:", output_file, "\n")
    cat("rows:", nrow(manifest), "\n")
    if (nrow(coverage) > 0) {
        cat("coverage groups:", nrow(coverage), "\n")
    }
}

main()
