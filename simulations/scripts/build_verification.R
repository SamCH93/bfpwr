support_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    script <- if (length(hit) > 0) {
        normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                      winslash = "/", mustWork = TRUE)
    } else {
        normalizePath("simulations/scripts/build_verification.R",
                      winslash = "/", mustWork = FALSE)
    }
    file.path(dirname(script), "simulation_support.R")
}

source(support_path())
source_simulation_library()

elapsed_step <- function(label, expr) {
    cat("\n== ", label, " ==\n", sep = "")
    t0 <- proc.time()[["elapsed"]]
    value <- force(expr)
    cat("elapsed seconds:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
    invisible(value)
}

append_arg <- function(args, flag, value) {
    if (is.null(value) || !nzchar(as.character(value))) {
        return(args)
    }
    c(args, flag, as.character(value))
}

main <- function() {
    args <- parse_args()
    root <- repo_root()
    corpus_root <- normalizePath(
        arg_value(args, "corpus-root",
                  file.path(root, "simulations", "corpus", "v1")),
        winslash = "/", mustWork = FALSE)
    download_dir <- normalizePath(
        arg_value(args, "download-dir",
                  file.path(root, "simulations", "downloads",
                            bfpwr_sim_corpus_release_tag())),
        winslash = "/", mustWork = FALSE)

    asset_set <- arg_value(args, "asset-set", "fixture-validation")
    repo <- arg_value(args, "repo", bfpwr_sim_corpus_release_repo())
    tag <- arg_value(args, "tag", bfpwr_sim_corpus_release_tag())

    skip_fetch <- arg_flag(args, "skip-fetch")
    skip_unpack <- arg_flag(args, "skip-unpack")
    skip_recompute <- arg_flag(args, "skip-recompute")
    skip_validation <- arg_flag(args, "skip-validation")
    skip_report <- arg_flag(args, "skip-report")
    skip_search <- arg_flag(args, "skip-search")
    skip_package_manifest <- arg_flag(args, "skip-package-manifest")
    refresh_manifests <- arg_flag(args, "refresh-manifests")
    refresh_package_manifest <- refresh_manifests ||
        arg_flag(args, "refresh-package-manifest")
    refresh_search_manifest <- arg_flag(args, "refresh-search-manifest")
    require_search <- arg_flag(args, "require-search")
    search_manifest <- normalizePath(
        arg_value(args, "search-manifest",
                  bfpwr_sim_search_manifest_file(root)),
        winslash = "/", mustWork = FALSE)
    package_manifest <- normalizePath(
        arg_value(args, "package-manifest",
                  bfpwr_sim_package_verification_manifest_file(root)),
        winslash = "/", mustWork = FALSE)

    cat("Repo:", root, "\n")
    cat("Corpus:", corpus_root, "\n")
    cat("Release:", paste0(repo, "@", tag), "\n")
    cat("Prepared asset set:", asset_set, "\n")
    cat("Note: raw Monte Carlo corpus rebuilding is cluster-only; this script\n")
    cat("      consumes prepared release assets and recomputes package-side\n")
    cat("      verification outputs from the current checkout.\n")

    if (!skip_fetch) {
        download_args <- c("--corpus-root", corpus_root,
                           "--download-dir", download_dir,
                           "--asset-set", asset_set,
                           "--repo", repo,
                           "--tag", tag)
        if (arg_flag(args, "skip-download")) {
            download_args <- c(download_args, "--skip-download")
        }
        if (skip_unpack) {
            download_args <- c(download_args, "--no-extract")
        }
        elapsed_step(
            "fetch and unpack prepared corpus",
            run_rscript(file.path(root, "simulations", "scripts",
                                  "download_corpus_release.R"),
                        download_args)
        )
    }

    if (!skip_recompute) {
        if (!skip_search && refresh_search_manifest) {
            manifest_args <- c("--corpus-root", corpus_root,
                               "--output-file", search_manifest)
            manifest_args <- append_arg(
                manifest_args, "--achieved-per-function",
                arg_value(args, "search-achieved-per-function", NULL))
            manifest_args <- append_arg(
                manifest_args, "--diagnostic-per-function",
                arg_value(args, "search-diagnostic-per-function", NULL))
            manifest_args <- append_arg(
                manifest_args, "--min-achieved-per-function",
                arg_value(args, "search-min-achieved-per-function", NULL))
            manifest_args <- append_arg(
                manifest_args, "--required-mcse-margin",
                arg_value(args, "search-required-mcse-margin", NULL))
            manifest_args <- append_arg(manifest_args, "--families",
                                        arg_value(args, "search-families",
                                                  NULL))
            elapsed_step(
                "prepare sequential sample-size search cases",
                run_rscript(file.path(root, "simulations", "scripts",
                                      "build_sequential_search_validation_manifest.R"),
                            manifest_args)
                )
        }
        if (!skip_package_manifest) {
            package_manifest_args <- c("--corpus-root", corpus_root,
                                       "--output-file", package_manifest)
            if (arg_flag(args, "use-legacy-search-manifest")) {
                if (!refresh_package_manifest) {
                    stop("--use-legacy-search-manifest only applies when ",
                         "refreshing manifests. Add --refresh-manifests or ",
                         "--refresh-package-manifest.", call. = FALSE)
                }
                package_manifest_args <- c(package_manifest_args,
                                           "--use-legacy-search-manifest",
                                           "--search-manifest",
                                           search_manifest)
            }
            if (refresh_package_manifest) {
                elapsed_step(
                    "prepare package checks",
                    run_rscript(file.path(root, "simulations", "scripts",
                                          "build_package_verification_manifest.R"),
                                package_manifest_args)
                )
            } else if (!file.exists(package_manifest)) {
                stop("package check manifest is missing: ", package_manifest,
                     "\nRestore the committed file or rerun with ",
                     "--refresh-manifests after fetching the corpus.",
                     call. = FALSE)
            } else {
                cat("\nUsing package check manifest:", package_manifest, "\n")
            }
        }

        recompute_args <- c("--corpus-root", corpus_root)
        recompute_args <- append_arg(recompute_args, "--fixture-set",
                                     arg_value(args, "fixture-set", NULL))
        recompute_args <- append_arg(recompute_args, "--case-set",
                                     arg_value(args, "case-set", NULL))
        recompute_args <- append_arg(recompute_args, "--max-search-rows",
                                     arg_value(args, "max-search-rows", NULL))
        recompute_args <- append_arg(recompute_args, "--package-manifest",
                                     if (skip_package_manifest) NULL else
                                         package_manifest)
        recompute_args <- append_arg(recompute_args, "--search-manifest",
                                     if (skip_search ||
                                         !skip_package_manifest) NULL else
                                             search_manifest)
        if (arg_flag(args, "skip-fixtures")) {
            recompute_args <- c(recompute_args, "--skip-fixtures")
        }
        if (skip_search) {
            recompute_args <- c(recompute_args, "--skip-search")
        }
        if (require_search) {
            recompute_args <- c(recompute_args, "--require-search")
        }
        if (arg_flag(args, "evaluate-search-diagnostics")) {
            recompute_args <- c(recompute_args,
                                "--evaluate-search-diagnostics")
        }
        elapsed_step(
            "recompute package results",
            run_rscript(file.path(root, "simulations", "scripts",
                                  "recompute_package_references.R"),
                        recompute_args)
        )
    }

    if (!skip_validation) {
        validation_args <- c("--corpus-root", corpus_root,
                             "--search-scope",
                             arg_value(args, "search-scope", "smoke"))
        validation_args <- append_arg(validation_args, "--fixture-set",
                                      arg_value(args, "fixture-set", NULL))
        if (skip_search) {
            validation_args <- c(validation_args, "--skip-search")
        }
        if (require_search) {
            validation_args <- c(validation_args, "--require-search")
        }
        elapsed_step(
            "validate simulation outputs",
            run_rscript(file.path(root, "simulations", "scripts",
                                  "run_verification.R"),
                        validation_args)
        )
    }

    if (!skip_report) {
        report_args <- c("--corpus-root", corpus_root)
        elapsed_step(
            "render simulation verification PDF",
            run_rscript(file.path(root, "simulations", "scripts",
                                  "render_simulation_verification_background.R"),
                        report_args)
        )
    }

    cat("\nVerification build complete.\n")
    cat("Report:", file.path(root, "simulations",
                             "simulation-verification.pdf"), "\n")
}

main()
