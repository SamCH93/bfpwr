## Download immutable versioned assets; source() only defines the entry point.
bfpwr_sim_download_moment_extension <- function(args, root) {
    source(file.path(root, "simulations/scripts/simulation_support.R"), local = .GlobalEnv)
    source(file.path(root, "simulations/R/moment-identity.R"), local = .GlobalEnv)
    args <- parse_args(args)
    version <- arg_value(args, "version", "moment-v2")
    stopifnot(version %in% c("moment-v1", "moment-v2"))
    output <- arg_value(args, "output", file.path(root, "simulations/corpus", version))
    downloads <- arg_value(args, "downloads", file.path(root, "simulations/downloads",
                                                       paste0("sim-corpus-v1-", version)))
    ensure_dir(downloads)
    ensure_dir(output)
    base <- arg_value(args, "release-base", "https://github.com/FBartos/bfpwr/releases/download/sim-corpus-v1/")
    if (!endsWith(base, "/")) base <- paste0(base, "/")
    fetch <- function(name, sha256 = NULL, url = paste0(base, name)) {
        path <- file.path(downloads, name)
        if (is.null(sha256) || !file.exists(path) ||
            !identical(unname(tools::sha256sum(path)), sha256)) {
            utils::download.file(url, path, mode = "wb", quiet = TRUE)
        }
        if (!is.null(sha256)) stopifnot(identical(unname(tools::sha256sum(path)), sha256))
        path
    }
    checksum_name <- paste0("SHA256SUMS-", version, ".txt")
    parts <- strsplit(trimws(readLines(fetch(checksum_name))), "[[:space:]]+")
    sha <- setNames(vapply(parts, `[`, "", 1L), vapply(parts, function(x) tail(x, 1L), ""))
    if (version == "moment-v2") {
        manifest_name <- paste0(version, "-manifest.csv")
        stopifnot(manifest_name %in% names(sha))
        assets <- read.csv(fetch(manifest_name, sha[[manifest_name]]), stringsAsFactors = FALSE)
        stopifnot(all(assets$version == version), !anyDuplicated(assets$name),
                  all(assets$historical | assets$sha256 == unname(sha[assets$name])))
        if (!is.null(args[["release-base"]])) {
            assets$url[!assets$historical] <- paste0(base, assets$name[!assets$historical])
        }
        assets <- assets[order(match(assets$asset_set, c("bf-chunks", "results", "sources", "metadata"))), ]
    } else {
        names <- names(sha)[grepl("[.]tar[.]zst$", names(sha))]
        assets <- data.frame(name = names, sha256 = unname(sha[names]),
            url = paste0(base, names), asset_set = ifelse(grepl("bf-chunks", names), "bf-chunks", "other"))
    }
    if (!arg_flag(args, "with-bf-chunks")) assets <- assets[assets$asset_set != "bf-chunks", ]
    extract <- function(path, target) {
        bfpwr_sim_moment_archive_members(path)
        stopifnot(system2("tar", c("-xf", shQuote(path), "-C", shQuote(target))) == 0L)
    }
    for (i in seq_len(nrow(assets))) {
        extract(fetch(assets$name[i], assets$sha256[i], assets$url[i]), output)
    }
    if (version == "moment-v1") {
        for (phase in c("bf-evaluation", "reference")) {
            archive <- file.path(output, paste0(phase, "-source.tar"))
            if (file.exists(archive)) {
                target <- file.path(output, "sources", paste0(phase, "-historical"))
                ensure_dir(target)
                extract(archive, target)
            }
        }
    }
    cat("Moment extension", version, "ready at", output, "\n")
    invisible(output)
}

if (sys.nframe() == 0L) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
    bfpwr_sim_download_moment_extension(commandArgs(TRUE), root)
}
