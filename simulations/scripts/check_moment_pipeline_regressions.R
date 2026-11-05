## Small end-to-end corpus engineering checks, isolated from all published data.
if (!file.exists("simulations/R/moment-reference-check.R")) {
    stop("Missing required simulations/R/moment-reference-check.R", call. = FALSE)
}
source("simulations/scripts/simulation_support.R")
repo_root <- function() normalizePath(getwd(), winslash = "/")
source_simulation_library()
source_package_checkout()
source("simulations/R/corpus-io.R")
source("simulations/R/corpus-validation.R")
source("simulations/R/moment-extension.R")
source("simulations/R/moment-identity.R")
scratch <- tempfile("engineering-", tmpdir = file.path(repo_root(), "simulations/logs/runner-fixes"))
if (dir.exists(scratch)) stop("use a fresh scratch destination; existing regression evidence is retained")
dir.create(scratch, recursive = TRUE)
scratch <- normalizePath(scratch, winslash = "/")
writeLines(scratch, "simulations/logs/runner-fixes/latest-engineering.txt")
root <- file.path(scratch, "source")
moment_root <- file.path(scratch, "moment inputs")
stopifnot(dir.create(moment_root))
missing_bundle <- tryCatch(resolve_moment_root(root, moment_root), error = conditionMessage)
stopifnot(grepl("requires a saved moment verification.rds bundle", missing_bundle, fixed = TRUE))
saveRDS(list(), file.path(moment_root, "verification.rds"))
moment_root <- resolve_moment_root(root, moment_root)
child <- file.path(scratch, "report argument probe.R")
child_log <- file.path(scratch, "report argument probe.rds")
writeLines("saveRDS(commandArgs(TRUE)[-1L], commandArgs(TRUE)[[1L]])", child)
run_rscript(child, c(child_log, "--moment-root", moment_root))
stopifnot(identical(readRDS(child_log), c("--moment-root", moment_root)))

files <- unique(c(bfpwr_sim_moment_source_files(repo_root(), "reference"),
    "simulations/scripts/package_moment_extension.R"))
for (file in files) {
    target <- file.path(root, file)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(repo_root(), file), target))
}
## Check collation only in this test process and restore it on return.
check_collation <- function() {
    previous <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", previous))
    stopifnot(nzchar(Sys.setlocale("LC_COLLATE", "C")))
    files <- bfpwr_sim_moment_source_files(root, "reference")
    expected <- bfpwr_sim_moment_source_manifest(root, files)
    for (locale in c("English_United States.utf8", "en_US.UTF-8", "C.UTF-8")) {
        if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) next
        stopifnot(identical(expected, bfpwr_sim_moment_source_manifest(root,
            bfpwr_sim_moment_source_files(root, "reference"))))
        cat("Source manifest ordering identical under C and", locale, "\n")
        return(locale)
    }
    cat("SKIP alternate collation: no English/UTF8 locale available.\n")
    "C"
}
alternate_collation <- check_collation()
case <- bfpwr_sim_moment_extension_cases()[[1]]
design <- bfpwr_sim_find_design_case(case$design_case_id)
design$look_grid <- seq(10, 50, 10)
design$nsim <- 20L
design$chunks <- bfpwr_sim_make_chunks(20L, 10L, design$master_seed)
## These private registries keep every production phase unchanged.
saveRDS(case, file.path(scratch, "case.rds"))
saveRDS(design, file.path(scratch, "design.rds"))
cat("\nbfpwr_sim_moment_extension_cases <- function() list(readRDS(",
    dQuote(file.path(scratch, "case.rds"), FALSE), "))\n", sep = "",
    file = file.path(root, "simulations/R/moment-extension.R"), append = TRUE)
cat("\nbfpwr_sim_find_design_case <- function(...) readRDS(",
    dQuote(file.path(scratch, "design.rds"), FALSE), ")\n", sep = "",
    file = file.path(root, "simulations/registry/design-cases.R"), append = TRUE)
input <- file.path(scratch, "input")
output <- file.path(scratch, "fresh-v2")
set.seed(812)
for (j in 1:2) {
    ids <- seq(design$chunks$replicate_start[j], design$chunks$replicate_end[j])
    rows <- do.call(rbind, lapply(ids, function(id) data.frame(
        design_case_id = design$design_case_id, replicate_id = id, chunk_id = j,
        look = seq_along(design$look_grid), n = design$look_grid,
        true_effect = 0, estimate = cumsum(rnorm(5, 0, sqrt(10)))/design$look_grid,
        se = design$generation$usd/sqrt(design$look_grid))))
    bfpwr_sim_write_rds(rows, bfpwr_sim_design_chunk_files(input, design)[j])
}
bfpwr_sim_write_rds(design, file.path(bfpwr_sim_design_dir(input, design), "design.rds"))
invisible(bfpwr_sim_write_design_chunk_index(input, design))
run <- function(phase, out = output, extra = character(), expect = 0L, script = "run_moment_extension.R") {
    log <- file.path(scratch, paste0(length(list.files(scratch, "[.]log$")) + 1L, "-", phase, ".log"))
    status <- system2(file.path(R.home("bin"), "Rscript"), c(
        shQuote(file.path(root, "simulations/scripts", script)),
        "--phase", phase, "--input-corpus", shQuote(input), "--output", shQuote(out), extra),
        stdout = log, stderr = log)
    cat(phase, "exit", status, "expected", expect, "log", basename(log), "\n")
    if (status != expect) { cat(readLines(log), sep = "\n"); stop("unexpected phase status") }
    invisible(log)
}
run("plan")
## Planning creates both phase identities and self-contained source snapshots.
stopifnot(all(file.exists(file.path(output, c("bf-evaluation-identity.rds",
    "bf-evaluation-provenance.rds", "sources/bf/source-manifest.rds",
    "sources/reference/source-manifest.rds")))))
package_manifest <- "simulations/registry/manifests/package-verification-cases.csv"
snapshot_root <- file.path(output, "sources/reference")
snapshot_manifest <- readRDS(file.path(snapshot_root, "source-manifest.rds"))
snapshot_csv <- bfpwr_sim_package_verification_manifest_file(snapshot_root)
original_csv <- bfpwr_sim_package_verification_manifest_file(root)
stopifnot(file.exists(snapshot_csv),
    identical(snapshot_manifest$sha256[match(package_manifest, snapshot_manifest$file)],
        unname(tools::sha256sum(original_csv))),
    identical(bfpwr_sim_read_package_verification_manifest(snapshot_csv, required = TRUE),
        bfpwr_sim_read_package_verification_manifest(original_csv, required = TRUE)),
    !package_manifest %in% bfpwr_sim_moment_dependency_files(root, "reference"),
    !package_manifest %in% bfpwr_sim_moment_dependency_files(root, "bf"))
## Each report caller must reject a missing bundle before reaching corpus work.
missing_moment_root <- file.path(scratch, "missing moment inputs")
missing_corpus <- file.path(scratch, "missing corpus")
prerequisite <- paste0("requires a saved moment verification.rds bundle at ",
    file.path(missing_moment_root, "verification.rds"))
for (script in c("build_verification.R", "refresh_simulation_grid.R", "verify_one_sided_z.R")) {
    log <- file.path(scratch, paste0("report-prerequisite-", script, ".txt"))
    status <- system2(file.path(R.home("bin"), "Rscript"), c(
        shQuote(file.path(snapshot_root, "simulations/scripts", script)),
        "--moment-root", shQuote(missing_moment_root),
        "--corpus-root", shQuote(missing_corpus), "--skip-fetch"), stdout = log, stderr = log)
    messages <- readLines(log)
    if (status != 1L || !any(grepl(prerequisite, messages, fixed = TRUE))) {
        cat(messages, sep = "\n")
        stop("Unexpected report prerequisite result: ", script)
    }
    cat(script, "rejected the missing report prerequisite before corpus work.\n")
}
run("chunks", extra = c("--workers", "2", "--worker", "1"))
chunk_paths <- file.path(output, "chunk-status", sprintf("%04d.rds", 1:2))
chunk_times <- file.info(chunk_paths)$mtime
run("chunks")
stopifnot(identical(chunk_times, file.info(chunk_paths)$mtime))
## A fresh-BF runtime mismatch is rejected before any completed or missing chunk.
bf_identity_path <- file.path(output, "bf-dependency-identity.rds")
bf_identity_raw <- readBin(bf_identity_path, "raw", file.info(bf_identity_path)$size)
changed_runtime <- readRDS(bf_identity_path)
stopifnot(identical(names(changed_runtime$values$runtime), c("R", "platform", "stats")))
changed_runtime$values$runtime$platform <- "different-platform"
saveRDS(changed_runtime, bf_identity_path)
run("chunks", expect = 1L)
writeBin(bf_identity_raw, bf_identity_path)
status_raw <- readBin(chunk_paths[1], "raw", file.info(chunk_paths[1])$size)
changed_runtime <- readRDS(chunk_paths[1])
changed_runtime$identity$values$runtime$stats <- "different-stats"
saveRDS(changed_runtime, chunk_paths[1])
run("chunks", expect = 1L)
writeBin(status_raw, chunk_paths[1])
chunk_times <- file.info(chunk_paths)$mtime
cat("Changed fresh BF dependency runtime and chunk-status runtime rejected.\n")
## A partial resume fills only the missing status, retaining an existing BF result.
unlink(chunk_paths[2])
run("chunks")
stopifnot(identical(chunk_times[1], file.info(chunk_paths[1])$mtime))
run("summaries")
run("finalize")
cat("Fresh planning, BF evaluation, partial/complete resume, references and finalization passed.\n")
## Presentation tools are retained in the snapshot without invalidating calculations.
presentation <- c("simulations/moment-extension-body.Rnw", "simulations/MOMENT-EXTENSION.md",
    "simulations/scripts/render_moment_extension.R", "simulations/scripts/prepare_moment_extension_report.R",
    "package/README.md", "package/NEWS.md",
    paste0("package/man/", list.files(file.path(root, "package/man"), "[.]Rd$")[1L]),
    package_manifest)
presentation_raw <- lapply(file.path(root, presentation), function(path) readBin(path, "raw", file.info(path)$size))
for (file in presentation) cat("\n# presentation-only edit\n", file = file.path(root, file), append = TRUE)
run("plan")
for (i in seq_along(presentation)) writeBin(presentation_raw[[i]], file.path(root, presentation[i]))
## Calculation and validation code both participate in reference compatibility.
for (file in c("simulations/R/moment-extension.R", "simulations/scripts/validate_moment_extension.R")) {
    path <- file.path(root, file)
    raw <- readBin(path, "raw", file.info(path)$size)
    cat("\n# calculation dependency edit\n", file = path, append = TRUE)
    run("plan", expect = 1L)
    writeBin(raw, path)
}
## Identity checks: changed registry code, defaults and input index bytes must fail.
registry <- file.path(root, "simulations/registry/fixtures/thresholds.R")
original <- readBin(registry, "raw", file.info(registry)$size)
writeLines(sub("c(3, 10, 30)", "c(3, 10, 31)", readLines(registry), fixed = TRUE), registry)
run("plan", expect = 1L)
writeBin(original, registry)
## Change the effective numerical option in a real fresh worker process.
profile <- file.path(scratch, "changed-options.R")
writeLines("options(bfpwr.ngrid = 10001)", profile)
old_profile <- Sys.getenv("R_PROFILE_USER", unset = NA_character_)
Sys.setenv(R_PROFILE_USER = profile)
run("plan", expect = 1L)
if (is.na(old_profile)) Sys.unsetenv("R_PROFILE_USER") else Sys.setenv(R_PROFILE_USER = old_profile)
index_path <- file.path(bfpwr_sim_design_dir(input, design), "chunk-index.rds")
index_raw <- readBin(index_path, "raw", file.info(index_path)$size)
index <- readRDS(index_path)
index$sha256[1] <- paste(rep("0", 64), collapse = "")
saveRDS(index, index_path)
run("plan", expect = 1L)
writeBin(index_raw, index_path)
## Reference-only source changes do not invalidate mathematically unchanged BF work.
reference_file <- file.path(root, "package/R/pbf01seq.R")
reference_raw <- readBin(reference_file, "raw", file.info(reference_file)$size)
cat("\n# changed reference source\n", file = reference_file, append = TRUE)
run("chunks")
run("plan", expect = 1L)
writeBin(reference_raw, reference_file)
## Design metadata must match the complete frozen specification, not just its ID.
design_path <- file.path(bfpwr_sim_design_dir(input, design), "design.rds")
design_raw <- readBin(design_path, "raw", file.info(design_path)$size)
different_design <- design
different_design$design_prior$mean <- different_design$design_prior$mean + .1
saveRDS(different_design, design_path)
run("plan", expect = 1L)
writeBin(design_raw, design_path)
## Raw input damage with an unchanged index is rejected even on a completed chunk resume.
input_path <- bfpwr_sim_design_chunk_files(input, design)[1]
input_raw <- readBin(input_path, "raw", file.info(input_path)$size)
writeBin(c(input_raw, as.raw(0)), input_path)
run("chunks", expect = 1L)
writeBin(input_raw, input_path)
## CRLF/LF and trailing empty lines do not stale computational identities.
writeBin(charToRaw(paste0(paste(readLines(registry), collapse = "\r\n"), "\r\n\r\n")), registry)
run("plan")
writeBin(original, registry)
## Source snapshots still use raw hashes, and corruption is never normalized away.
snapshot <- file.path(output, "sources/reference", "package/R/nmbf01.R")
snapshot_raw <- readBin(snapshot, "raw", file.info(snapshot)$size)
writeBin(c(snapshot_raw, as.raw(0)), snapshot)
run("plan", expect = 1L)
writeBin(snapshot_raw, snapshot)
for (path in c(file.path(snapshot_root, "simulations/moment-extension-body.Rnw"), snapshot_csv)) {
    snapshot_raw <- readBin(path, "raw", file.info(path)$size)
    cat("\n% damaged presentation snapshot\n", file = path, append = TRUE)
    run("plan", expect = 1L)
    writeBin(snapshot_raw, path)
}
cat("Presentation edits accepted; calculation/registry/validation/default/index/input changes and all raw snapshot damage rejected; harmless line endings accepted.\n")




run("assemble")
bundle <- readRDS(file.path(output, "verification.rds"))
stopifnot(bundle$timing_schema_version == 2L, length(bundle$designs) == 1L,
          nrow(bundle$numerical_diagnostics$moment_rows) == nrow(bundle$moments))
saveRDS(bundle, file.path(scratch, "small-bundle-v2.rds"))
## A synthetic v1 fixture retains the old status identity and frozen provenance.
old <- file.path(scratch, "synthetic-v1")
dir.create(old)
stopifnot(file.copy(output, old, recursive = TRUE))
old <- file.path(old, basename(output))
historical <- file.path(scratch, "historical-source")
for (file in bfpwr_sim_moment_source_files(root, "bf")) {
    target <- file.path(historical, file)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(root, file), target))
}
files_bf <- bfpwr_sim_moment_source_files(root, "bf")
provenance <- readRDS(file.path(old, "bf-evaluation-provenance.rds"))
provenance$source_checksums <- setNames(unname(tools::md5sum(file.path(historical, files_bf))), files_bf)
historical_identity <- list(md5 = unname(provenance$source_checksums), defaults = bfpwrOptions())
saveRDS(provenance, file.path(old, "bf-evaluation-provenance.rds"))
saveRDS(historical_identity, file.path(old, "bf-evaluation-identity.rds"))
for (path in file.path(old, "chunk-status", sprintf("%04d.rds", 1:2))) {
    status <- readRDS(path)
    status$identity <- historical_identity
    status$input <- if (endsWith(path, "0001.rds")) {
        "/historical/corpus/designs/original/chunk-0001.rds"
    } else "designs/original/chunk-0002.rds"
    saveRDS(status, path)
}
empirical_columns <- c("max_observed_n", "tolerance_EN", "tolerance_VarN", "error_EN",
    "error_VarN", "relative_error_EN", "relative_error_VarN", "z_EN", "z_VarN")
bundle$moments <- bundle$moments[setdiff(names(bundle$moments), empirical_columns)]
saveRDS(bundle, file.path(old, "verification.rds"))
imported <- file.path(scratch, "imported-v2")
run("plan", imported)
run("import", imported, c("--reuse-from", shQuote(old), "--bf-source-root", shQuote(historical)))
stopifnot(identical(readRDS(file.path(imported, "bf-evaluation-identity.rds")), historical_identity))
for (i in 1:2) {
    file <- file.path("chunk-status", sprintf("%04d.rds", i))
    original_status <- readRDS(file.path(old, file))
    imported_status <- readRDS(file.path(imported, file))
    stopifnot(identical(imported_status$historical_input, original_status$input),
        imported_status$input_path_base == "provenance$source_corpus",
        imported_status$historical_input_path_base == if (i == 1L) {
            "absolute historical path"
        } else "bf-evaluation-provenance$source_corpus")
}
cat("Verbatim historical absolute/relative inputs retained with accurate path bases.\n")
run("summaries", imported)
run("finalize", imported)
run("assemble", imported)
imported_bundle <- readRDS(file.path(imported, "verification.rds"))
cat("Compatible historical BF reuse, fresh references, diagnostics and exact empirical equality passed.\n")
## Unknown BF changes are rejected, even when the new output has its own valid identity.
bf_file <- file.path(root, "package/R/nmbf01.R")
bf_original <- readBin(bf_file, "raw", file.info(bf_file)$size)
cat("\n# unknown BF change\n", file = bf_file, append = TRUE)
incompatible <- file.path(scratch, "incompatible-v2")
run("plan", incompatible)
run("import", incompatible, c("--reuse-from", shQuote(old), "--bf-source-root", shQuote(historical)), expect = 1L)
writeBin(bf_original, bf_file)
## A BF output hash change is rejected by import, preserving v1 as read-only evidence.
old_bf <- bfpwr_sim_bf_chunk_files(old, case, design)[1]
bf_raw <- readBin(old_bf, "raw", file.info(old_bf)$size)
writeBin(c(bf_raw, as.raw(0)), old_bf)
damaged <- file.path(scratch, "damaged-v2")
run("plan", damaged)
run("import", damaged, c("--reuse-from", shQuote(old), "--bf-source-root", shQuote(historical)), expect = 1L)
writeBin(bf_raw, old_bf)
## Package/download the fresh output without preparation/render files or manual tar steps.
search_checks <- bfpwr_sim_moment_search_failures(imported_bundle$searches)
checks <- data.frame(check = names(search_checks), failures = unname(search_checks))
stopifnot(!any(checks$failures))
write.csv(checks, file.path(imported, "validation-checks.csv"), row.names = FALSE)
numerical <- imported_bundle$numerical_diagnostics$moment_rows
write.csv(data.frame(check = "Reference convergence", failures = sum(!numerical$reference_converged)),
    file.path(imported, "numerical-validation-checks.csv"), row.names = FALSE)
release <- file.path(scratch, "release")
## Retried-task evidence remains local and never enters an archive or its identity.
evidence <- c("chunk-status/0001.rds.error", "case-results/retry.tmp",
    "case-results/retry.lock/evidence", "chunk-status/stale.lock")
for (file in evidence) {
    dir.create(dirname(file.path(imported, file)), recursive = TRUE, showWarnings = FALSE)
    writeLines("retained local evidence", file.path(imported, file))
}
## Repeat packaging under another collation without changing the user's locale.
package_collated <- function(locale) {
    profile <- file.path(scratch, "packaging-locale.R")
    writeLines(paste0("Sys.setlocale('LC_COLLATE', '", locale, "')"), profile)
    previous <- Sys.getenv("R_PROFILE_USER", unset = NA_character_)
    on.exit(if (is.na(previous)) Sys.unsetenv("R_PROFILE_USER") else Sys.setenv(R_PROFILE_USER = previous))
    Sys.setenv(R_PROFILE_USER = profile)
    run("package", imported, c("--release", shQuote(release)), script = "package_moment_extension.R")
}
numerical_path <- file.path(imported, "numerical-diagnostics.rds")
numerical_sha256 <- unname(tools::sha256sum(numerical_path))
package_collated("C")
manifest <- read.csv(file.path(release, "moment-v2-manifest.csv"))
stopifnot(!any(manifest$historical), sum(manifest$asset_set == "bf-chunks") == 1L,
          any(grepl("moment-v2-bf-chunks", list.files(release), fixed = TRUE)))
package_collated(alternate_collation)
payload_identity <- readRDS(file.path(release, "moment-v2-packaging-identity.rds"))
stopifnot(!any(evidence %in% payload_identity$file),
    all(file.exists(file.path(imported, evidence))),
    identical(numerical_sha256, unname(tools::sha256sum(numerical_path))))
for (name in manifest$name[!manifest$historical]) {
    members <- bfpwr_sim_moment_archive_members(file.path(release, name))
    stopifnot(all(members %in% payload_identity$file), !any(evidence %in% members))
}
## Existing directories do not prove that all expected per-case indexes exist.
index_path <- file.path(imported, "index-status/0001.rds")
index_raw <- readBin(index_path, "raw", file.info(index_path)$size)
unlink(index_path)
run("package", imported, c("--release", shQuote(file.path(scratch, "missing-index-release"))),
    script = "package_moment_extension.R", expect = 1L)
writeBin(index_raw, index_path)
cat("C/alternate-locale packaging identity, complete archive coverage, retained errors, unchanged diagnostics and missing-index rejection passed.\n")
## A present but mismatched index must identify the case and failed condition.
index_status <- readRDS(index_path)
index_status$sha256[1L] <- paste0(index_status$sha256[1L], "-corrupted")
saveRDS(index_status, index_path)
invalid_release <- file.path(scratch, "invalid-index-release")
index_log <- run("invalid-index-status", imported, c("--release", shQuote(invalid_release)),
    script = "package_moment_extension.R", expect = 1L)
writeBin(index_raw, index_path)
stopifnot(any(grepl(paste0("Invalid BF metadata indexes for case ", case$bf_prior_id,
    ": index_status_matches_chunk_index"), readLines(index_log), fixed = TRUE)),
    !length(list.files(invalid_release, all.files = TRUE, no.. = TRUE)))
cat("Mismatched index status identified its case and condition before release publication.\n")
downloaded <- file.path(scratch, "downloaded-v2")
run("download", downloaded, c("--release-base", shQuote(paste0("file:///", release, "/")),
    "--downloads", shQuote(file.path(scratch, "downloads"))), script = "download_moment_extension.R")
stopifnot(identical(unname(tools::sha256sum(file.path(imported, "verification.rds"))),
                    unname(tools::sha256sum(file.path(downloaded, "verification.rds")))),
          file.exists(file.path(downloaded, "sources/reference/simulations/scripts/run_moment_extension.R")),
          !dir.exists(file.path(bfpwr_sim_bf_prior_dir(downloaded, case), "bf_chunks")))
cat("Generic tiny v1 import packages its actual BF chunks; repeat packaging and fresh result-only download passed.\n")
generic_download <- file.path(scratch, "generic-download-with-bf")
run("download", generic_download, c("--with-bf-chunks",
    "--release-base", shQuote(paste0("file:///", release, "/")),
    "--downloads", shQuote(file.path(scratch, "generic-downloads"))), script = "download_moment_extension.R")
stopifnot(identical(unname(tools::sha256sum(bfpwr_sim_bf_chunk_files(imported, case, design))),
                    unname(tools::sha256sum(bfpwr_sim_bf_chunk_files(generic_download, case, design)))))
## A correctly checksummed archive containing a hard link is still rejected.
link_release <- file.path(scratch, "link-release")
link_input <- file.path(scratch, "link-input")
dir.create(link_release)
dir.create(link_input)
writeLines("regular file", file.path(link_input, "target"))
stopifnot(file.link(file.path(link_input, "target"), file.path(link_input, "link")))
link_name <- "rejected-link.tar.zst"
link_archive <- file.path(link_release, link_name)
stopifnot(system2("tar", c("--zstd", "-cf", shQuote(link_archive), "-C", shQuote(link_input),
    "target", "link")) == 0L)
link_types <- system2("tar", c("-tvf", shQuote(link_archive)), stdout = TRUE)
stopifnot(any(substr(link_types, 1L, 1L) %in% c("h", "l")))
link_manifest <- data.frame(name = link_name, asset_set = "results",
    sha256 = unname(tools::sha256sum(link_archive)), bytes = file.info(link_archive)$size,
    url = paste0("file:///", link_archive), historical = FALSE, version = "moment-v2")
link_manifest_path <- file.path(link_release, "moment-v2-manifest.csv")
write.csv(link_manifest, link_manifest_path, row.names = FALSE)
writeLines(paste(unname(tools::sha256sum(c(link_archive, link_manifest_path))),
    c(link_name, "moment-v2-manifest.csv")), file.path(link_release, "SHA256SUMS-moment-v2.txt"))
link_output <- file.path(scratch, "rejected-link-output")
run("download-link", link_output, c("--release-base", shQuote(paste0("file:///", link_release, "/")),
    "--downloads", shQuote(file.path(scratch, "link-downloads"))),
    script = "download_moment_extension.R", expect = 1L)
stopifnot(!length(list.files(link_output)))
cat("Downloader rejected a checksummed link archive before extraction.\n")
## The fresh-evaluation branch packages its own BF data, and optional BF downloads work.
for (file in c("validation-checks.csv", "numerical-validation-checks.csv")) {
    stopifnot(file.copy(file.path(imported, file), file.path(output, file)))
}
fresh_release <- file.path(scratch, "fresh-release")
run("package", output, c("--release", shQuote(fresh_release)), script = "package_moment_extension.R")
fresh_download <- file.path(scratch, "fresh-download-with-bf")
run("download", fresh_download, c("--with-bf-chunks",
    "--release-base", shQuote(paste0("file:///", fresh_release, "/")),
    "--downloads", shQuote(file.path(scratch, "fresh-downloads"))), script = "download_moment_extension.R")
stopifnot(identical(unname(tools::sha256sum(bfpwr_sim_bf_chunk_files(output, case, design))),
                    unname(tools::sha256sum(bfpwr_sim_bf_chunk_files(fresh_download, case, design)))))
## Fresh releases must reject damaged data even without an existing release identity.
reject_fresh_payload <- function(label, path, damage) {
    original <- readBin(path, "raw", file.info(path)$size)
    on.exit(writeBin(original, path))
    damage(path)
    destination <- file.path(scratch, paste0(label, "-fresh-release"))
    stopifnot(dir.create(destination))
    run(label, output, c("--release", shQuote(destination)),
        script = "package_moment_extension.R", expect = 1L)
    stopifnot(!length(list.files(destination, all.files = TRUE, no.. = TRUE)))
}
chunk <- bfpwr_sim_bf_chunk_files(output, case, design)[2L]
reject_fresh_payload("missing-bf", chunk, unlink)
reject_fresh_payload("corrupt-bf", chunk, function(path) writeLines("damaged BF", path))
reject_fresh_payload("corrupt-bf-index",
    file.path(bfpwr_sim_bf_prior_dir(output, case), "bf-chunk-index.rds"), function(path) {
        index <- readRDS(path)
        index$sha256[2L] <- strrep("0", 64)
        saveRDS(index, path)
    })
reject_fresh_payload("incomplete-plan", file.path(output, "chunk-plan.csv"), function(path) {
    write.csv(read.csv(path)[1L, , drop = FALSE], path, row.names = FALSE)
})
## A hard-linked file has valid bytes but would be rejected by the downloader.
reject_linked_payload <- function() {
    target <- file.path(output, "case-results", "0001.rds")
    link <- file.path(output, "case-results", "0001-link.rds")
    stopifnot(file.link(target, link))
    on.exit(unlink(link))
    destination <- file.path(scratch, "linked-payload-release")
    run("linked-payload", output, c("--release", shQuote(destination)),
        script = "package_moment_extension.R", expect = 1L)
    stopifnot(!length(list.files(destination, "[.](tar[.]zst|sha256|csv|tmp)$")))
}
reject_linked_payload()
## Repackaging never blesses a damaged archive or changed result payload.
archive_path <- file.path(release, manifest$name[!manifest$historical][1L])
archive_raw <- readBin(archive_path, "raw", file.info(archive_path)$size)
writeBin(c(archive_raw, as.raw(0)), archive_path)
run("package", imported, c("--release", shQuote(release)), script = "package_moment_extension.R", expect = 1L)
writeBin(archive_raw, archive_path)
payload <- file.path(imported, "timings.csv")
payload_raw <- readBin(payload, "raw", file.info(payload)$size)
cat("\n", file = payload, append = TRUE)
run("package", imported, c("--release", shQuote(release)), script = "package_moment_extension.R", expect = 1L)
writeBin(payload_raw, payload)
cat("Fresh BF packaging/download and archive/payload corruption rejection passed.\n")

## A private metadata/catalog fixture exercises the same optional official proof.
## The production catalog has no command-line override and remains immutable.
archive <- file.path(scratch, "private-metadata.tar.zst")
prior <- file.path("bf-priors", case$test_family, case$bf_prior_id)
stopifnot(system2("tar", c("--zstd", "-cf", shQuote(archive), "-C", shQuote(old),
    shQuote(file.path(prior, "bf-prior.rds")), shQuote(file.path(prior, "bf-chunk-index.rds")))) == 0L)
catalog <- bfpwr_sim_moment_v1_catalog()
catalog$metadata_sha256 <- unname(tools::sha256sum(archive))
catalog$assets <- manifest[manifest$asset_set == "bf-chunks", setdiff(names(manifest), "version")]
rownames(catalog$assets) <- NULL
catalog$assets$url <- paste0("file:///", release, "/", catalog$assets$name)
catalog$assets$historical <- TRUE
proof <- readRDS(file.path(imported, "bf-import.rds"))
expect_error <- function(expr, label) {
    result <- tryCatch({ force(expr); FALSE }, error = function(e) {
        cat(label, "rejected:", conditionMessage(e), "\n"); TRUE
    })
    stopifnot(result)
}
expect_error(bfpwr_sim_moment_official_bf_proof(archive, proof$chunks), "Unpublished metadata digest")
official <- bfpwr_sim_moment_official_bf_proof(archive, proof$chunks, catalog)
bad <- proof$chunks
bad$output_sha256[1] <- paste(rep("0", 64), collapse = "")
expect_error(bfpwr_sim_moment_official_bf_proof(archive, bad, catalog), "Different BF output")
expect_error(bfpwr_sim_moment_official_bf_proof(archive, proof$chunks[1, , drop = FALSE], catalog),
    "Incomplete BF key set")
bad <- proof$chunks
bad$case_id[1] <- "other-case"
expect_error(bfpwr_sim_moment_official_bf_proof(archive, bad, catalog), "Different BF keys")
archive_raw <- readBin(archive, "raw", file.info(archive)$size)
writeBin(c(archive_raw, as.raw(0)), archive)
expect_error(bfpwr_sim_moment_official_bf_proof(archive, proof$chunks, catalog), "Damaged metadata archive")
writeBin(archive_raw, archive)
## Override the fixed catalog only inside this disposable source checkout.
catalog_path <- file.path(scratch, "private-catalog.rds")
saveRDS(catalog, catalog_path)
cat("\nbfpwr_sim_moment_v1_catalog <- function() readRDS(", dQuote(catalog_path, FALSE), ")\n", sep = "",
    file = file.path(root, "simulations/R/moment-identity.R"), append = TRUE)
official_output <- file.path(scratch, "official-fixture-v2")
run("plan", official_output)
run("import", official_output, c("--reuse-from", shQuote(old), "--bf-source-root", shQuote(historical),
    "--official-metadata", shQuote(archive)))
stopifnot(identical(readRDS(file.path(official_output, "bf-import.rds"))$official_payload, official))
run("summaries", official_output)
run("finalize", official_output)
run("assemble", official_output)
for (file in c("validation-checks.csv", "numerical-validation-checks.csv")) {
    stopifnot(file.copy(file.path(imported, file), file.path(official_output, file)))
}
official_release <- file.path(scratch, "official-fixture-release")
run("package", official_output, c("--release", shQuote(official_release)), script = "package_moment_extension.R")
official_manifest <- read.csv(file.path(official_release, "moment-v2-manifest.csv"))
stopifnot(sum(official_manifest$historical) == 1L, sum(!official_manifest$historical) == 3L,
    !any(grepl("bf-chunks", list.files(official_release), fixed = TRUE)))
official_download <- file.path(scratch, "official-fixture-download")
run("download", official_download, c("--with-bf-chunks",
    "--release-base", shQuote(paste0("file:///", official_release, "/")),
    "--downloads", shQuote(file.path(scratch, "official-fixture-downloads"))), script = "download_moment_extension.R")
stopifnot(identical(unname(tools::sha256sum(bfpwr_sim_bf_chunk_files(official_output, case, design))),
                    unname(tools::sha256sum(bfpwr_sim_bf_chunk_files(official_download, case, design)))))
## Even a new release destination must reject damage after a successful proof.
path <- bfpwr_sim_bf_chunk_files(official_output, case, design)[1]
raw <- readBin(path, "raw", file.info(path)$size)
writeBin(c(raw, as.raw(0)), path)
run("package", official_output, c("--release", shQuote(file.path(scratch, "damaged-official-release"))),
    script = "package_moment_extension.R", expect = 1L)
status_path <- file.path(official_output, "chunk-status/0001.rds")
status_raw <- readBin(status_path, "raw", file.info(status_path)$size)
status <- readRDS(status_path)
status$output_sha256 <- unname(tools::sha256sum(path))
saveRDS(status, status_path)
run("package", official_output, c("--release", shQuote(file.path(scratch, "different-official-release"))),
    script = "package_moment_extension.R", expect = 1L)
writeBin(raw, path)
writeBin(status_raw, status_path)
cat("Optional official proof: digest, complete keys, output hashes, historical links/download and post-proof damage checks passed.\n")
cat("Engineering regression checks passed at", scratch, "\n")
