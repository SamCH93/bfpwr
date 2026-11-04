## Render the pooled full report and new-priors-only incremental report from cached inputs.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
args <- parse_args()
root <- repo_root()
scope <- arg_value(args, "scope", "both")
stopifnot(scope %in% c("full", "incremental", "both"))
baseline <- arg_value(args, "baseline-root", "")
if (!nzchar(baseline)) {
    stop("Supply --baseline-root with the complete saved baseline report inputs")
}
baseline <- normalizePath(baseline, winslash = "/", mustWork = TRUE)
moment <- resolve_moment_root(root, arg_value(args, "moment-root", arg_value(args, "results")))
destination <- normalizePath(arg_value(args, "pdf-dir", file.path(root, "output/pdf")),
                             winslash = "/", mustWork = FALSE)
build <- normalizePath(arg_value(args, "build-dir",
    file.path(destination, "build/moment-extension")), winslash = "/", mustWork = FALSE)
Sys.setenv(BFPWR_SIM_REPO = root, BFPWR_SIM_CORPUS = baseline, BFPWR_SIM_RESULTS = baseline,
           BFPWR_SIM_MOMENT_RESULTS = moment)
template <- readLines(file.path(root, "simulations/simulation-verification.Rnw"), warn = FALSE)
start <- grep("^<<fixture-setup,", template)
end <- which(seq_along(template) > start & template == "@")[1]
setup <- template[(start + 1):(end - 1)]
eval(parse(text = setup), envir = .GlobalEnv)
cat("Report input preflight passed:", length(fixtures), "report fixtures;",
    length(moment_bundle$cases), "new moment cases\n")
if (isTRUE(arg_flag(args, "preflight-only"))) quit(status = 0L)

ensure_dir(destination)
destination <- normalizePath(destination, winslash = "/", mustWork = TRUE)
ensure_dir(build)
build <- normalizePath(build, winslash = "/", mustWork = TRUE)
inventory <- function(directory, files) {
    files <- sort(unique(files), method = "radix")
    hashes <- tools::sha256sum(file.path(directory, files))
    stopifnot(!anyNA(hashes))
    setNames(unname(hashes), files)
}
baseline_files <- substring(baseline_inputs, nchar(baseline) + 2L)
for (directory in c("fixtures", "fixture-validation", "search-validation")) {
    files <- list.files(file.path(baseline, directory), "[.](rds|csv|json)$",
                        all.files = TRUE, recursive = TRUE)
    baseline_files <- c(baseline_files, file.path(directory, files))
}
moment_files <- c("verification.rds", "numerical-diagnostics.rds", "provenance.rds")
moment_paths <- file.path(moment, moment_files)
moment_files <- moment_files[file.exists(moment_paths) & !dir.exists(moment_paths)]
source_files <- c("package/DESCRIPTION", "simulations/simulation-verification.Rnw",
    "simulations/moment-extension-body.Rnw",
    "simulations/registry/manifests/package-verification-cases.csv",
    "simulations/scripts/simulation_support.R", "simulations/scripts/moment_report_adapters.R",
    "simulations/scripts/render_moment_extension.R")
input_provenance <- list(baseline = reference_provenance, normal = one_sided$provenance,
    moment = moment_bundle$provenance, baseline_root = baseline, moment_root = moment,
    scope = scope, baseline_input_sha256 = inventory(baseline, baseline_files),
    moment_input_sha256 = inventory(moment, moment_files), source_sha256 = inventory(root, source_files),
    renderer = list(R = R.version.string, platform = R.version$platform,
                    libraries = .libPaths(), source_root = root))
old <- setwd(build)
on.exit(setwd(old), add = TRUE)
options(knitr.include_graphics.ext = TRUE)
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, error = FALSE,
    fig.height = 6.8, fig.width = 6.8, fig.align = "center", dev = "png", dpi = 300,
    out.width = "\\textwidth")
if (scope != "incremental") {
    writeLines(template, "simulation-verification-full.Rnw")
    knitr::knit2pdf("simulation-verification-full.Rnw", clean = FALSE, quiet = TRUE, envir = .GlobalEnv)
    stopifnot(file.copy("simulation-verification-full.pdf",
        file.path(destination, "simulation-verification-full.pdf"), overwrite = TRUE))
}

if (scope != "full") {
    .fixture_table_counter <- 0L
    preamble <- template[seq_len(grep("^\\\\begin\\{document\\}", template) - 1L)]
    body <- readLines(file.path(root, "simulations/moment-extension-body.Rnw"), warn = FALSE)
    settings <- c("<<moment-render-options, include=FALSE>>=",
        "knitr::opts_chunk$set(dev='png', dpi=300, echo=FALSE, message=FALSE, warning=FALSE, error=FALSE, fig.width=6.8, fig.height=6.8, fig.align='center', out.width='\\\\textwidth')", "@")
    writeLines(c(preamble, "\\begin{document}", settings, "\\maketitle",
        "\\noindent This report compares \\textbf{bfpwr} predictions with simulations",
        "across the listed priors, evidence thresholds, and sampling schedules.",
        body, "\\end{document}"), "simulation-verification-incremental.Rnw")
    knitr::knit2pdf("simulation-verification-incremental.Rnw", clean = FALSE, quiet = TRUE, envir = .GlobalEnv)
    stopifnot(file.copy("simulation-verification-incremental.pdf",
        file.path(destination, "simulation-verification-incremental.pdf"), overwrite = TRUE))
}
scopes <- if (scope == "both") c("full", "incremental") else scope
pdf_files <- paste0("simulation-verification-", scopes, ".pdf")
input_provenance$output_sha256 <- inventory(destination, pdf_files)
saveRDS(input_provenance, file.path(destination, "report-input-provenance.rds"))
cat("Rendered", scope, "report scope in", destination, "\n")
