support_path <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    marker <- "--file="
    hit <- grep(paste0("^", marker), cmd, value = TRUE)
    script <- if (length(hit) > 0) {
        normalizePath(sub(paste0("^", marker), "", hit[[1]]),
                      winslash = "/", mustWork = TRUE)
    } else {
        normalizePath("simulations/scripts/render_simulation_verification_background.R",
                      winslash = "/", mustWork = FALSE)
    }
    file.path(dirname(script), "simulation_support.R")
}

source(support_path())
source_simulation_library()

args <- parse_args()
root <- repo_root()
corpus_root <- normalizePath(arg_value(args, "corpus-root",
                                       file.path(root, "simulations",
                                                 "corpus", "v1")),
                             winslash = "/", mustWork = FALSE)
corpus_root <- bfpwr_sim_require_corpus_root(
    corpus_root,
    required_paths = c("fixtures", "fixture-validation"),
    asset_set = "fixture-validation",
    purpose = "simulation verification report rendering")

results_root <- normalizePath(arg_value(args, "results-root", corpus_root),
                              winslash = "/", mustWork = TRUE)
stopifnot(dir.exists(file.path(results_root, "fixture-validation")))
Sys.setenv(
    BFPWR_SIM_CORPUS = corpus_root,
    BFPWR_SIM_RESULTS = results_root,
    BFPWR_SIM_REPO = root
)

cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Repo:", getwd(), "\n")
cat("Corpus:", Sys.getenv("BFPWR_SIM_CORPUS"), "\n")
cat("Package results:", results_root, "\n")

t0 <- proc.time()[["elapsed"]]
old <- setwd(file.path(root, "simulations"))
on.exit(setwd(old), add = TRUE)
knitr::knit2pdf("simulation-verification.Rnw", clean = TRUE, quiet = TRUE)

cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Elapsed seconds:", proc.time()[["elapsed"]] - t0, "\n")
