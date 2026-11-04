## Compatibility entry point for the standard BF chunk validator/index writer.
if (sys.nframe() == 0L) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
    root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
    source(file.path(root, "simulations/scripts/run_moment_extension.R"))
    bfpwr_sim_run_moment_extension(c("--phase", "finalize", commandArgs(TRUE)), root)
}
