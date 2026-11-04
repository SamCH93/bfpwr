script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source_simulation_library()
args <- parse_args()
output <- arg_value(args, "output", file.path(repo_root(), "simulations/corpus/moment-v2"))
files <- list.files(file.path(output,"case-results"), pattern="[.]rds$", full.names=TRUE)
results <- lapply(files, readRDS)
cat("Completed cases:",length(results),"\n")
for (mode in c("fixed","sequential")) {
    rows <- do.call(rbind,lapply(results,`[[`,mode))
    diagnostics <- bfpwr_sim_mc_reference_diagnostics(rows,label=mode,
        group_cols=if(mode=="fixed") "tail" else "event")
    print(diagnostics$summary[c("label","rows_checked","invalid_probabilities",
        "severe_interval_failures","holm_failures","max_abs_error")])
}
moments <- do.call(rbind,lapply(results,`[[`,"moments"))
searches <- do.call(rbind,lapply(results,`[[`,"searches"))
cat("Moment failures:",sum(!moments$EN_ok),sum(!moments$VarN_ok),"\n")
print(do.call(rbind,lapply(split(searches,searches$mode),function(x) data.frame(
    mode=x$mode[1], calls=nrow(x), errors=sum(nzchar(x$error)),
    returned_below_target=sum(is.finite(x$n_found)&!x$attains_target),
    later_crossing=sum(is.finite(x$n_found)&is.finite(x$minimum_n)&x$n_found>x$minimum_n),
    lower_bound_warning=sum(!is.finite(x$n_found)&x$minimum_n==10,na.rm=TRUE),
    missed_interior_crossing=sum(!is.finite(x$n_found)&is.finite(x$minimum_n)&x$minimum_n>10)))))
