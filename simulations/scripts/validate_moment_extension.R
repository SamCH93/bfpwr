script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(script), "simulation_support.R"))
source_simulation_library()
source(file.path(repo_root(),"simulations/R/moment-extension.R"))
args <- parse_args()
output <- arg_value(args,"output",file.path(repo_root(),"simulations/corpus/moment-v2"))
x <- readRDS(file.path(output,"verification.rds"))
expected_cases <- bfpwr_sim_moment_extension_cases()
ids <- vapply(expected_cases,`[[`,"","bf_prior_id")
stopifnot(identical(vapply(x$cases,`[[`,"","bf_prior_id"),ids))
stopifnot(identical(x$cases, expected_cases))
designs <- lapply(expected_cases,function(p) bfpwr_sim_find_design_case(p$design_case_id))
schedules <- bfpwr_sim_z_fixture_schedules_core()
compatible <- lapply(designs,function(d) schedules[vapply(schedules,
    bfpwr_sim_schedule_is_compatible,logical(1),design=d)])
expected_fixed <- sum(vapply(designs,function(d) 6*length(d$look_grid),numeric(1)))
expected_seq <- sum(vapply(compatible,function(s) 9*sum(vapply(s,function(y) length(y$n),integer(1))),numeric(1)))
expected_moments <- sum(lengths(compatible))*3
stopifnot(nrow(x$fixed)==expected_fixed,nrow(x$sequential)==expected_seq,
          nrow(x$moments)==expected_moments,nrow(x$searches)==4200,
          x$provenance$bf_rows==244000000,
          setequal(unique(x$fixed$case_id),ids),setequal(unique(x$sequential$case_id),ids),
          all(x$fixed$nsim==10000),all(x$sequential$nsim==10000))
search_failures <- bfpwr_sim_moment_search_failures(x$searches)
checks <- data.frame(check=c("Fixed probabilities valid", "Sequential probabilities valid",
    "Fixed severe MC intervals", "Sequential severe MC intervals", "Expected sample size",
    "Sample-size variance", "Search runtime errors", "Returned size attains target",
    "Returned size is not later than first crossing", "No missed crossing (fixed lower-bound exception)"),
    failures=c(x$fixed_diagnostics$summary$invalid_probabilities,
        x$sequential_diagnostics$summary$invalid_probabilities,
        x$fixed_diagnostics$summary$severe_interval_failures,
        x$sequential_diagnostics$summary$severe_interval_failures,
        sum(!x$moments$EN_ok),sum(!x$moments$VarN_ok),unname(search_failures)))
write.csv(checks,file.path(output,"validation-checks.csv"),row.names=FALSE)
print(checks,row.names=FALSE)
cat("Coverage:",length(ids),"conditions;",expected_fixed,"fixed probabilities;",
    expected_seq,"sequential probabilities;",expected_moments,"stopping distributions;",
    nrow(x$searches),"searches.\n")
## Reference convergence is separate from the preserved MC/search screening.
if (!is.null(x$numerical_diagnostics)) {
    numerical <- x$numerical_diagnostics$moment_rows
    numerical_checks <- data.frame(check = "Independent reference convergence",
        failures = sum(is.na(numerical$reference_converged) | !numerical$reference_converged))
    write.csv(numerical_checks, file.path(output, "numerical-validation-checks.csv"), row.names = FALSE)
    print(numerical_checks, row.names = FALSE)
    if (any(numerical_checks$failures > 0)) stop("Independent reference did not converge; MC screening is separate.")
} else if (identical(x$timing_schema_version, 2L)) stop("v2 numerical diagnostics are missing")
if(any(checks$failures>0)) stop("Validation discrepancies require inspection; see saved checks.")
