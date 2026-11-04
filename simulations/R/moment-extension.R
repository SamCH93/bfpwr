## PR #8: add directions to the existing 35 moment-prior cases.
## Trajectories, BF evaluation, fixture summaries, and MC diagnostics are reused.
bfpwr_sim_moment_extension_cases <- function() {
    base <- bfpwr_sim_bf_prior_cases_z(bfpwr_sim_design_cases())
    base <- base[vapply(base, function(x) x$bf_type == "moment", logical(1))]
    stopifnot(length(base) == 35L)
    cases <- list()
    for (original in base) for (alternative in c("less", "greater")) {
        x <- original
        x$analysis_prior$alternative <- alternative
        x$bf_prior_id <- sub("-on-", paste0("-", alternative, "-on-"),
                             sub("-(less|greater)-on-", "-on-", x$bf_prior_id))
        x$tags <- unique(c(x$tags, "one-sided-moment", alternative))
        cases[[length(cases) + 1L]] <- x
    }
    stopifnot(length(cases) == 70L, !anyDuplicated(vapply(cases, `[[`, "", "bf_prior_id")))
    cases
}

bfpwr_sim_moment_extension_plan <- function(cases = bfpwr_sim_moment_extension_cases()) {
    do.call(rbind, lapply(seq_along(cases), function(i) {
        design <- bfpwr_sim_find_design_case(cases[[i]]$design_case_id)
        data.frame(case_index = i, bf_prior_id = cases[[i]]$bf_prior_id,
                   design_case_id = design$design_case_id,
                   chunk_id = design$chunks$chunk_id)
    }))
}

bfpwr_sim_moment_extension_chunk <- function(input, output, case, chunk_id) {
    design <- bfpwr_sim_find_design_case(case$design_case_id)
    path <- bfpwr_sim_design_chunk_files(input, design)[chunk_id]
    source_record <- bfpwr_sim_moment_verify_input(input, design, chunk_id)
    trajectories <- readRDS(path)
    chunk <- design$chunks[design$chunks$chunk_id == chunk_id, ]
    stopifnot(nrow(trajectories) == (chunk$replicate_end - chunk$replicate_start + 1)*length(design$look_grid),
              all(trajectories$design_case_id == design$design_case_id),
              setequal(trajectories$replicate_id, seq(chunk$replicate_start, chunk$replicate_end)),
              all(is.finite(trajectories$estimate)),
              all(abs(trajectories$se - design$generation$usd/sqrt(trajectories$n)) < 1e-12))
    expected <- expand.grid(n = design$look_grid,
                            replicate_id = seq(chunk$replicate_start, chunk$replicate_end))
    key <- function(x) paste(x$replicate_id, x$n)
    stopifnot(!anyDuplicated(key(trajectories)), setequal(key(trajectories), key(expected)))
    rows <- bfpwr_sim_bf_rows(case, trajectories)
    stopifnot(!anyNA(rows$log_bf01))
    target <- bfpwr_sim_bf_chunk_files(output, case, design)[chunk_id]
    bfpwr_sim_write_rds(rows, paste0(target, ".tmp"))
    stopifnot(file.rename(paste0(target, ".tmp"), target))
    list(input = file.path("designs", design$test_family, design$design_case_id,
                           source_record$file), input_path_base = "provenance$source_corpus",
         input_sha256 = source_record$sha256[[1]],
         input_md5 = unname(tools::md5sum(path)),
         output_sha256 = unname(tools::sha256sum(target)),
         output_md5 = unname(tools::md5sum(target)), rows = nrow(rows))
}

## A single prediction supplies every nested schedule in its family.
bfpwr_sim_moment_prediction_key <- function(family, threshold) paste(family, threshold, sep = ":")

bfpwr_sim_moment_mc_fields <- function(m, cumulative, absolute_denominator = FALSE) {
    if (!"max_observed_n" %in% names(m)) {
        key <- function(x) paste(x$case_id, x$schedule_id, x$evidence_threshold, sep = "|")
        observed <- cumulative[cumulative$stop_H1_at_look + cumulative$stop_H0_at_look > 0, ]
        maximum <- tapply(observed$n, key(observed), max)
        m$max_observed_n <- unname(maximum[key(m)])
        m$max_observed_n[m$nInc > 0] <- m$max_n[m$nInc > 0]
        if (anyNA(m$max_observed_n)) stop("Cannot recover legacy MC tolerance scale")
    }
    for (quantity in c("EN", "VarN")) {
        power <- if (quantity == "EN") 1 else 2
        if (!paste0("tolerance_", quantity) %in% names(m)) {
            m[[paste0("tolerance_", quantity)]] <- pmax(5*m[[paste0("mcse_", quantity)]],
                                                      .002*m$max_observed_n^power)
        }
        observed <- m[[quantity]]
        error <- m[[paste0("reference_", quantity)]] - observed
        ## Preserve absolute denominators for runner rows and signed ones for diagnostics,
        ## including when the observed VarN is slightly negative.
        denominator <- if (absolute_denominator) abs(observed) else observed
        m[[paste0("error_", quantity)]] <- error
        m[[paste0("relative_error_", quantity)]] <- ifelse(observed != 0,
                                                               error/denominator, NA_real_)
        m[[paste0("z_", quantity)]] <- ifelse(m[[paste0("mcse_", quantity)]] > 0,
            error/m[[paste0("mcse_", quantity)]], NA_real_)
    }
    m
}

bfpwr_sim_moment_schedule_row <- function(final, cumulative, cumpH1, cumpH0, n) {
    counts <- cumulative$stop_H1_at_look + cumulative$stop_H0_at_look
    counts[length(counts)] <- counts[length(counts)] + final$nInc
    values <- rep(cumulative$n, counts)
    ## The package moment helper assigns unresolved mass to this prefix's last n.
    reference <- .bfseq_sample_size_moments(diff(c(0, cumpH1)), diff(c(0, cumpH0)), n)
    final$reference_EN <- reference$EN
    final$reference_VarN <- reference$VarN
    final$mcse_EN <- sd(values)/sqrt(length(values))
    final$mcse_VarN <- sd((values - mean(values))^2)/sqrt(length(values))
    final$max_observed_n <- max(values)
    mc <- bfpwr_sim_moment_mc_fields(final, cumulative, absolute_denominator = TRUE)
    ## Assign fields explicitly to preserve the legacy runner column order.
    fields <- c("tolerance_EN", "tolerance_VarN", "error_EN", "error_VarN",
        "relative_error_EN", "relative_error_VarN", "z_EN", "z_VarN")
    final[fields] <- mc[fields]
    final$EN_ok <- abs(final$error_EN) <= final$tolerance_EN
    final$VarN_ok <- abs(final$error_VarN) <= final$tolerance_VarN
    final
}

bfpwr_sim_moment_extension_case <- function(output, case, max_looks = Inf) {
    started <- Sys.time()
    case_start <- proc.time()[["elapsed"]]
    phases <- numeric()
    timings <- reuse <- list()
    record <- function(operation, seconds, phase, schedule_id = NA_character_,
                       evidence_threshold = NA_real_, n_looks = NA_integer_, rows_covered = 1L) {
        timings[[length(timings) + 1L]] <<- data.frame(operation = operation,
            seconds = seconds, phase = phase, schedule_id = schedule_id,
            evidence_threshold = evidence_threshold, n_looks = n_looks,
            rows_covered = rows_covered, stringsAsFactors = FALSE)
    }
    t0 <- proc.time()[["elapsed"]]
    design <- bfpwr_sim_find_design_case(case$design_case_id)
    prior <- case$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    arguments <- list(usd = design$generation$usd, null = prior$null,
                      psd = prior$psd, dpm = dp$dpm, dpsd = dp$dpsd,
                      alternative = prior$alternative)
    thresholds <- bfpwr_sim_fixture_evidence_thresholds_core()
    schedules <- bfpwr_sim_z_fixture_schedules_core()
    schedules <- schedules[vapply(schedules, function(s) {
        s$n_looks <= max_looks && bfpwr_sim_schedule_is_compatible(s, design)
    }, logical(1))]
    phases["setup_seconds"] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    fixed <- bfpwr_sim_summarize_fixed_tails(output, case, design, thresholds,
        fixture_set_id = "z-nmbf01-one-sided-fixed-v1")
    decisions <- bfpwr_sim_summarize_fixed_decisions(fixed,
        bfpwr_sim_threshold_pairs(thresholds))
    phases["fixed_empirical_seconds"] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    fixed$reference_prob <- do.call(pnmbf01, c(list(k = fixed$threshold,
        n = fixed$n, lower.tail = fixed$tail == "H1"), arguments))
    record("pnmbf01", proc.time()[["elapsed"]] - t0, "fixed_reference",
           n_looks = length(design$look_grid), rows_covered = nrow(fixed))
    phases["fixed_reference_seconds"] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    sequential <- bfpwr_sim_summarize_sequential(output, case, design, schedules,
        thresholds, fixture_set_id = "z-nmbf01-one-sided-sequential-v1")
    phases["sequential_empirical_seconds"] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    predictions <- probabilities <- moments <- list()
    families <- split(schedules, vapply(schedules, `[[`, "", "schedule_family_id"))
    for (family in families) {
        longest <- family[[which.max(vapply(family, function(s) length(s$n), integer(1)))]]
        for (schedule in family) {
            if (!identical(schedule$n, longest$n[seq_along(schedule$n)])) {
                stop("schedule family is not nested: ", schedule$schedule_family_id)
            }
        }
        for (threshold in thresholds) {
            call_start <- proc.time()[["elapsed"]]
            prediction <- do.call(pbf01seq, c(list(k1 = 1/threshold, k0 = threshold,
                n = longest$n, se = arguments$usd/sqrt(longest$n), type = "moment"),
                arguments[setdiff(names(arguments), "usd")]))
            record("pbf01seq", proc.time()[["elapsed"]] - call_start, "sequential_reference",
                longest$schedule_id, threshold, length(longest$n), 3L*length(longest$n))
            key <- bfpwr_sim_moment_prediction_key(longest$schedule_family_id, threshold)
            predictions[[key]] <- list(prediction = prediction, schedule = longest)
            for (schedule in family) {
                index <- seq_along(schedule$n)
                x <- sequential$cumulative_summary[
                    sequential$cumulative_summary$schedule_id == schedule$schedule_id &
                    sequential$cumulative_summary$evidence_threshold == threshold, ]
                x <- x[order(x$n), ]
                stopifnot(nrow(x) == length(index), identical(x$n, schedule$n))
                for (event in c("H1", "H0", "Inc")) {
                    rows <- x
                    rows$event <- event
                    rows$prob <- x[[paste0("cum_p", event)]]
                    rows$n_event <- round(rows$prob*x$nsim)
                    rows$reference_prob <- prediction[[paste0("cump", event)]][index]
                    probabilities[[length(probabilities) + 1L]] <- rows
                }
                if (schedule$schedule_id != longest$schedule_id) {
                    reuse[[length(reuse) + 1L]] <- data.frame(phase = "sequential_reference",
                        schedule_id = schedule$schedule_id, source_schedule_id = longest$schedule_id,
                        evidence_threshold = threshold, n_looks = length(index),
                        rows_covered = 3L*length(index), stringsAsFactors = FALSE)
                }
                final <- sequential$final_summary[
                    sequential$final_summary$schedule_id == schedule$schedule_id &
                    sequential$final_summary$evidence_threshold == threshold, ]
                moments[[length(moments) + 1L]] <- bfpwr_sim_moment_schedule_row(
                    final, x, prediction$cumpH1[index], prediction$cumpH0[index], schedule$n)
            }
        }
    }
    phases["sequential_reference_seconds"] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    searches <- bfpwr_sim_moment_extension_search(case, design, fixed, sequential,
                                                 predictions = predictions, record = record)
    search_reuse <- attr(searches, "reuse")
    attr(searches, "reuse") <- NULL
    reuse <- c(reuse, search_reuse)
    phases["search_seconds"] <- proc.time()[["elapsed"]] - t0
    case_timings <- data.frame(started = started, completed = Sys.time(),
        elapsed_seconds = proc.time()[["elapsed"]] - case_start, as.list(phases))
    list(case = case, design = design, fixed = fixed, decisions = decisions,
         cumulative = sequential$cumulative_summary, final = sequential$final_summary,
         sequential = do.call(rbind, probabilities), moments = do.call(rbind, moments),
         searches = searches, timings = do.call(rbind, timings),
         reuse = if (length(reuse)) do.call(rbind, reuse) else data.frame(),
         case_timings = case_timings, timing_schema_version = 2L)
}

bfpwr_sim_moment_extension_search <- function(case, design, fixed, sequential,
                                             predictions = list(), record = function(...) NULL) {
    p <- case$analysis_prior
    dp <- bfpwr_sim_design_prior_mean_sd(design)
    args <- list(usd = design$generation$usd, null = p$null, psd = p$psd,
                 dpm = dp$dpm, dpsd = dp$dpsd, alternative = p$alternative)
    rows <- reuse <- list()
    targets <- bfpwr_sim_fixture_target_probabilities_core()
    max_n <- if (design$look_grid_name == "long") 1000 else 500
    for (mode in c("fixed", "sequential")) for (threshold in c(3, 10, 30)) {
        for (event in c("H1", "H0")) {
            k <- if (event == "H1") 1/threshold else threshold
            schedule_id <- NA_character_
            if (mode == "fixed") {
                candidate_n <- 10:max_n
                scan_start <- proc.time()[["elapsed"]]
                candidate_power <- do.call(pnmbf01, c(list(k = k, n = candidate_n,
                    lower.tail = event == "H1"), args))
                record("pnmbf01", proc.time()[["elapsed"]] - scan_start, "search_reference_scan",
                    evidence_threshold = threshold, n_looks = length(candidate_n),
                    rows_covered = length(candidate_n))
                empirical <- fixed[fixed$evidence_threshold == threshold &
                    fixed$tail == event & fixed$n <= max_n, ]
                empirical <- empirical[order(empirical$n), ]
            } else {
                schedule_id <- if (max_n == 1000) "start10-by10-looks100" else "start10-by10-looks50"
                empirical <- sequential$cumulative_summary[
                    sequential$cumulative_summary$schedule_id == schedule_id &
                    sequential$cumulative_summary$evidence_threshold == threshold, ]
                empirical <- empirical[order(empirical$n), ]
                if (!nrow(empirical)) next
                candidate_n <- empirical$n
                key <- bfpwr_sim_moment_prediction_key("start10-by10", threshold)
                cached <- predictions[[key]]
                if (is.null(cached) || !identical(candidate_n, cached$schedule$n[seq_along(candidate_n)])) {
                    stop("missing compatible sequential reference curve for search: ", schedule_id)
                }
                candidate_power <- cached$prediction[[paste0("cump", event)]][seq_along(candidate_n)]
                reuse[[length(reuse) + 1L]] <- data.frame(phase = "search_reference_scan",
                    schedule_id = schedule_id, source_schedule_id = cached$schedule$schedule_id,
                    evidence_threshold = threshold, n_looks = length(candidate_n),
                    rows_covered = length(candidate_n), stringsAsFactors = FALSE)
                empirical$prob <- empirical[[paste0("cum_p", event)]]
                empirical$mcse <- empirical[[paste0("mcse_cum_p", event)]]
            }
            for (target in targets) {
                warning <- error <- ""
                t0 <- proc.time()[["elapsed"]]
                found <- tryCatch(withCallingHandlers({
                    if (mode == "fixed") do.call(nnmbf01, c(list(k = k, power = target,
                        nrange = c(10, max_n), lower.tail = event == "H1"), args))
                    else do.call(nbf01seq, c(list(k1 = 1/threshold, k0 = threshold,
                        power = target, target = event, type = "moment", minN = 10,
                        by = 10, nrange = c(10, max_n)), args))
                }, warning = function(w) {
                    warning <<- paste(warning, conditionMessage(w)); invokeRestart("muffleWarning")
                }), error = function(e) { error <<- conditionMessage(e); NA_real_ })
                record(if (mode == "fixed") "nnmbf01" else "nbf01seq",
                    proc.time()[["elapsed"]] - t0, "search_solver", schedule_id,
                    threshold, length(candidate_n))
                hit <- which(candidate_power >= target)
                first_n <- if (length(hit)) candidate_n[hit[1]] else NA_real_
                sim_hit <- which(empirical$prob >= target)
                grid_n <- if (is.finite(found)) ceiling(found/10)*10 else NA_real_
                j <- match(grid_n, empirical$n)
                actual <- NA_real_
                if (is.finite(found) && mode == "fixed") {
                    attainment_start <- proc.time()[["elapsed"]]
                    actual <- do.call(pnmbf01, c(list(k = k, n = found,
                        lower.tail = event == "H1"), args))
                    record("pnmbf01", proc.time()[["elapsed"]] - attainment_start,
                        "search_attainment", evidence_threshold = threshold, n_looks = 1L)
                } else if (is.finite(found)) actual <- candidate_power[match(found, candidate_n)]
                rows[[length(rows) + 1L]] <- data.frame(mode = mode,
                    threshold = threshold, event = event, target = target, n_found = found,
                    minimum_n = first_n, attains_target = is.finite(actual) && actual >= target - 1e-7,
                    minimum_correct = if (is.finite(first_n)) isTRUE(found == first_n) else !is.finite(found),
                    actual_power = actual, simulation_n = if (length(sim_hit)) empirical$n[sim_hit[1]] else NA_real_,
                    grid_n = grid_n, simulation_power = empirical$prob[j], simulation_mcse = empirical$mcse[j],
                    seconds = proc.time()[["elapsed"]] - t0, warning = trimws(warning), error = error)
            }
        }
    }
    result <- do.call(rbind, rows)
    attr(result, "reuse") <- reuse
    result
}

## The lower-bound warning exception belongs only to the fixed search.
bfpwr_sim_moment_search_failures <- function(searches) {
    c(runtime_errors = sum(nzchar(searches$error)),
      target_not_attained = sum(is.finite(searches$n_found) & !searches$attains_target),
      returned_after_first_crossing = sum(is.finite(searches$n_found) &
          is.finite(searches$minimum_n) & searches$n_found > searches$minimum_n),
      missed_crossing = sum(!is.finite(searches$n_found) & is.finite(searches$minimum_n) &
          (searches$mode != "fixed" | searches$minimum_n > 10)))
}
