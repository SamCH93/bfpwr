bfpwr_sim_fixture_dir <- function(corpus_root, fixture) {
    file.path(corpus_root, "fixtures", fixture$family, fixture$mode,
              fixture$fixture_set_id)
}

bfpwr_sim_fixture_dir_from_id <- function(corpus_root,
                                          fixture_set_id,
                                          family,
                                          mode = c("fixed", "sequential")) {
    mode <- match.arg(mode)
    file.path(corpus_root, "fixtures", family, mode, fixture_set_id)
}

bfpwr_sim_read_fixture_summary <- function(corpus_root,
                                           fixture_set_id,
                                           family = "z",
                                           mode = c("fixed", "sequential")) {
    mode <- match.arg(mode)
    fixture_dir <- bfpwr_sim_fixture_dir_from_id(
        corpus_root = corpus_root,
        fixture_set_id = fixture_set_id,
        family = family,
        mode = mode)
    if (!dir.exists(fixture_dir)) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = fixture_dir,
            asset_set = "fixture-tests",
            label = "fixture summary directories")
    }
    read_file <- function(name) {
        file <- file.path(fixture_dir, name)
        if (!file.exists(file)) {
            bfpwr_sim_stop_missing_corpus_files(
                corpus_root,
                missing = file,
                asset_set = "fixture-tests",
                label = "fixture summary artifacts")
        }
        readRDS(file)
    }
    out <- list(
        fixture_dir = fixture_dir,
        spec = read_file("spec.rds"),
        manifest = utils::read.csv(file.path(fixture_dir, "manifest.csv"),
                                   stringsAsFactors = FALSE)
    )
    if (mode == "fixed") {
        out$tail_summary <- read_file("tail-summary.rds")
        out$decision_summary <- read_file("decision-summary.rds")
        out$search_summary <- read_file("search-summary.rds")
    } else {
        out$cumulative_summary <- read_file("cumulative-summary.rds")
        out$final_summary <- read_file("final-summary.rds")
        out$search_summary <- read_file("search-summary.rds")
    }
    out
}

bfpwr_sim_mcse <- function(prob, nsim) {
    sqrt(prob * (1 - prob) / nsim)
}

bfpwr_sim_mc_tolerance <- function(reference_prob,
                                   nsim,
                                   observed_mcse = NULL,
                                   z = 4,
                                   floor = 0.002) {
    ref_mcse <- bfpwr_sim_mcse(reference_prob, nsim)
    if (is.null(observed_mcse)) {
        observed_mcse <- ref_mcse
    }
    pmax(floor, z * pmax(ref_mcse, observed_mcse))
}

bfpwr_sim_mc_close <- function(observed,
                               reference,
                               nsim,
                               observed_mcse = NULL,
                               z = 4,
                               floor = 0.002) {
    abs(observed - reference) <= bfpwr_sim_mc_tolerance(
        reference_prob = reference,
        nsim = nsim,
        observed_mcse = observed_mcse,
        z = z,
        floor = floor)
}

bfpwr_sim_tail_thresholds <- function(evidence_thresholds) {
    evidence_thresholds <- as.numeric(evidence_thresholds)
    if (length(evidence_thresholds) == 0 ||
        any(!is.finite(evidence_thresholds)) ||
        any(evidence_thresholds <= 1)) {
        stop("evidence_thresholds must contain finite values greater than 1")
    }
    evidence_thresholds <- sort(unique(evidence_thresholds))

    rows <- lapply(evidence_thresholds, function(k) {
        k_id <- bfpwr_sim_decimal_id(k)
        rbind(
            data.frame(
                threshold_id = paste0("h1-bf", k_id),
                evidence_threshold = k,
                tail = "H1",
                threshold = 1 / k,
                log_threshold = log(1 / k),
                stringsAsFactors = FALSE
            ),
            data.frame(
                threshold_id = paste0("h0-bf", k_id),
                evidence_threshold = k,
                tail = "H0",
                threshold = k,
                log_threshold = log(k),
                stringsAsFactors = FALSE
            )
        )
    })
    do.call(rbind, rows)
}

bfpwr_sim_threshold_pairs <- function(evidence_thresholds) {
    evidence_thresholds <- as.numeric(evidence_thresholds)
    if (length(evidence_thresholds) == 0 ||
        any(!is.finite(evidence_thresholds)) ||
        any(evidence_thresholds <= 1)) {
        stop("evidence_thresholds must contain finite values greater than 1")
    }
    evidence_thresholds <- sort(unique(evidence_thresholds))

    data.frame(
        threshold_pair_id = paste0("bf", bfpwr_sim_decimal_id(evidence_thresholds)),
        evidence_threshold = evidence_thresholds,
        k1 = 1 / evidence_thresholds,
        k0 = evidence_thresholds,
        h1_threshold_id = paste0("h1-bf",
                                 bfpwr_sim_decimal_id(evidence_thresholds)),
        h0_threshold_id = paste0("h0-bf",
                                 bfpwr_sim_decimal_id(evidence_thresholds)),
        stringsAsFactors = FALSE
    )
}

bfpwr_sim_search_targets <- function(evidence_thresholds = c(3, 10, 30),
                                     target_prob = c(0.1, 0.3, 0.8, 0.9, 0.95),
                                     evidence = c("H1", "H0")) {
    evidence <- match.arg(evidence, c("H1", "H0"), several.ok = TRUE)
    grid <- expand.grid(
        evidence_threshold = sort(unique(as.numeric(evidence_thresholds))),
        target_prob = sort(unique(as.numeric(target_prob))),
        evidence = evidence,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    if (any(!is.finite(grid$evidence_threshold)) ||
        any(grid$evidence_threshold <= 1) ||
        any(!is.finite(grid$target_prob)) ||
        any(grid$target_prob <= 0 | grid$target_prob >= 1)) {
        stop("invalid evidence thresholds or target probabilities")
    }
    if (nrow(grid) == 0) {
        grid$search_target_id <- character(0)
        return(grid[c("search_target_id", "evidence", "evidence_threshold",
                      "target_prob")])
    }
    grid$search_target_id <- paste0(
        tolower(grid$evidence), "-bf",
        bfpwr_sim_decimal_id(grid$evidence_threshold),
        "-p", bfpwr_sim_decimal_id(grid$target_prob)
    )
    grid[c("search_target_id", "evidence", "evidence_threshold",
           "target_prob")]
}

bfpwr_sim_schedule <- function(schedule_id,
                               n = NULL,
                               start_n = NULL,
                               increment = NULL,
                               n_looks = NULL,
                               max_n = NULL,
                               schedule_family_id = NULL,
                               tags = character()) {
    if (is.null(n)) {
        if (is.null(start_n) || is.null(increment)) {
            stop("start_n and increment are required when n is not supplied")
        }
        if (is.null(n_looks) && is.null(max_n)) {
            stop("either n_looks or max_n is required when n is not supplied")
        }
        if (!is.null(n_looks)) {
            n <- start_n + increment * seq.int(0, n_looks - 1L)
            if (!is.null(max_n) && max(n) > max_n) {
                stop("schedule exceeds max_n")
            }
        } else {
            n <- seq(from = start_n, to = max_n, by = increment)
        }
    } else {
        n <- as.numeric(n)
        start_n <- n[[1]]
        n_looks <- length(n)
        max_n <- max(n)
        diffs <- unique(diff(n))
        increment <- if (length(diffs) == 1) diffs[[1]] else NA_real_
    }

    n <- as.integer(n)
    if (length(n) == 0 ||
        any(!is.finite(n)) ||
        any(n <= 0) ||
        any(n != round(n)) ||
        any(diff(n) <= 0)) {
        stop("schedule n must be a positive increasing integer vector")
    }
    if (is.null(schedule_family_id)) {
        schedule_family_id <- if (is.finite(increment)) {
            paste0("start", start_n, "-by", increment)
        } else {
            schedule_id
        }
    }

    list(
        schedule_id = schedule_id,
        schedule_family_id = schedule_family_id,
        n = n,
        start_n = as.integer(start_n),
        increment = if (is.finite(increment)) as.integer(increment) else NA_integer_,
        n_looks = as.integer(length(n)),
        max_n = as.integer(max(n)),
        tags = bfpwr_sim_normalize_tags(tags)
    )
}

bfpwr_sim_validate_schedule <- function(schedule, design) {
    required <- c("schedule_id", "schedule_family_id", "n", "start_n",
                  "increment", "n_looks", "max_n")
    missing <- setdiff(required, names(schedule))
    if (length(missing) > 0) {
        stop("schedule is missing fields: ", paste(missing, collapse = ", "))
    }
    bfpwr_sim_validate_look_grid(schedule$n)
    if (!all(schedule$n %in% design$look_grid)) {
        stop("schedule contains n values not present in design look grid: ",
             schedule$schedule_id)
    }
    invisible(TRUE)
}

bfpwr_sim_schedule_is_compatible <- function(schedule, design) {
    isTRUE(tryCatch({
        bfpwr_sim_validate_schedule(schedule, design)
        TRUE
    }, error = function(e) FALSE))
}

bfpwr_sim_fixed_fixture_spec <- function(fixture_set_id,
                                         family = "z",
                                         bf_types = NULL,
                                         bf_prior_ids = NULL,
                                         evidence_thresholds = c(3, 10, 30),
                                         n = NULL,
                                         search_targets = bfpwr_sim_search_targets(evidence_thresholds),
                                         tags = character(),
                                         rationale = "") {
    if (!family %in% c("z", "t", "binomial")) {
        stop("unsupported fixture family: ", family)
    }
    list(
        fixture_set_id = fixture_set_id,
        family = family,
        mode = "fixed",
        bf_types = bf_types,
        bf_prior_ids = bf_prior_ids,
        evidence_thresholds = sort(unique(as.numeric(evidence_thresholds))),
        tail_thresholds = bfpwr_sim_tail_thresholds(evidence_thresholds),
        threshold_pairs = bfpwr_sim_threshold_pairs(evidence_thresholds),
        n = if (is.null(n)) NULL else sort(unique(as.integer(n))),
        search_targets = search_targets,
        tags = bfpwr_sim_normalize_tags(tags),
        rationale = rationale
    )
}

bfpwr_sim_sequential_fixture_spec <- function(fixture_set_id,
                                              family = "z",
                                              bf_types = NULL,
                                              bf_prior_ids = NULL,
                                              evidence_thresholds = c(3, 10, 30),
                                              schedules,
                                              search_targets = bfpwr_sim_search_targets(evidence_thresholds),
                                              tags = character(),
                                              rationale = "") {
    if (!family %in% c("z", "t", "binomial")) {
        stop("unsupported fixture family: ", family)
    }
    if (!is.list(schedules) || length(schedules) == 0) {
        stop("sequential fixture specs require at least one schedule")
    }
    list(
        fixture_set_id = fixture_set_id,
        family = family,
        mode = "sequential",
        bf_types = bf_types,
        bf_prior_ids = bf_prior_ids,
        evidence_thresholds = sort(unique(as.numeric(evidence_thresholds))),
        threshold_pairs = bfpwr_sim_threshold_pairs(evidence_thresholds),
        schedules = schedules,
        search_targets = search_targets,
        tags = bfpwr_sim_normalize_tags(tags),
        rationale = rationale
    )
}

bfpwr_sim_select_fixture_bf_priors <- function(fixture, bf_priors,
                                               corpus_root = NULL,
                                               available_only = FALSE) {
    keep <- vapply(bf_priors, function(x) x$test_family == fixture$family,
                   logical(1))
    if (!is.null(fixture$bf_types)) {
        keep <- keep & vapply(bf_priors, function(x) x$bf_type %in% fixture$bf_types,
                              logical(1))
    }
    if (!is.null(fixture$bf_prior_ids)) {
        keep <- keep & vapply(bf_priors, function(x) {
            x$bf_prior_id %in% fixture$bf_prior_ids
        }, logical(1))
    }
    selected <- bf_priors[keep]
    if (available_only) {
        if (is.null(corpus_root)) {
            stop("corpus_root is required when available_only = TRUE")
        }
        selected <- selected[vapply(selected, function(x) {
            dir.exists(bfpwr_sim_bf_prior_dir(corpus_root, x))
        }, logical(1))]
    }
    selected
}

bfpwr_sim_find_design_for_bf_prior <- function(bf_prior, designs) {
    ids <- vapply(designs, function(x) x$design_case_id, character(1))
    hit <- match(bf_prior$design_case_id, ids)
    if (is.na(hit)) {
        stop("design case not found for BF prior: ", bf_prior$bf_prior_id)
    }
    designs[[hit]]
}

bfpwr_sim_summarize_fixed_tails <- function(corpus_root,
                                            bf_prior,
                                            design,
                                            thresholds,
                                            n = NULL,
                                            fixture_set_id = NA_character_,
                                            validate = FALSE) {
    thresholds <- if (is.data.frame(thresholds)) {
        thresholds
    } else {
        bfpwr_sim_tail_thresholds(thresholds)
    }
    required_threshold_cols <- c("threshold_id", "evidence_threshold", "tail",
                                 "threshold", "log_threshold")
    missing <- setdiff(required_threshold_cols, names(thresholds))
    if (length(missing) > 0) {
        stop("thresholds are missing columns: ", paste(missing, collapse = ", "))
    }
    if (is.null(n)) {
        n_values <- design$look_grid
    } else {
        n_values <- sort(unique(as.integer(n)))
        missing_n <- setdiff(n_values, design$look_grid)
        if (length(missing_n) > 0) {
            stop("requested n values are not in the design grid: ",
                 paste(missing_n, collapse = ", "))
        }
    }

    chunk_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior, design)
    missing_files <- chunk_files[!file.exists(chunk_files)]
    if (length(missing_files) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing_files,
            asset_set = "chunk-validation",
            label = "BF chunk files")
    }

    n_count <- length(n_values)
    threshold_count <- nrow(thresholds)
    event_counts <- matrix(0, nrow = n_count, ncol = threshold_count)
    total_counts <- rep(0, n_count)
    finite_counts <- rep(0, n_count)
    infinite_counts <- rep(0, n_count)
    nan_counts <- rep(0, n_count)
    na_counts <- rep(0, n_count)

    for (file in chunk_files) {
        bf_rows <- readRDS(file)
        if (validate) {
            bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
        }
        keep <- bf_rows$n %in% n_values
        if (!any(keep)) next
        n_index <- match(bf_rows$n[keep], n_values)
        log_bf01 <- bf_rows$log_bf01[keep]

        total_counts <- total_counts + tabulate(n_index, nbins = n_count)
        finite_counts <- finite_counts +
            tabulate(n_index[is.finite(log_bf01)], nbins = n_count)
        infinite_counts <- infinite_counts +
            tabulate(n_index[is.infinite(log_bf01)], nbins = n_count)
        nan_counts <- nan_counts +
            tabulate(n_index[is.nan(log_bf01)], nbins = n_count)
        na_counts <- na_counts +
            tabulate(n_index[is.na(log_bf01) & !is.nan(log_bf01)],
                     nbins = n_count)

        for (j in seq_len(threshold_count)) {
            hit <- if (thresholds$tail[[j]] == "H1") {
                log_bf01 <= thresholds$log_threshold[[j]]
            } else if (thresholds$tail[[j]] == "H0") {
                log_bf01 >= thresholds$log_threshold[[j]]
            } else {
                stop("unsupported threshold tail: ", thresholds$tail[[j]])
            }
            hit[is.na(hit)] <- FALSE
            event_counts[, j] <- event_counts[, j] +
                tabulate(n_index[hit], nbins = n_count)
        }
    }

    rows <- vector("list", threshold_count)
    for (j in seq_len(threshold_count)) {
        n_event <- event_counts[, j]
        prob <- n_event / total_counts
        rows[[j]] <- data.frame(
            fixture_set_id = fixture_set_id,
            family = bf_prior$test_family,
            bf_type = bf_prior$bf_type,
            bf_prior_id = bf_prior$bf_prior_id,
            design_case_id = design$design_case_id,
            look_grid_name = design$look_grid_name,
            threshold_id = thresholds$threshold_id[[j]],
            evidence_threshold = thresholds$evidence_threshold[[j]],
            tail = thresholds$tail[[j]],
            threshold = thresholds$threshold[[j]],
            log_threshold = thresholds$log_threshold[[j]],
            n = n_values,
            nsim = total_counts,
            n_event = n_event,
            prob = prob,
            mcse = bfpwr_sim_mcse(prob, total_counts),
            n_finite_log_bf01 = finite_counts,
            n_infinite_log_bf01 = infinite_counts,
            n_nan_log_bf01 = nan_counts,
            n_na_log_bf01 = na_counts,
            source_chunk_count = length(chunk_files),
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, rows)
}

bfpwr_sim_summarize_fixed_decisions <- function(tail_summary,
                                                threshold_pairs) {
    if (!is.data.frame(tail_summary) || nrow(tail_summary) == 0) {
        return(data.frame())
    }
    threshold_pairs <- if (is.data.frame(threshold_pairs)) {
        threshold_pairs
    } else {
        bfpwr_sim_threshold_pairs(threshold_pairs)
    }

    rows <- vector("list", nrow(threshold_pairs))
    for (i in seq_len(nrow(threshold_pairs))) {
        pair <- threshold_pairs[i, , drop = FALSE]
        h1 <- tail_summary[
            tail_summary$tail == "H1" &
                tail_summary$evidence_threshold == pair$evidence_threshold,
            , drop = FALSE]
        h0 <- tail_summary[
            tail_summary$tail == "H0" &
                tail_summary$evidence_threshold == pair$evidence_threshold,
            , drop = FALSE]
        key_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                      "design_case_id", "look_grid_name", "n")
        key <- do.call(paste, c(h1[key_cols], sep = "\r"))
        h0_key <- do.call(paste, c(h0[key_cols], sep = "\r"))
        h0 <- h0[match(key, h0_key), , drop = FALSE]
        if (any(is.na(h0$prob))) {
            stop("could not align H0 and H1 tail summaries")
        }
        n_inc <- h1$nsim - h1$n_event - h0$n_event
        if (any(n_inc < 0)) {
            stop("fixed decision counts imply negative inconclusive counts")
        }
        p_inc <- n_inc / h1$nsim
        rows[[i]] <- data.frame(
            h1[key_cols],
            threshold_pair_id = pair$threshold_pair_id,
            evidence_threshold = pair$evidence_threshold,
            k1 = pair$k1,
            k0 = pair$k0,
            pH1 = h1$prob,
            pH0 = h0$prob,
            pInc = p_inc,
            nH1 = h1$n_event,
            nH0 = h0$n_event,
            nInc = n_inc,
            nsim = h1$nsim,
            mcse_pH1 = h1$mcse,
            mcse_pH0 = h0$mcse,
            mcse_pInc = bfpwr_sim_mcse(p_inc, h1$nsim),
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, rows)
}

bfpwr_sim_search_fixed_summary <- function(tail_summary, search_targets) {
    if (!is.data.frame(tail_summary) || nrow(tail_summary) == 0 ||
        !is.data.frame(search_targets) || nrow(search_targets) == 0) {
        return(data.frame())
    }
    rows <- list()
    row_id <- 0L
    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "tail",
                    "evidence_threshold", "threshold", "log_threshold")
    groups <- split(tail_summary, do.call(paste, c(tail_summary[group_cols],
                                                   sep = "\r")))
    for (group in groups) {
        group <- group[order(group$n), , drop = FALSE]
        targets <- search_targets[
            search_targets$evidence == group$tail[[1]] &
                search_targets$evidence_threshold == group$evidence_threshold[[1]],
            , drop = FALSE]
        if (nrow(targets) == 0) next
        for (i in seq_len(nrow(targets))) {
            hit <- which(group$prob >= targets$target_prob[[i]])
            achieved <- length(hit) > 0
            first <- if (achieved) hit[[1]] else NA_integer_
            previous <- if (achieved && first > 1) first - 1L else NA_integer_
            best <- which.max(group$prob)
            row_id <- row_id + 1L
            rows[[row_id]] <- data.frame(
                group[1, group_cols, drop = FALSE],
                search_target_id = targets$search_target_id[[i]],
                target_prob = targets$target_prob[[i]],
                criterion = "first_n_with_prob_ge_target",
                achieved = achieved,
                n_found = if (achieved) group$n[[first]] else NA_integer_,
                prob_found = if (achieved) group$prob[[first]] else NA_real_,
                mcse_found = if (achieved) group$mcse[[first]] else NA_real_,
                n_previous = if (is.na(previous)) NA_integer_ else group$n[[previous]],
                prob_previous = if (is.na(previous)) NA_real_ else group$prob[[previous]],
                best_n = group$n[[best]],
                best_prob = group$prob[[best]],
                max_n_searched = max(group$n),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

bfpwr_sim_bf_matrix_for_schedule <- function(bf_rows, schedule) {
    keep <- bf_rows$n %in% schedule$n
    x <- bf_rows[keep, , drop = FALSE]
    if (nrow(x) == 0) {
        stop("no BF rows match schedule: ", schedule$schedule_id)
    }
    stage <- match(x$n, schedule$n)
    ord <- order(x$replicate_id, stage)
    x <- x[ord, , drop = FALSE]
    stage <- stage[ord]
    replicate_ids <- unique(x$replicate_id)
    n_stage <- length(schedule$n)
    expected_rows <- length(replicate_ids) * n_stage
    if (nrow(x) != expected_rows) {
        stop("schedule rows are incomplete for schedule: ",
             schedule$schedule_id)
    }
    expected_stage <- rep(seq_len(n_stage), times = length(replicate_ids))
    if (!identical(stage, expected_stage)) {
        stop("schedule rows are not complete by replicate for schedule: ",
             schedule$schedule_id)
    }
    list(
        log_bf01 = matrix(x$log_bf01, nrow = length(replicate_ids),
                          ncol = n_stage, byrow = TRUE),
        replicate_ids = replicate_ids,
        n_finite_log_bf01 = sum(is.finite(x$log_bf01)),
        n_infinite_log_bf01 = sum(is.infinite(x$log_bf01)),
        n_nan_log_bf01 = sum(is.nan(x$log_bf01)),
        n_na_log_bf01 = sum(is.na(x$log_bf01) & !is.nan(x$log_bf01))
    )
}

bfpwr_sim_sequential_counts <- function(log_bf01, threshold_pair) {
    log_k1 <- log(threshold_pair$k1[[1]])
    log_k0 <- log(threshold_pair$k0[[1]])
    hit_h1 <- log_bf01 <= log_k1
    hit_h0 <- log_bf01 >= log_k0
    hit_h1[is.na(hit_h1)] <- FALSE
    hit_h0[is.na(hit_h0)] <- FALSE
    hit_any <- hit_h1 | hit_h0

    first <- apply(hit_any, 1L, function(x) {
        hit <- which(x)
        if (length(hit) == 0) NA_integer_ else hit[[1]]
    })
    n_stage <- ncol(log_bf01)
    stop_h1 <- rep(0, n_stage)
    stop_h0 <- rep(0, n_stage)
    hit_rows <- which(!is.na(first))
    if (length(hit_rows) > 0) {
        for (i in hit_rows) {
            j <- first[[i]]
            if (hit_h0[i, j]) {
                stop_h0[[j]] <- stop_h0[[j]] + 1L
            } else {
                stop_h1[[j]] <- stop_h1[[j]] + 1L
            }
        }
    }
    list(stop_h1 = stop_h1, stop_h0 = stop_h0, total = nrow(log_bf01))
}

bfpwr_sim_summarize_sequential <- function(corpus_root,
                                           bf_prior,
                                           design,
                                           schedules,
                                           threshold_pairs,
                                           fixture_set_id = NA_character_,
                                           validate = FALSE) {
    compatible <- vapply(schedules, bfpwr_sim_schedule_is_compatible,
                         logical(1), design = design)
    schedules <- schedules[compatible]
    if (length(schedules) == 0) {
        return(list(cumulative_summary = data.frame(),
                    final_summary = data.frame()))
    }
    threshold_pairs <- if (is.data.frame(threshold_pairs)) {
        threshold_pairs
    } else {
        bfpwr_sim_threshold_pairs(threshold_pairs)
    }

    chunk_files <- bfpwr_sim_bf_chunk_files(corpus_root, bf_prior, design)
    missing_files <- chunk_files[!file.exists(chunk_files)]
    if (length(missing_files) > 0) {
        bfpwr_sim_stop_missing_corpus_files(
            corpus_root,
            missing = missing_files,
            asset_set = "chunk-validation",
            label = "BF chunk files")
    }

    accum <- list()
    for (s in schedules) {
        for (p in seq_len(nrow(threshold_pairs))) {
            key <- paste(s$schedule_id,
                         threshold_pairs$threshold_pair_id[[p]], sep = "\r")
            accum[[key]] <- list(
                schedule = s,
                threshold_pair = threshold_pairs[p, , drop = FALSE],
                stop_h1 = rep(0, length(s$n)),
                stop_h0 = rep(0, length(s$n)),
                total = 0,
                n_finite_log_bf01 = 0,
                n_infinite_log_bf01 = 0,
                n_nan_log_bf01 = 0,
                n_na_log_bf01 = 0
            )
        }
    }

    for (file in chunk_files) {
        bf_rows <- readRDS(file)
        if (validate) {
            bfpwr_sim_validate_bf_rows(bf_prior, design, bf_rows)
        }
        matrices <- lapply(schedules, function(s) {
            bfpwr_sim_bf_matrix_for_schedule(bf_rows, s)
        })
        names(matrices) <- vapply(schedules, function(s) s$schedule_id,
                                  character(1))
        for (s in schedules) {
            matrix_info <- matrices[[s$schedule_id]]
            for (p in seq_len(nrow(threshold_pairs))) {
                pair <- threshold_pairs[p, , drop = FALSE]
                key <- paste(s$schedule_id, pair$threshold_pair_id[[1]],
                             sep = "\r")
                counts <- bfpwr_sim_sequential_counts(matrix_info$log_bf01,
                                                       pair)
                accum[[key]]$stop_h1 <- accum[[key]]$stop_h1 + counts$stop_h1
                accum[[key]]$stop_h0 <- accum[[key]]$stop_h0 + counts$stop_h0
                accum[[key]]$total <- accum[[key]]$total + counts$total
                accum[[key]]$n_finite_log_bf01 <-
                    accum[[key]]$n_finite_log_bf01 +
                    matrix_info$n_finite_log_bf01
                accum[[key]]$n_infinite_log_bf01 <-
                    accum[[key]]$n_infinite_log_bf01 +
                    matrix_info$n_infinite_log_bf01
                accum[[key]]$n_nan_log_bf01 <-
                    accum[[key]]$n_nan_log_bf01 +
                    matrix_info$n_nan_log_bf01
                accum[[key]]$n_na_log_bf01 <-
                    accum[[key]]$n_na_log_bf01 +
                    matrix_info$n_na_log_bf01
            }
        }
    }

    cumulative_rows <- list()
    final_rows <- list()
    for (key in names(accum)) {
        a <- accum[[key]]
        s <- a$schedule
        pair <- a$threshold_pair
        total <- a$total
        cum_h1 <- cumsum(a$stop_h1)
        cum_h0 <- cumsum(a$stop_h0)
        cum_pH1 <- cum_h1 / total
        cum_pH0 <- cum_h0 / total
        cum_inc <- total - cum_h1 - cum_h0
        if (any(cum_inc < 0)) {
            stop("sequential counts imply negative inconclusive counts")
        }
        cum_pInc <- cum_inc / total
        stopped_by_look <- cum_h1 + cum_h0
        EN_to_look <- vapply(seq_along(s$n), function(j) {
            stopped_n <- sum((a$stop_h1[seq_len(j)] +
                              a$stop_h0[seq_len(j)]) * s$n[seq_len(j)])
            active_n <- (total - stopped_by_look[[j]]) * s$n[[j]]
            (stopped_n + active_n) / total
        }, numeric(1))

        cumulative_rows[[length(cumulative_rows) + 1L]] <- data.frame(
            fixture_set_id = fixture_set_id,
            family = bf_prior$test_family,
            bf_type = bf_prior$bf_type,
            bf_prior_id = bf_prior$bf_prior_id,
            design_case_id = design$design_case_id,
            look_grid_name = design$look_grid_name,
            schedule_id = s$schedule_id,
            schedule_family_id = s$schedule_family_id,
            start_n = s$start_n,
            increment = s$increment,
            n_looks = s$n_looks,
            max_n = s$max_n,
            threshold_pair_id = pair$threshold_pair_id,
            evidence_threshold = pair$evidence_threshold,
            k1 = pair$k1,
            k0 = pair$k0,
            look = seq_along(s$n),
            n = s$n,
            nsim = total,
            stop_H1_at_look = a$stop_h1,
            stop_H0_at_look = a$stop_h0,
            cum_pH1 = cum_pH1,
            cum_pH0 = cum_pH0,
            cum_pInc = cum_pInc,
            mcse_cum_pH1 = bfpwr_sim_mcse(cum_pH1, total),
            mcse_cum_pH0 = bfpwr_sim_mcse(cum_pH0, total),
            mcse_cum_pInc = bfpwr_sim_mcse(cum_pInc, total),
            EN_to_look = EN_to_look,
            n_finite_log_bf01 = a$n_finite_log_bf01,
            n_infinite_log_bf01 = a$n_infinite_log_bf01,
            n_nan_log_bf01 = a$n_nan_log_bf01,
            n_na_log_bf01 = a$n_na_log_bf01,
            source_chunk_count = length(chunk_files),
            stringsAsFactors = FALSE
        )

        stop_counts <- a$stop_h1 + a$stop_h0
        n_inc <- total - sum(stop_counts)
        stop_n_sum <- sum(stop_counts * s$n) + n_inc * s$max_n
        stop_n2_sum <- sum(stop_counts * s$n^2) + n_inc * s$max_n^2
        EN <- stop_n_sum / total
        VarN <- stop_n2_sum / total - EN^2
        pH1 <- sum(a$stop_h1) / total
        pH0 <- sum(a$stop_h0) / total
        pInc <- n_inc / total

        stop_n_values <- rep(s$max_n, n_inc)
        if (sum(stop_counts) > 0) {
            stop_n_values <- c(rep(s$n, stop_counts), stop_n_values)
        }
        final_rows[[length(final_rows) + 1L]] <- data.frame(
            fixture_set_id = fixture_set_id,
            family = bf_prior$test_family,
            bf_type = bf_prior$bf_type,
            bf_prior_id = bf_prior$bf_prior_id,
            design_case_id = design$design_case_id,
            look_grid_name = design$look_grid_name,
            schedule_id = s$schedule_id,
            schedule_family_id = s$schedule_family_id,
            start_n = s$start_n,
            increment = s$increment,
            n_looks = s$n_looks,
            max_n = s$max_n,
            threshold_pair_id = pair$threshold_pair_id,
            evidence_threshold = pair$evidence_threshold,
            k1 = pair$k1,
            k0 = pair$k0,
            nsim = total,
            pH1 = pH1,
            pH0 = pH0,
            pInc = pInc,
            nH1 = sum(a$stop_h1),
            nH0 = sum(a$stop_h0),
            nInc = n_inc,
            mcse_pH1 = bfpwr_sim_mcse(pH1, total),
            mcse_pH0 = bfpwr_sim_mcse(pH0, total),
            mcse_pInc = bfpwr_sim_mcse(pInc, total),
            EN = EN,
            VarN = VarN,
            q25_N = unname(stats::quantile(stop_n_values, 0.25,
                                           type = 1)),
            median_N = unname(stats::quantile(stop_n_values, 0.5,
                                              type = 1)),
            q75_N = unname(stats::quantile(stop_n_values, 0.75,
                                           type = 1)),
            p_stop_by_final = 1 - pInc,
            n_finite_log_bf01 = a$n_finite_log_bf01,
            n_infinite_log_bf01 = a$n_infinite_log_bf01,
            n_nan_log_bf01 = a$n_nan_log_bf01,
            n_na_log_bf01 = a$n_na_log_bf01,
            source_chunk_count = length(chunk_files),
            stringsAsFactors = FALSE
        )
    }

    list(
        cumulative_summary = if (length(cumulative_rows) == 0) {
            data.frame()
        } else {
            do.call(rbind, cumulative_rows)
        },
        final_summary = if (length(final_rows) == 0) {
            data.frame()
        } else {
            do.call(rbind, final_rows)
        }
    )
}

bfpwr_sim_search_sequential_summary <- function(cumulative_summary,
                                                search_targets) {
    empty_summary <- function() {
        data.frame(
            fixture_set_id = character(),
            family = character(),
            bf_type = character(),
            bf_prior_id = character(),
            design_case_id = character(),
            look_grid_name = character(),
            schedule_id = character(),
            schedule_family_id = character(),
            start_n = integer(),
            increment = integer(),
            n_looks = integer(),
            max_n = integer(),
            threshold_pair_id = character(),
            evidence_threshold = numeric(),
            k1 = numeric(),
            k0 = numeric(),
            search_target_id = character(),
            evidence = character(),
            target_prob = numeric(),
            criterion = character(),
            achieved = logical(),
            look_found = integer(),
            n_found = integer(),
            prob_found = numeric(),
            mcse_found = numeric(),
            look_previous = integer(),
            n_previous = integer(),
            prob_previous = numeric(),
            best_look = integer(),
            best_n = integer(),
            best_prob = numeric(),
            max_n_searched = integer(),
            stringsAsFactors = FALSE
        )
    }
    if (!is.data.frame(cumulative_summary) || nrow(cumulative_summary) == 0 ||
        !is.data.frame(search_targets) || nrow(search_targets) == 0) {
        return(empty_summary())
    }
    rows <- list()
    row_id <- 0L
    group_cols <- c("fixture_set_id", "family", "bf_type", "bf_prior_id",
                    "design_case_id", "look_grid_name", "schedule_id",
                    "schedule_family_id", "start_n", "increment", "n_looks",
                    "max_n", "threshold_pair_id", "evidence_threshold",
                    "k1", "k0")
    groups <- split(cumulative_summary,
                    do.call(paste, c(cumulative_summary[group_cols],
                                     sep = "\r")))
    for (group in groups) {
        group <- group[order(group$look), , drop = FALSE]
        targets <- search_targets[
            search_targets$evidence_threshold == group$evidence_threshold[[1]],
            , drop = FALSE]
        if (nrow(targets) == 0) next
        for (i in seq_len(nrow(targets))) {
            prob_col <- if (targets$evidence[[i]] == "H1") "cum_pH1" else "cum_pH0"
            mcse_col <- if (targets$evidence[[i]] == "H1") {
                "mcse_cum_pH1"
            } else {
                "mcse_cum_pH0"
            }
            prob <- group[[prob_col]]
            hit <- which(prob >= targets$target_prob[[i]])
            achieved <- length(hit) > 0
            first <- if (achieved) hit[[1]] else NA_integer_
            previous <- if (achieved && first > 1) first - 1L else NA_integer_
            best <- which.max(prob)
            row_id <- row_id + 1L
            rows[[row_id]] <- data.frame(
                group[1, group_cols, drop = FALSE],
                search_target_id = targets$search_target_id[[i]],
                evidence = targets$evidence[[i]],
                target_prob = targets$target_prob[[i]],
                criterion = "first_look_with_cumulative_prob_ge_target",
                achieved = achieved,
                look_found = if (achieved) group$look[[first]] else NA_integer_,
                n_found = if (achieved) group$n[[first]] else NA_integer_,
                prob_found = if (achieved) prob[[first]] else NA_real_,
                mcse_found = if (achieved) group[[mcse_col]][[first]] else NA_real_,
                look_previous = if (is.na(previous)) NA_integer_ else group$look[[previous]],
                n_previous = if (is.na(previous)) NA_integer_ else group$n[[previous]],
                prob_previous = if (is.na(previous)) NA_real_ else prob[[previous]],
                best_look = group$look[[best]],
                best_n = group$n[[best]],
                best_prob = prob[[best]],
                max_n_searched = max(group$n),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows) == 0) empty_summary() else do.call(rbind, rows)
}

bfpwr_sim_write_fixture_manifest <- function(fixture_dir, files) {
    files <- normalizePath(files[file.exists(files)], winslash = "/",
                           mustWork = TRUE)
    rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1",
                                normalizePath(fixture_dir, winslash = "/")),
               "/?"), "", files)
    manifest <- data.frame(
        file = rel,
        bytes = file.info(files)$size,
        sha256 = unname(tools::sha256sum(files)),
        stringsAsFactors = FALSE
    )
    utils::write.csv(manifest, file.path(fixture_dir, "manifest.csv"),
                     row.names = FALSE)
    bfpwr_sim_write_hashes(c(files, file.path(fixture_dir, "manifest.csv")),
                           file.path(fixture_dir, "sha256.txt"))
    manifest
}

bfpwr_sim_materialize_fixed_fixture <- function(corpus_root,
                                                fixture,
                                                bf_priors = bfpwr_sim_bf_prior_case_set("production"),
                                                designs = bfpwr_sim_design_case_set("production"),
                                                available_only = TRUE,
                                                validate = FALSE) {
    if (!identical(fixture$mode, "fixed")) {
        stop("fixture is not a fixed fixture spec")
    }
    selected <- bfpwr_sim_select_fixture_bf_priors(
        fixture, bf_priors, corpus_root = corpus_root,
        available_only = available_only)
    if (length(selected) == 0) {
        stop("no BF prior cases selected for fixture: ",
             fixture$fixture_set_id)
    }

    tail_rows <- vector("list", length(selected))
    for (i in seq_along(selected)) {
        bf_prior <- selected[[i]]
        design <- bfpwr_sim_find_design_for_bf_prior(bf_prior, designs)
        tail_rows[[i]] <- bfpwr_sim_summarize_fixed_tails(
            corpus_root = corpus_root,
            bf_prior = bf_prior,
            design = design,
            thresholds = fixture$tail_thresholds,
            n = fixture$n,
            fixture_set_id = fixture$fixture_set_id,
            validate = validate)
    }
    tail_summary <- do.call(rbind, tail_rows)
    decision_summary <- bfpwr_sim_summarize_fixed_decisions(
        tail_summary, fixture$threshold_pairs)
    search_summary <- bfpwr_sim_search_fixed_summary(
        tail_summary, fixture$search_targets)

    fixture_dir <- bfpwr_sim_fixture_dir(corpus_root, fixture)
    dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
    files <- c(
        bfpwr_sim_write_rds(fixture, file.path(fixture_dir, "spec.rds")),
        bfpwr_sim_write_rds(tail_summary,
                            file.path(fixture_dir, "tail-summary.rds")),
        bfpwr_sim_write_rds(decision_summary,
                            file.path(fixture_dir, "decision-summary.rds")),
        bfpwr_sim_write_rds(search_summary,
                            file.path(fixture_dir, "search-summary.rds"))
    )
    overview <- data.frame(
        fixture_set_id = fixture$fixture_set_id,
        family = fixture$family,
        mode = fixture$mode,
        bf_priors = length(selected),
        tail_rows = nrow(tail_summary),
        decision_rows = nrow(decision_summary),
        search_rows = nrow(search_summary),
        stringsAsFactors = FALSE
    )
    overview_file <- file.path(fixture_dir, "overview.csv")
    utils::write.csv(overview, overview_file, row.names = FALSE)
    files <- c(files, overview_file)
    manifest <- bfpwr_sim_write_fixture_manifest(fixture_dir, files)
    list(
        fixture_dir = fixture_dir,
        selected_bf_priors = vapply(selected, function(x) x$bf_prior_id,
                                    character(1)),
        overview = overview,
        manifest = manifest,
        tail_summary = tail_summary,
        decision_summary = decision_summary,
        search_summary = search_summary
    )
}

bfpwr_sim_materialize_sequential_fixture <- function(corpus_root,
                                                     fixture,
                                                     bf_priors = bfpwr_sim_bf_prior_case_set("production"),
                                                     designs = bfpwr_sim_design_case_set("production"),
                                                     available_only = TRUE,
                                                     validate = FALSE) {
    if (!identical(fixture$mode, "sequential")) {
        stop("fixture is not a sequential fixture spec")
    }
    selected <- bfpwr_sim_select_fixture_bf_priors(
        fixture, bf_priors, corpus_root = corpus_root,
        available_only = available_only)
    if (length(selected) == 0) {
        stop("no BF prior cases selected for fixture: ",
             fixture$fixture_set_id)
    }

    cumulative_rows <- list()
    final_rows <- list()
    for (i in seq_along(selected)) {
        bf_prior <- selected[[i]]
        design <- bfpwr_sim_find_design_for_bf_prior(bf_prior, designs)
        summary <- bfpwr_sim_summarize_sequential(
            corpus_root = corpus_root,
            bf_prior = bf_prior,
            design = design,
            schedules = fixture$schedules,
            threshold_pairs = fixture$threshold_pairs,
            fixture_set_id = fixture$fixture_set_id,
            validate = validate)
        if (nrow(summary$cumulative_summary) > 0) {
            cumulative_rows[[length(cumulative_rows) + 1L]] <-
                summary$cumulative_summary
        }
        if (nrow(summary$final_summary) > 0) {
            final_rows[[length(final_rows) + 1L]] <- summary$final_summary
        }
    }
    cumulative_summary <- if (length(cumulative_rows) == 0) {
        data.frame()
    } else {
        do.call(rbind, cumulative_rows)
    }
    final_summary <- if (length(final_rows) == 0) {
        data.frame()
    } else {
        do.call(rbind, final_rows)
    }
    search_summary <- bfpwr_sim_search_sequential_summary(
        cumulative_summary, fixture$search_targets)

    fixture_dir <- bfpwr_sim_fixture_dir(corpus_root, fixture)
    dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
    files <- c(
        bfpwr_sim_write_rds(fixture, file.path(fixture_dir, "spec.rds")),
        bfpwr_sim_write_rds(cumulative_summary,
                            file.path(fixture_dir, "cumulative-summary.rds")),
        bfpwr_sim_write_rds(final_summary,
                            file.path(fixture_dir, "final-summary.rds")),
        bfpwr_sim_write_rds(search_summary,
                            file.path(fixture_dir, "search-summary.rds"))
    )
    overview <- data.frame(
        fixture_set_id = fixture$fixture_set_id,
        family = fixture$family,
        mode = fixture$mode,
        bf_priors = length(selected),
        cumulative_rows = nrow(cumulative_summary),
        final_rows = nrow(final_summary),
        search_rows = nrow(search_summary),
        stringsAsFactors = FALSE
    )
    overview_file <- file.path(fixture_dir, "overview.csv")
    utils::write.csv(overview, overview_file, row.names = FALSE)
    files <- c(files, overview_file)
    manifest <- bfpwr_sim_write_fixture_manifest(fixture_dir, files)
    list(
        fixture_dir = fixture_dir,
        selected_bf_priors = vapply(selected, function(x) x$bf_prior_id,
                                    character(1)),
        overview = overview,
        manifest = manifest,
        cumulative_summary = cumulative_summary,
        final_summary = final_summary,
        search_summary = search_summary
    )
}
