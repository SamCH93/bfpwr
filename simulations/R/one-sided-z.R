## One-sided normal-prior validation on the existing simulation trajectories.
## These cases deliberately exclude numerical extremes covered by unit tests.
bfpwr_sim_one_sided_z_cases <- function() {
    priors <- data.frame(
        prior_id = c("narrow", "standard", "diffuse", "positive", "negative"),
        pm = c(0, 0, 0, 0.2, -0.2),
        psd = c(0.1, 1/sqrt(2), 2, 0.3, 0.3))
    designs <- data.frame(
        design_id = c("z-dpoint-0-usd1-short", "z-dpoint-0p2-usd1-short",
                      "z-dpoint-0p5-usd1-short", "z-dpoint-m0p3-usd1-short",
                      "z-dnorm-0p5-s0p1-usd1-short",
                      "z-dnorm-0p2-s0p3-usd1-short",
                      "z-dnorm-0-s0p5-usd1-short"),
        dpm = c(0, 0.2, 0.5, -0.3, 0.5, 0.2, 0),
        dpsd = c(0, 0, 0, 0, 0.1, 0.3, 0.5))
    cases <- merge(merge(priors, designs, by = NULL),
                   data.frame(alternative = c("greater", "less")), by = NULL)
    cases$group <- "core"
    long <- cases[cases$prior_id == "standard" & cases$dpsd == 0 &
                  cases$dpm %in% c(0, 0.2), ]
    long$design_id <- sub("short$", "long", long$design_id)
    long$group <- "long"
    cases <- rbind(cases, long)
    cases$case_id <- paste(cases$prior_id, cases$alternative, cases$design_id,
                           sep = "--")
    cases
}

bfpwr_sim_one_sided_z_thresholds <- function() {
    data.frame(pair = c("3", "10", "30", "H1-10-H0-3", "H1-3-H0-10"),
               k1 = c(1/3, 0.1, 1/30, 0.1, 1/3),
               k0 = c(3, 10, 30, 3, 10))
}

bfpwr_sim_one_sided_z_read <- function(corpus_root, case) {
    directory <- file.path(corpus_root, "designs", "z", case$design_id)
    spec <- readRDS(file.path(directory, "design.rds"))
    index <- readRDS(file.path(directory, "chunk-index.rds"))
    n <- if (case$group == "long") seq(10, 2000, 10) else seq(10, 500, 10)
    stopifnot(spec$nsim == 10000, spec$generation$usd == 1,
              spec$design_prior$mean == case$dpm,
              spec$design_prior$sd == case$dpsd, all(n %in% spec$look_grid))
    chunks <- lapply(index$file, function(file) {
        x <- readRDS(file.path(directory, file))
        x[x$n %in% n, c("replicate_id", "n", "estimate", "se", "true_effect")]
    })
    x <- do.call(rbind, chunks)
    x <- x[order(x$replicate_id, x$n), ]
    ids <- sort(unique(x$replicate_id))
    stopifnot(length(ids) == 10000, nrow(x) == 10000*length(n),
              identical(x$n, rep(n, 10000)),
              identical(x$replicate_id, rep(ids, each = length(n))),
              all(abs(x$se - 1/sqrt(x$n)) < 1e-12),
              all(is.finite(x$estimate)))
    truth <- matrix(x$true_effect, nrow = 10000, byrow = TRUE)
    stopifnot(all(truth == truth[, 1]))
    list(n = n, replicate_id = ids, true_effect = truth[, 1],
         estimate = matrix(x$estimate, nrow = 10000, byrow = TRUE),
         files = file.path(directory, index$file), sha256 = index$sha256)
}

## Table 1's prior/posterior mass correction, evaluated independently of
## bf01() and zcrit(). The approved priors have moderate truncation masses.
bfpwr_sim_one_sided_z_logbf <- function(estimate, se, pm, psd, alternative) {
    direction <- if (alternative == "greater") 1 else -1
    post_sd <- 1/sqrt(1/se^2 + 1/psd^2)
    post_mean <- (estimate/se^2 + pm/psd^2)*post_sd^2
    stats::dnorm(estimate, 0, se, log = TRUE) -
        stats::dnorm(estimate, pm, sqrt(se^2 + psd^2), log = TRUE) +
        stats::pnorm(direction*pm/psd, log.p = TRUE) -
        stats::pnorm(direction*post_mean/post_sd, log.p = TRUE)
}

## Record first stopping times directly from simulated BFs. Sampling continues
## in the stored data, but a stopped replicate cannot contribute another event.
bfpwr_sim_one_sided_z_stopping <- function(logbf, n, k1, k0) {
    active <- rep(TRUE, nrow(logbf))
    stop_n <- rep(tail(n, 1), nrow(logbf))
    counts <- matrix(0L, nrow = length(n), ncol = 3,
                      dimnames = list(NULL, c("H1", "H0", "Inc")))
    h1_total <- h0_total <- 0L
    for (i in seq_along(n)) {
        h1 <- active & logbf[, i] <= log(k1)
        h0 <- active & logbf[, i] >= log(k0)
        stop_n[h1 | h0] <- n[i]
        active[h1 | h0] <- FALSE
        h1_total <- h1_total + sum(h1)
        h0_total <- h0_total + sum(h0)
        counts[i, ] <- c(h1_total, h0_total, sum(active))
    }
    list(counts = counts, EN = mean(stop_n), VarN = stats::var(stop_n),
         mcse_EN = stats::sd(stop_n)/sqrt(length(stop_n)),
         mcse_VarN = stats::sd((stop_n - mean(stop_n))^2)/sqrt(length(stop_n)))
}

bfpwr_sim_one_sided_z_probability_rows <- function(case, mode, schedule, pair,
                                                   n, counts, prediction) {
    stopifnot(all(is.finite(prediction)), all(prediction >= -1e-10),
              all(prediction <= 1 + 1e-10), all(rowSums(counts) == 10000))
    prediction <- pmin(1, pmax(0, prediction))
    rows <- data.frame(case_id = case$case_id, mode = mode, schedule = schedule,
                       pair = pair$pair, n = rep(n, 3),
                       outcome = rep(c("H1", "H0", "Inc"), each = length(n)),
                       nsim = 10000L, n_event = as.vector(counts),
                       reference_prob = as.vector(prediction))
    rows$prob <- rows$n_event/rows$nsim
    rows$mcse <- bfpwr_sim_mcse(rows$prob, rows$nsim)
    rows$error <- rows$reference_prob - rows$prob
    rows
}

bfpwr_sim_one_sided_z_compare <- function(case, logbf, n) {
    thresholds <- bfpwr_sim_one_sided_z_thresholds()
    probabilities <- moments <- timings <- list()
    if (case$group == "core") {
        for (i in 1:3) {
            pair <- thresholds[i, ]
            t0 <- proc.time()[["elapsed"]]
            h1 <- pbf01(k = pair$k1, n = n, usd = 1, pm = case$pm,
                         psd = case$psd, dpm = case$dpm, dpsd = case$dpsd,
                         alternative = case$alternative)
            h0 <- pbf01(k = pair$k0, n = n, usd = 1, pm = case$pm,
                         psd = case$psd, dpm = case$dpm, dpsd = case$dpsd,
                         alternative = case$alternative, lower.tail = FALSE)
            timings[[length(timings) + 1L]] <- data.frame(
                case_id = case$case_id, operation = "fixed", schedule = "fixed",
                pair = pair$pair, seconds = proc.time()[["elapsed"]] - t0)
            counts <- cbind(colSums(logbf <= log(pair$k1)),
                             colSums(logbf >= log(pair$k0)),
                             colSums(logbf > log(pair$k1) & logbf < log(pair$k0)))
            probabilities[[length(probabilities) + 1L]] <-
                bfpwr_sim_one_sided_z_probability_rows(
                    case, "fixed", "fixed", pair, n, counts,
                    cbind(h1, h0, 1 - h1 - h0))
        }
        schedules <- list(single = 300, equal3 = c(100, 200, 300),
                          uneven3 = c(20, 100, 300), dense50 = n)
    } else {
        schedules <- list(dense100 = n[1:100], dense200 = n)
        thresholds <- thresholds[thresholds$pair == "10", ]
    }
    for (schedule in names(schedules)) {
        looks <- schedules[[schedule]]
        for (i in seq_len(nrow(thresholds))) {
            pair <- thresholds[i, ]
            sim <- bfpwr_sim_one_sided_z_stopping(
                logbf[, match(looks, n), drop = FALSE], looks, pair$k1, pair$k0)
            t0 <- proc.time()[["elapsed"]]
            result <- pbf01seq(k1 = pair$k1, k0 = pair$k0, n = looks,
                               se = 1/sqrt(looks), pm = case$pm, psd = case$psd,
                               dpm = case$dpm, dpsd = case$dpsd,
                               alternative = case$alternative)
            timings[[length(timings) + 1L]] <- data.frame(
                case_id = case$case_id, operation = "sequential",
                schedule = schedule, pair = pair$pair,
                seconds = proc.time()[["elapsed"]] - t0)
            probabilities[[length(probabilities) + 1L]] <-
                bfpwr_sim_one_sided_z_probability_rows(
                    case, "sequential", schedule, pair, looks, sim$counts,
                    cbind(result$cumpH1, result$cumpH0, result$cumpInc))
            moments[[length(moments) + 1L]] <- data.frame(
                case_id = case$case_id, schedule = schedule, pair = pair$pair,
                n = tail(looks, 1), EN = sim$EN, predicted_EN = result$EN,
                mcse_EN = sim$mcse_EN, VarN = sim$VarN,
                predicted_VarN = result$VarN, mcse_VarN = sim$mcse_VarN)
        }
    }
    list(probabilities = do.call(rbind, probabilities),
         moments = do.call(rbind, moments), timings = do.call(rbind, timings))
}

## Fixed searches cover the whole core grid. Sequential searches use selected
## priors and all generating conditions on the stored by-10 schedule. Validate
## the integer result separately from the simulation's coarser n grid.
bfpwr_sim_one_sided_z_search <- function(cases, probabilities) {
    rows <- list()
    for (i in which(cases$group == "core")) {
        case <- cases[i, ]
        for (mode in c("fixed", "sequential")) {
            if (mode == "sequential" &&
                !case$prior_id %in% c("standard", "positive")) next
            pairs <- if (mode == "fixed") c("3", "10", "30") else "10"
            for (pair_id in pairs) for (outcome in c("H1", "H0")) {
                pair <- bfpwr_sim_one_sided_z_thresholds()
                pair <- pair[pair$pair == pair_id, ]
                curve <- probabilities[probabilities$case_id == case$case_id &
                    probabilities$mode == mode & probabilities$pair == pair_id &
                    probabilities$outcome == outcome &
                    probabilities$schedule == if (mode == "fixed") "fixed" else "dense50", ]
                curve <- curve[order(curve$n), ]
                for (target in c(0.8, 0.9)) {
                    warnings <- character()
                    t0 <- proc.time()[["elapsed"]]
                    n_found <- withCallingHandlers({
                        if (mode == "fixed") {
                            nbf01(k = if (outcome == "H1") pair$k1 else pair$k0,
                                  power = target, usd = 1, pm = case$pm,
                                  psd = case$psd, dpm = case$dpm, dpsd = case$dpsd,
                                  alternative = case$alternative,
                                  lower.tail = outcome == "H1", nrange = c(10, 500))
                        } else {
                            nbf01seq(k1 = pair$k1, k0 = pair$k0, power = target,
                                     usd = 1, pm = case$pm, psd = case$psd,
                                     dpm = case$dpm, dpsd = case$dpsd,
                                     alternative = case$alternative, target = outcome,
                                     minN = 10, by = 10, nrange = c(10, 500))
                        }
                    }, warning = function(w) {
                        warnings <<- c(warnings, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    })
                    elapsed <- proc.time()[["elapsed"]] - t0
                    reached <- is.finite(n_found)
                    expected <- which(curve$reference_prob >= target)
                    sim_hit <- which(curve$prob >= target)
                    grid_n <- if (reached) ceiling(n_found/10)*10 else NA_real_
                    hit <- if (reached) match(grid_n, curve$n) else NA_integer_
                    if (mode == "fixed" && reached) {
                        check_n <- c(n_found, max(10, n_found - 1))
                        check_power <- pbf01(
                            k = if (outcome == "H1") pair$k1 else pair$k0,
                            n = check_n, usd = 1, pm = case$pm, psd = case$psd,
                            dpm = case$dpm, dpsd = case$dpsd,
                            alternative = case$alternative, lower.tail = outcome == "H1")
                        valid <- check_power[1] >= target - 1e-7 &&
                            (n_found == 10 || check_power[2] < target + 1e-7)
                    } else {
                        valid <- if (length(expected)) reached &&
                            grid_n == curve$n[expected[1]] else !reached
                    }
                    ## Fixed searchN() warns and returns NaN when the target is
                    ## already exceeded at the lower bound. This is expected
                    ## range handling, distinct from a missed interior crossing.
                    below_range <- mode == "fixed" && !reached &&
                        curve$reference_prob[1] > target
                    if (below_range) valid <- TRUE
                    status <- if (reached) "returned" else if (below_range)
                        "below_range" else if (length(expected))
                        "missed_crossing" else "not_reached"
                    rows[[length(rows) + 1L]] <- data.frame(
                        case_id = case$case_id, mode = mode, pair = pair_id,
                        outcome = outcome, target = target, n_found = n_found,
                        reached = reached, search_valid = valid,
                        status = status,
                        predicted_grid_n = if (length(expected)) curve$n[expected[1]] else NA_real_,
                        grid_n = grid_n,
                        simulation_n = if (length(sim_hit)) curve$n[sim_hit[1]] else NA_real_,
                        simulation_power = if (reached) curve$prob[hit] else NA_real_,
                        simulation_mcse = if (reached) curve$mcse[hit] else NA_real_,
                        predicted_grid_power = if (reached) curve$reference_prob[hit] else NA_real_,
                        seconds = elapsed, warning = paste(unique(warnings), collapse = " | "))
                }
            }
        }
    }
    do.call(rbind, rows)
}

bfpwr_sim_one_sided_z_diagnostics <- function(probabilities, moments) {
    diagnostics <- bfpwr_sim_mc_reference_diagnostics(
        probabilities, label = "one-sided normal z", group_cols = c("mode", "outcome"))
    summaries <- lapply(split(probabilities, probabilities$mode), function(x) {
        data.frame(mode = x$mode[1], comparisons = nrow(x),
                   median_abs_error = stats::median(abs(x$error)),
                   max_abs_error = max(abs(x$error)),
                   within_mc_tolerance = sum(bfpwr_sim_mc_close(
                       x$prob, x$reference_prob, x$nsim)))
    })
    moments$EN_error <- moments$predicted_EN - moments$EN
    moments$VarN_error <- moments$predicted_VarN - moments$VarN
    ## Moment tolerances combine sampling error with a small allowance for
    ## deterministic integration, scaled to the range of possible N.
    moments$EN_ok <- abs(moments$EN_error) <= pmax(5*moments$mcse_EN, 0.002*moments$n)
    moments$VarN_ok <- abs(moments$VarN_error) <=
        pmax(5*moments$mcse_VarN, 0.002*moments$n^2)
    list(probability_summary = do.call(rbind, summaries),
         mc_diagnostics = diagnostics, moments = moments)
}

## Follow up representative discrepancies without changing the default results
## or drawing more simulated data. Record alternative-backend errors explicitly.
bfpwr_sim_one_sided_z_convergence <- function(cases, probabilities) {
    selected <- data.frame(
        case_id = c("narrow--greater--z-dpoint-0-usd1-short",
                    "negative--greater--z-dpoint-0p2-usd1-short",
                    "narrow--less--z-dnorm-0-s0p5-usd1-short",
                    "standard--less--z-dnorm-0-s0p5-usd1-short"),
        n = c(310, 500, 260, 330), pair = c("10", "30", "H1-10-H0-3", "30"),
        outcome = c("H0", "H1", "Inc", "Inc"))
    settings <- data.frame(method = c("lpmvnorm", "lpmvnorm", "pmvnorm"),
                            ngrid = .bfseq_integration_settings(list())$ngrid*c(1L, 10L, 1L))
    rows <- list()
    for (i in seq_len(nrow(selected))) {
        x <- selected[i, ]
        case <- cases[cases$case_id == x$case_id, ]
        pair <- bfpwr_sim_one_sided_z_thresholds()
        pair <- pair[pair$pair == x$pair, ]
        reference <- probabilities[probabilities$case_id == x$case_id &
            probabilities$schedule == "dense50" & probabilities$pair == x$pair &
            probabilities$n == x$n & probabilities$outcome == x$outcome, ]
        stopifnot(nrow(case) == 1, nrow(reference) == 1)
        for (j in seq_len(nrow(settings))) {
            setting <- settings[j, ]
            warning <- error <- ""
            seconds <- NA_real_
            predicted <- reference$reference_prob
            if (j > 1) {
                looks <- seq(10, x$n, 10)
                t0 <- proc.time()[["elapsed"]]
                result <- tryCatch(withCallingHandlers(
                    pbf01seq(k1 = pair$k1, k0 = pair$k0, n = looks,
                        se = 1/sqrt(looks), pm = case$pm, psd = case$psd,
                        dpm = case$dpm, dpsd = case$dpsd,
                        alternative = case$alternative,
                        method = setting$method, ngrid = setting$ngrid),
                    warning = function(w) {
                        warning <<- paste(warning, conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }), error = function(e) {
                        error <<- conditionMessage(e)
                        NULL
                    })
                seconds <- proc.time()[["elapsed"]] - t0
                predicted <- if (is.null(result)) NA_real_ else
                    tail(result[[paste0("cump", x$outcome)]], 1)
            }
            rows[[length(rows) + 1L]] <- data.frame(
                x, setting, predicted = predicted, simulated = reference$prob,
                mcse = reference$mcse, seconds = seconds,
                warning = trimws(warning), error = error)
        }
    }
    do.call(rbind, rows)
}

bfpwr_sim_one_sided_z_plot <- function(bundle, kind = c("curves", "agreement")) {
    kind <- match.arg(kind)
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old))
    colors <- c(H1 = "#2374AB", H0 = "#C65335", Inc = "#777777")
    if (kind == "curves") {
        graphics::par(mfrow = c(3, 2), mar = c(3.5, 4, 2.7, 0.7),
                      oma = c(0, 0, 4, 0), mgp = c(2.2, 0.7, 0), las = 1)
        design_ids <- c("z-dpoint-0-usd1-short", "z-dpoint-0p2-usd1-short",
                        "z-dnorm-0p2-s0p3-usd1-short")
        labels <- c("True effect = 0", "True effect = 0.2",
                    "Effect mean = 0.2, SD = 0.3")
        for (i in seq_along(design_ids)) for (alternative in c("greater", "less")) {
            id <- paste("standard", alternative, design_ids[i], sep = "--")
            x <- bundle$probabilities[bundle$probabilities$case_id == id &
                bundle$probabilities$schedule == "dense50" &
                bundle$probabilities$pair == "10", ]
            graphics::plot(NA, xlim = c(10, 500), ylim = c(0, 100),
                xlab = if (i == 3) "Sample size" else "",
                ylab = "Cumulative probability (%)",
                main = paste0(labels[i], "\nH1: effect ",
                              if (alternative == "greater") "> 0" else "< 0"),
                cex.main = 0.95)
            graphics::grid(col = "#EEEEEE")
            for (outcome in c("Inc", "H0", "H1")) {
                y <- x[x$outcome == outcome, ]
                graphics::polygon(c(y$n, rev(y$n)),
                    100*c(pmax(0, y$prob - 1.96*y$mcse),
                          rev(pmin(1, y$prob + 1.96*y$mcse))),
                    col = grDevices::adjustcolor(colors[outcome], alpha.f = 0.12),
                    border = NA)
                graphics::lines(y$n, 100*y$reference_prob, col = colors[outcome], lwd = 1.8)
                selected <- seq(1, nrow(y), by = 3)
                graphics::points(y$n[selected], 100*y$prob[selected],
                                  col = colors[outcome], pch = 16, cex = 0.45)
            }
        }
        graphics::mtext("One-sided normal z-tests: package predictions and simulations",
                         side = 3, outer = TRUE, line = 2.6, font = 2, cex = 1.1)
        graphics::mtext("10,000 original replicates | Half-normal scale 0.707 | BF thresholds 1/10 and 10",
                         side = 3, outer = TRUE, line = 1.3, cex = 0.82)
        graphics::mtext("Blue: H1   Red: H0   Grey: inconclusive | Lines: package   Points: simulation   Bands: 95% MC uncertainty",
                         side = 3, outer = TRUE, line = 0.1, cex = 0.73)
    } else {
        graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.4, 0.7),
                      mgp = c(2.2, 0.7, 0), las = 1)
        for (mode in c("fixed", "sequential")) for (outcome in c("H1", "H0")) {
            x <- bundle$probabilities[bundle$probabilities$mode == mode &
                                      bundle$probabilities$outcome == outcome, ]
            graphics::plot(100*x$reference_prob, 100*x$prob,
                xlab = "Package probability (%)", ylab = "Simulated probability (%)",
                main = paste(if (mode == "fixed") "Fixed" else "Sequential", outcome),
                xlim = c(0, 100), ylim = c(0, 100), pch = 16, cex = 0.4,
                col = grDevices::adjustcolor(colors[outcome], alpha.f = 0.12))
            graphics::abline(0, 1, lty = 2, col = "#333333")
        }
    }
    invisible(NULL)
}
