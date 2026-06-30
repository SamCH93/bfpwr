## Helper functions for sequential BF sample-size searches
## -----------------------------------------------------------------------------

.bfseq_default_search_lower <- function(nrange, default = 10) {
    lower <- as.integer(ceiling(nrange[1]))
    upper <- as.integer(ceiling(nrange[2]))
    default <- as.integer(ceiling(default))
    if (lower <= default && default <= upper) default else lower
}

.bfseq_schedule_spec <- function(looks = 1, timing = NULL, minN = NULL,
                                 by = NULL, nrange, lookMinN = 2) {
    stopifnot(
        length(lookMinN) == 1,
        is.numeric(lookMinN),
        is.finite(lookMinN),
        lookMinN >= 2
    )
    lookMinN <- as.integer(ceiling(lookMinN))

    if (!is.null(by)) {
        if (!is.null(timing)) {
            stop("only one of 'timing' and 'by' may be specified")
        }
        stopifnot(
            length(by) == 1,
            is.numeric(by),
            is.finite(by),
            by > 0
        )
        if (is.null(minN)) {
            minN <- max(ceiling(nrange[1]), lookMinN)
        }
        stopifnot(
            length(minN) == 1,
            is.numeric(minN),
            is.finite(minN),
            minN >= lookMinN
        )
        return(list(type = "increase", minN = as.integer(ceiling(minN)),
                    by = as.integer(ceiling(by)), lookMinN = lookMinN))
    }

    if (is.null(timing)) {
        stopifnot(
            length(looks) == 1,
            is.numeric(looks),
            is.finite(looks),
            looks >= 1
        )
        looks <- as.integer(ceiling(looks))
        timing <- if (looks == 1) {
            1
        } else {
            seq(1/looks, 1, length.out = looks)
        }
    } else {
        stopifnot(
            is.numeric(timing),
            length(timing) >= 1,
            all(is.finite(timing)),
            all(timing > 0),
            all(timing <= 1)
        )
        if (any(diff(timing) <= 0)) {
            stop("information fractions in 'timing' must be strictly increasing")
        }
        if (!isTRUE(all.equal(utils::tail(timing, 1), 1,
                              tolerance = sqrt(.Machine$double.eps)))) {
            stop("information fractions in 'timing' must end at 1")
        }
        looks <- length(timing)
    }

    list(type = "timing", looks = looks, timing = timing,
         lookMinN = lookMinN)
}

.bfseq_fixed_schedule_range <- function(n, nrange, lookMinN = 2) {
    stopifnot(
        length(n) == 1,
        is.numeric(n),
        is.finite(n),
        length(nrange) == 2,
        is.numeric(nrange),
        all(is.finite(nrange)),
        nrange[1] >= 2,
        nrange[2] >= nrange[1],
        length(lookMinN) == 1,
        is.numeric(lookMinN),
        is.finite(lookMinN),
        lookMinN >= 2
    )
    n <- as.integer(ceiling(n))
    lower <- max(as.integer(ceiling(nrange[1])),
                 as.integer(ceiling(lookMinN)))
    c(min(lower, n), n)
}

.bfseq_schedule_n <- function(maxN, schedule) {
    maxN <- as.integer(ceiling(maxN))

    if (schedule$type == "increase") {
        if (schedule$minN > maxN) {
            stop("the first-look sample size must be smaller than or equal to the maximum sample size")
        }
        n <- seq(schedule$minN, maxN, by = schedule$by)
        n <- unique(as.integer(c(n[n < maxN], maxN)))
    } else {
        n <- as.integer(ceiling(maxN*schedule$timing - sqrt(.Machine$double.eps)))
        n[length(n)] <- maxN
    }

    .bfseq_validate_schedule(n, minN = schedule$lookMinN)
    n
}

.bfseq_validate_schedule <- function(n, minN = 2) {
    if (any(n < minN)) {
        stop("the generated sample-size schedule contains sample sizes smaller than the minimum required look size")
    }
    if (length(n) > 1 && any(diff(n) <= 0)) {
        stop("the generated sample-size schedule is not strictly increasing")
    }
    invisible(TRUE)
}

.bfseq_minimum_max_n <- function(schedule) {
    if (schedule$type == "increase") {
        return(schedule$minN)
    }

    timing <- schedule$timing
    lower <- max(schedule$lookMinN,
                 floor((schedule$lookMinN - 1)/timing[1]) + 1)
    if (length(timing) > 1) {
        lower <- max(lower, floor(max(1/diff(timing))))
    }
    lower <- as.integer(lower)

    ## The formula above is a conservative lower bound. Verify because ceiling
    ## at exact integer information levels can still create duplicate looks.
    while (TRUE) {
        valid <- try(.bfseq_schedule_n(lower, schedule), silent = TRUE)
        if (!inherits(valid, "try-error")) {
            return(lower)
        }
        lower <- lower + 1L
    }
}

.bfseq_search_bounds <- function(nrange, schedule) {
    stopifnot(
        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1],
        nrange[1] > 0
    )

    lower <- max(as.integer(ceiling(nrange[1])), .bfseq_minimum_max_n(schedule))
    if (!identical(schedule$type, "increase")) {
        lower <- max(lower, .bfseq_default_search_lower(nrange))
    }
    upper <- as.integer(ceiling(nrange[2]))
    if (lower > upper) {
        stop("the lower sample-size search bound exceeds the upper bound after applying the look schedule")
    }

    c(lower, upper)
}

.bfseq_target_probability <- function(design, target) {
    if (target == "H1") {
        return(utils::tail(design$cumpH1, 1))
    }
    utils::tail(design$cumpH0, 1)
}

.bfseq_next_search_candidate <- function(currentN, maximumN) {
    min(maximumN, max(currentN + 1L, as.integer(ceiling(2*currentN))))
}

.bfseq_candidate_invalid <- function(message, reason = "invalid",
                                     terminal = TRUE, power = NA_real_,
                                     result = NULL) {
    stopifnot(
        length(message) == 1,
        is.character(message),
        !is.na(message),
        length(reason) == 1,
        is.character(reason),
        !is.na(reason),
        length(terminal) == 1,
        is.logical(terminal),
        !is.na(terminal)
    )
    condition <- structure(
        list(message = message, reason = reason, terminal = terminal,
             power = power, result = result),
        class = c("bfseq_candidate_invalid", "error", "condition")
    )
    stop(condition)
}

.bfseq_make_candidate <- function(n, criterion, power, error = NULL,
                                  result = NULL, status = "ok",
                                  reason = NULL, terminal = FALSE) {
    list(n = n, criterion = criterion, power = power, error = error,
         result = result, status = status, reason = reason,
         terminal = terminal)
}

.bfseq_validate_evaluation <- function(value) {
    if (!is.list(value)) {
        stop("sequential search evaluator must return a list", call. = FALSE)
    }
    if (!("power" %in% names(value))) {
        stop("sequential search evaluator result is missing 'power'",
             call. = FALSE)
    }
    if (!("result" %in% names(value))) {
        stop("sequential search evaluator result is missing 'result'",
             call. = FALSE)
    }
    if (!is.numeric(value$power) || length(value$power) != 1) {
        stop("sequential search evaluator 'power' must be a scalar numeric value",
             call. = FALSE)
    }
    invisible(TRUE)
}

.bfseq_candidate_is_finite <- function(candidate) {
    is.finite(candidate$criterion)
}

.bfseq_candidate_reached <- function(candidate) {
    .bfseq_candidate_is_finite(candidate) && candidate$criterion >= 0
}

.bfseq_candidate_is_invalid <- function(candidate) {
    !.bfseq_candidate_is_finite(candidate)
}

.bfseq_candidate_is_terminal <- function(candidate) {
    .bfseq_candidate_is_invalid(candidate) && isTRUE(candidate$terminal)
}

.bfseq_adaptive_invalid_candidate <- function(candidate) {
    if (.bfseq_candidate_is_invalid(candidate) && !isTRUE(candidate$terminal)) {
        candidate$error <- paste0(
            candidate$error,
            "; adaptive sample-size bracketing encountered a transient ",
            "invalid candidate and cannot certify the search beyond that ",
            "point. Use search = \"exhaustive\" to scan past transient ",
            "invalid candidates."
        )
    }
    candidate
}

.bfseq_invalid_scan_error <- function(invalids) {
    terminal <- vapply(invalids, function(x) isTRUE(x$terminal), logical(1))
    transientCount <- sum(!terminal)
    terminalCount <- sum(terminal)
    first <- invalids[[1]]
    last <- invalids[[length(invalids)]]

    paste0(
        "sample-size search encountered ", length(invalids),
        " invalid candidate(s) while scanning the candidate range",
        " (transient: ", transientCount, ", terminal: ", terminalCount, ")",
        "; first invalid n = ", first$n,
        if (!is.null(first$reason)) paste0(" [", first$reason, "]") else "",
        ": ", first$error,
        if (length(invalids) > 1) {
            paste0(
                "; last invalid n = ", last$n,
                if (!is.null(last$reason)) paste0(" [", last$reason, "]") else "",
                ": ", last$error
            )
        } else {
            ""
        },
        "; finite candidates did not reach the target"
    )
}

.bfseq_invalid_scan_candidate <- function(limit, invalids) {
    if (is.null(limit) || .bfseq_candidate_is_invalid(limit)) {
        limit <- invalids[[length(invalids)]]
    }
    .bfseq_make_candidate(
        n = limit$n,
        criterion = NA_real_,
        power = limit$power,
        error = .bfseq_invalid_scan_error(invalids),
        result = limit$result,
        status = "invalid",
        reason = "invalid_scan",
        terminal = any(vapply(invalids, function(x) isTRUE(x$terminal),
                              logical(1)))
    )
}

.bfseq_search <- function(power, target, nrange, schedule, evaluate,
                          progress = NULL,
                          search = c("adaptive", "exhaustive")) {
    search <- match.arg(search)
    progress <- .bfseq_validate_progress(progress)
    bounds <- .bfseq_search_bounds(nrange = nrange, schedule = schedule)
    lowerN <- bounds[1]
    upperLimit <- bounds[2]
    increaseCandidates <- if (identical(schedule$type, "increase")) {
        .bfseq_increase_candidates(bounds = bounds, schedule = schedule)
    } else {
        NULL
    }
    maxEvaluations <- if (is.null(increaseCandidates)) {
        upperLimit - lowerN + 1L
    } else {
        length(increaseCandidates)
    }
    cache <- new.env(parent = emptyenv())
    evaluations <- 0L
    phase <- "initial"
    setPhase <- function(value) {
        phase <<- value
    }

    evalN <- function(n) {
        n <- as.integer(ceiling(n))
        key <- as.character(n)
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }

        evaluations <<- evaluations + 1L
        value <- tryCatch(
            evaluate(n),
            bfseq_candidate_invalid = function(e) e
        )
        if (inherits(value, "bfseq_candidate_invalid")) {
            invalidPower <- if (is.null(value$power)) NA_real_ else value$power
            invalidResult <- if (is.null(value$result)) NULL else value$result
            invalidReason <- if (is.null(value$reason)) "invalid" else value$reason
            invalidTerminal <- if (is.null(value$terminal)) TRUE else
                isTRUE(value$terminal)
            out <- .bfseq_make_candidate(
                n = n, criterion = NA_real_, power = invalidPower,
                error = conditionMessage(value), result = invalidResult,
                status = "invalid", reason = invalidReason,
                terminal = invalidTerminal
            )
        } else {
            .bfseq_validate_evaluation(value)
            achieved <- value$power
            if (!is.finite(achieved)) {
                out <- .bfseq_make_candidate(
                    n = n, criterion = NA_real_, power = achieved,
                    error = "non-finite sequential stopping probability",
                    result = value$result, status = "invalid",
                    reason = "nonfinite_power", terminal = TRUE
                )
            } else {
                out <- .bfseq_make_candidate(
                    n = n, criterion = achieved - power, power = achieved,
                    error = NULL, result = value$result
                )
            }
        }
        assign(key, out, envir = cache)
        .bfseq_call_progress(
            progress = progress,
            info = list(
                event = "evaluate",
                phase = phase,
                n = n,
                evaluations = evaluations,
                maxEvaluations = maxEvaluations,
                target = target,
                targetPower = power,
                actualPower = out$power,
                criterion = out$criterion,
                reached = .bfseq_candidate_reached(out),
                nrange = bounds,
                schedule = .bfseq_schedule_summary(schedule),
                search = search,
                error = out$error,
                status = out$status,
                reason = out$reason,
                terminal = out$terminal
            )
        )
        out
    }

    if (identical(search, "exhaustive")) {
        candidates <- if (is.null(increaseCandidates)) {
            lowerN:upperLimit
        } else {
            increaseCandidates
        }
        return(.bfseq_search_full_range(
            power = power, target = target, nrange = bounds,
            schedule = schedule, evalN = evalN, candidates = candidates,
            search = search, getEvaluations = function() evaluations,
            setPhase = setPhase
        ))
    }

    if (!is.null(increaseCandidates)) {
        return(.bfseq_search_full_range(
            power = power, target = target, nrange = bounds,
            schedule = schedule, evalN = evalN,
            candidates = increaseCandidates, search = search,
            getEvaluations = function() evaluations,
            setPhase = setPhase
        ))
    }

    phase <- "lower"
    lower <- evalN(lowerN)
    if (.bfseq_candidate_reached(lower)) {
        return(.bfseq_solver_result(candidate = lower,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations,
                                    reached = TRUE, search = search,
                                    firstCrossingCertified = TRUE))
    }
    if (.bfseq_candidate_is_invalid(lower)) {
        limit <- .bfseq_adaptive_invalid_candidate(lower)
        return(.bfseq_solver_result(candidate = limit,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations, reached = FALSE,
                                    search = search,
                                    firstCrossingCertified = FALSE))
    }

    phase <- "bracket"
    bracket <- .bfseq_find_bracket(evalN = evalN, lower = lower,
                                   lowerN = lowerN, upperLimit = upperLimit)

    if (is.null(bracket$upper)) {
        limit <- if (.bfseq_candidate_is_invalid(bracket$limit)) {
            .bfseq_adaptive_invalid_candidate(bracket$limit)
        } else {
            bracket$limit
        }
        .bfseq_warn_search_limit(limit)
        return(.bfseq_solver_result(candidate = limit,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations, reached = FALSE,
                                    search = search,
                                    firstCrossingCertified = FALSE))
    }

    phase <- "binary"
    firstScan <- .bfseq_first_scan(search, schedule)
    found <- .bfseq_binary_search(evalN = evalN, lowerN = bracket$lowerN,
                                  upperN = bracket$upper$n,
                                  minimumN = lowerN,
                                  firstScan = firstScan, setPhase = setPhase)

    .bfseq_solver_result(candidate = found$candidate,
                         target = target, targetPower = power,
                         nrange = bounds, schedule = schedule,
                         evaluations = evaluations, reached = found$reached,
                         search = search, firstCrossingCertified =
                             found$firstCrossingCertified)
}

.bfseq_find_bracket <- function(evalN, lower, lowerN, upperLimit) {
    currentN <- lowerN
    lastFinite <- if (.bfseq_candidate_is_finite(lower)) lower else NULL

    while (currentN < upperLimit) {
        candidateN <- .bfseq_next_search_candidate(currentN, upperLimit)
        candidate <- evalN(candidateN)

        if (.bfseq_candidate_is_finite(candidate)) {
            lastFinite <- candidate
            if (.bfseq_candidate_reached(candidate)) {
                return(list(lowerN = currentN, upper = candidate,
                            limit = NULL))
            }
            if (candidateN == upperLimit) {
                return(list(lowerN = currentN, upper = NULL,
                            limit = candidate))
            }
            currentN <- candidateN
            next
        }

        if (!.bfseq_candidate_is_terminal(candidate)) {
            return(list(lowerN = currentN, upper = NULL, limit = candidate))
        }

        boundary <- .bfseq_search_before_invalid(evalN = evalN,
                                                 validN = currentN,
                                                 invalidN = candidateN)
        if (!is.null(boundary$upper)) {
            return(list(lowerN = currentN, upper = boundary$upper,
                        limit = NULL))
        }
        if (!is.null(boundary$limit)) {
            return(list(lowerN = currentN, upper = NULL,
                        limit = boundary$limit))
        }
        if (!is.null(lastFinite)) {
            return(list(lowerN = currentN, upper = NULL, limit = lastFinite))
        }
        return(list(lowerN = currentN, upper = NULL, limit = candidate))
    }

    list(lowerN = lowerN, upper = NULL,
         limit = if (!is.null(lastFinite)) lastFinite else lower)
}

.bfseq_search_before_invalid <- function(evalN, validN, invalidN) {
    lastFinite <- evalN(validN)
    lastInvalid <- evalN(invalidN)
    solved <- if (.bfseq_candidate_reached(lastFinite)) lastFinite else NULL

    while ((invalidN - validN) > 1) {
        midpoint <- floor((validN + invalidN)/2)
        current <- evalN(midpoint)

        if (.bfseq_candidate_is_finite(current)) {
            validN <- midpoint
            lastFinite <- current
            if (.bfseq_candidate_reached(current)) {
                solved <- current
            }
        } else {
            invalidN <- midpoint
            lastInvalid <- current
        }
    }

    if (!is.null(solved)) {
        return(list(upper = solved, limit = NULL))
    }
    if (!is.null(lastInvalid) && .bfseq_candidate_is_invalid(lastInvalid)) {
        return(list(upper = NULL, limit = lastInvalid))
    }
    if (.bfseq_candidate_is_finite(lastFinite)) {
        return(list(upper = NULL, limit = lastFinite))
    }
    list(upper = NULL, limit = NULL)
}

.bfseq_increase_candidates <- function(bounds, schedule) {
    stopifnot(identical(schedule$type, "increase"))
    lowerN <- as.integer(bounds[1])
    upperN <- as.integer(bounds[2])
    grid <- seq(schedule$minN, upperN, by = schedule$by)
    grid <- grid[grid >= lowerN & grid <= upperN]
    sort(unique(as.integer(c(lowerN, grid, upperN))))
}

.bfseq_search_full_range <- function(power, target, nrange, schedule, evalN,
                                     candidates, search, getEvaluations,
                                     setPhase) {
    setPhase("scan")
    found <- NULL
    limit <- NULL
    invalids <- list()
    skippedBeforeFound <- FALSE
    for (candidateN in candidates) {
        current <- evalN(candidateN)
        limit <- current
        if (.bfseq_candidate_is_invalid(current)) {
            invalids[[length(invalids) + 1L]] <- current
            if (.bfseq_candidate_is_terminal(current)) {
                break
            }
            next
        }
        if (.bfseq_candidate_reached(current)) {
            found <- current
            skippedBeforeFound <- length(invalids) > 0
            break
        }
    }

    if (is.null(found)) {
        if (length(invalids) > 0) {
            limit <- .bfseq_invalid_scan_candidate(limit = limit,
                                                   invalids = invalids)
        }
        .bfseq_warn_search_limit(limit)
        return(.bfseq_solver_result(candidate = limit,
                                    target = target, targetPower = power,
                                    nrange = nrange, schedule = schedule,
                                    evaluations = getEvaluations(),
                                    reached = FALSE, search = search,
                                    firstCrossingCertified = FALSE))
    }

    .bfseq_solver_result(candidate = found,
                         target = target, targetPower = power,
                         nrange = nrange, schedule = schedule,
                         evaluations = getEvaluations(),
                         reached = TRUE, search = search,
                         firstCrossingCertified = !skippedBeforeFound)
}

.bfseq_warn_search_limit <- function(limit) {
    if (is.null(limit)) {
        return(invisible(NULL))
    }
    if (!is.null(limit$error)) {
        return(invisible(NULL))
    }
    warning("upper bound of sample size search range ('nrange') leads to lower power than specified",
            call. = FALSE)
    invisible(NULL)
}

.bfseq_binary_search <- function(evalN, lowerN, upperN, minimumN,
                                 firstScan = c("none", "adaptive", "exhaustive"),
                                 setPhase = NULL) {
    firstScan <- match.arg(firstScan)
    if (!is.null(setPhase)) {
        setPhase("binary")
    }
    while ((upperN - lowerN) > 1) {
        midpoint <- floor((lowerN + upperN)/2)
        current <- evalN(midpoint)
        if (.bfseq_candidate_is_invalid(current)) {
            return(list(
                candidate = .bfseq_adaptive_invalid_candidate(current),
                reached = FALSE,
                firstCrossingCertified = FALSE
            ))
        }
        if (.bfseq_candidate_reached(current)) {
            upperN <- midpoint
        } else {
            lowerN <- midpoint
        }
    }

    if (!is.null(setPhase)) {
        setPhase("backward")
    }
    foundN <- .bfseq_local_first_success(evalN = evalN, foundN = upperN,
                                         minimumN = minimumN)

    if (firstScan == "adaptive" && foundN > minimumN) {
        if (!is.null(setPhase)) {
            setPhase("probe")
        }
        foundN <- .bfseq_adaptive_first_success(evalN = evalN,
                                                foundN = foundN,
                                                minimumN = minimumN)
    } else if (firstScan == "exhaustive" && foundN > minimumN) {
        if (!is.null(setPhase)) {
            setPhase("scan")
        }
        for (candidateN in minimumN:foundN) {
            current <- evalN(candidateN)
            if (.bfseq_candidate_reached(current)) {
                foundN <- candidateN
                break
            }
        }
    }

    list(candidate = evalN(foundN), reached = TRUE,
         firstCrossingCertified = firstScan != "adaptive")
}

.bfseq_local_first_success <- function(evalN, foundN, minimumN) {
    while (foundN > minimumN) {
        previous <- evalN(foundN - 1L)
        if (!.bfseq_candidate_reached(previous)) {
            break
        }
        foundN <- foundN - 1L
    }
    foundN
}

.bfseq_adaptive_first_success <- function(evalN, foundN, minimumN) {
    repeat {
        earlierN <- .bfseq_probe_earlier_success(evalN = evalN,
                                                 foundN = foundN,
                                                 minimumN = minimumN)
        if (is.null(earlierN) || earlierN >= foundN) {
            break
        }
        foundN <- .bfseq_local_first_success(evalN = evalN,
                                             foundN = earlierN,
                                             minimumN = minimumN)
        if (foundN <= minimumN) {
            break
        }
    }
    foundN
}

.bfseq_probe_earlier_success <- function(evalN, foundN, minimumN) {
    step <- 1L
    repeat {
        candidateN <- max(minimumN, foundN - step)
        current <- evalN(candidateN)
        if (.bfseq_candidate_reached(current)) {
            return(candidateN)
        }
        if (candidateN <= minimumN) {
            return(NULL)
        }
        step <- min(foundN - minimumN, step*2L)
    }
}

.bfseq_solver_result <- function(candidate, target, targetPower, nrange,
                                 schedule, evaluations, reached, search,
                                 firstCrossingCertified) {
    n <- if (isTRUE(reached)) candidate$n else NaN
    status <- if (is.null(candidate$status)) {
        if (.bfseq_candidate_is_invalid(candidate)) "invalid" else "ok"
    } else {
        candidate$status
    }
    reason <- if (is.null(candidate$reason)) NULL else candidate$reason
    terminal <- isTRUE(candidate$terminal)
    result <- candidate$result
    if (!is.null(result)) {
        result$solver <- list(
            n = n,
            maximumN = candidate$n,
            target = target,
            targetPower = targetPower,
            actualPower = candidate$power,
            reached = isTRUE(reached),
            nrange = nrange,
            schedule = .bfseq_schedule_summary(schedule),
            evaluations = evaluations,
            search = search,
            firstCrossingCertified = isTRUE(reached) && isTRUE(firstCrossingCertified),
            status = status,
            reason = reason,
            terminal = terminal,
            error = candidate$error
        )
    }

    list(
        n = n,
        maximumN = candidate$n,
        target = target,
        targetPower = targetPower,
        actualPower = candidate$power,
        reached = isTRUE(reached),
        nrange = nrange,
        schedule = .bfseq_schedule_summary(schedule),
        evaluations = evaluations,
        search = search,
        firstCrossingCertified = isTRUE(reached) && isTRUE(firstCrossingCertified),
        status = status,
        reason = reason,
        terminal = terminal,
        error = candidate$error,
        result = result
    )
}

.bfseq_fixed_solver <- function(n, target, design, schedule) {
    n <- ceiling(n)
    list(
        n = n,
        maximumN = n,
        target = target,
        targetPower = NA_real_,
        actualPower = .bfseq_target_probability(design, target),
        reached = NA,
        nrange = c(n, n),
        schedule = .bfseq_schedule_summary(schedule),
        evaluations = 1L,
        search = NA_character_,
        firstCrossingCertified = NA,
        status = "ok",
        reason = NULL,
        terminal = FALSE,
        error = NULL
    )
}

.bfseq_z_schedule_evaluator <- function(k1, k0, usd, pm, psd, dpm, dpsd,
                                         type, target, schedule, strict,
                                         dots) {
    relpm <- if (type == "moment") NULL else pm
    oneCritical <- (type == "normal" && psd == 0) || type == "directional"
    boundaryCache <- new.env(parent = emptyenv())
    stageCache <- new.env(parent = emptyenv())

    getBoundary <- function(n) {
        key <- as.character(n)
        if (exists(key, envir = boundaryCache, inherits = FALSE)) {
            return(get(key, envir = boundaryCache, inherits = FALSE))
        }
        se <- usd/sqrt(n)
        out <- list(
            n = n,
            se = se,
            zk0 = zcrit(k = k0, se = se, mu = relpm, tau = psd,
                        type = type),
            zk1 = zcrit(k = k1, se = se, mu = relpm, tau = psd,
                        type = type)
        )
        assign(key, out, envir = boundaryCache)
        out
    }

    evalStage <- function(n, bounds = NULL) {
        ## The whole prefix is the cache key because terminal-stage regions
        ## depend on all previous continuation regions.
        key <- paste(n, collapse = "\r")
        if (exists(key, envir = stageCache, inherits = FALSE)) {
            return(get(key, envir = stageCache, inherits = FALSE))
        }

        if (is.null(bounds)) {
            bounds <- lapply(n, getBoundary)
        }
        out <- .bfseq_stage_probabilities_from_bounds(
            bounds = bounds, oneCritical = oneCritical, strict = strict,
            direction = NULL, dpm = dpm, dpsd = dpsd, dots = dots
        )
        assign(key, out, envir = stageCache)
        out
    }

    function(maxN) {
        n <- .bfseq_schedule_n(maxN = maxN, schedule = schedule)
        design <- .bfseq_build_z_design(
            k1 = k1, k0 = k0, se = usd/sqrt(n), n = n, pm = relpm,
            psd = psd, dpm = dpm, dpsd = dpsd, type = type,
            strict = strict, dots = dots,
            getBoundary = function(i) getBoundary(n[[i]]),
            evalStage = function(i, bounds) evalStage(n[seq_len(i)], bounds)
        )

        list(result = design,
             power = .bfseq_target_probability(design = design,
                                                target = target))
    }
}

.bfseq_t_schedule_evaluator <- function(k1, k0, plocation, pscale, pdf,
                                         dpm, dpsd, type, alternative, target,
                                         ratio, schedule, strict, trange,
                                         tail.eps = 1e-3,
                                         tail.nquad = .tbf01_tail_nquad_default,
                                         dots) {
    oneCritical <- alternative != "two.sided"
    regionDirection <- if (alternative == "greater") {
        "positive"
    } else if (alternative == "less") {
        "negative"
    } else {
        NULL
    }
    adaptiveOneSided <- oneCritical && !is.numeric(trange) &&
        trange == "adaptive"
    boundaryCache <- new.env(parent = emptyenv())
    stageCache <- new.env(parent = emptyenv())

    getBoundary <- function(n1, look = NA_integer_) {
        n2 <- if (type == "two.sample") {
            as.integer(ceiling(n1*ratio))
        } else {
            n1
        }
        key <- paste(n1, n2, sep = "\r")
        if (exists(key, envir = boundaryCache, inherits = FALSE)) {
            return(get(key, envir = boundaryCache, inherits = FALSE))
        }

        if (type == "two.sample") {
            neff <- 1/(1/n1 + 1/n2)
        } else {
            neff <- n1
        }
        se <- 1/sqrt(neff)
        searchLimit <- if (adaptiveOneSided) {
            .bfpwr_one_sided_tail_limits(
                origin = 0, step_scale = 1, mean = dpm/se,
                sd = sqrt(1 + (dpsd/se)^2), tail.eps = tail.eps
            )
        } else {
            NULL
        }
        zk0Result <- .bfpwr_tcrit_result(
            k = k0, n1 = n1, n2 = n2, plocation = plocation,
            pscale = pscale, pdf = pdf, alternative = alternative,
            type = type, trange = trange, search_limit = searchLimit,
            tail.nquad = tail.nquad
        )
        zk1Result <- .bfpwr_tcrit_result(
            k = k1, n1 = n1, n2 = n2, plocation = plocation,
            pscale = pscale, pdf = pdf, alternative = alternative,
            type = type, trange = trange, search_limit = searchLimit,
            tail.nquad = tail.nquad
        )
        zk0Message <- .bfseq_t_boundary_status_message(
            list(zk0Result), boundary = "H0", looks = look
        )
        if (!is.null(zk0Message)) {
            .bfseq_candidate_invalid(
                zk0Message, reason = "t_boundary", terminal = TRUE
            )
        }
        zk1Message <- .bfseq_t_boundary_status_message(
            list(zk1Result), boundary = "H1", looks = look
        )
        if (!is.null(zk1Message)) {
            .bfseq_candidate_invalid(
                zk1Message, reason = "t_boundary", terminal = TRUE
            )
        }
        unhandled <- .bfpwr_tcrit_unhandled_warnings(list(zk0Result,
                                                          zk1Result))
        for (msg in unhandled) {
            warning(msg, call. = FALSE)
        }

        out <- list(
            n1 = n1,
            n2 = n2,
            se = se,
            zk0 = zk0Result$value,
            zk1 = zk1Result$value,
            warnings = sum(c(zk0Result$status, zk1Result$status) ==
                               "tail_cutoff"),
            impossible = as.integer(zk0Result$status == "impossible")
        )
        assign(key, out, envir = boundaryCache)
        out
    }

    evalStage <- function(n1, bounds = NULL) {
        ## The whole prefix is the cache key because terminal-stage regions
        ## depend on all previous continuation regions.
        key <- paste(n1, collapse = "\r")
        if (exists(key, envir = stageCache, inherits = FALSE)) {
            return(get(key, envir = stageCache, inherits = FALSE))
        }

        if (is.null(bounds)) {
            bounds <- lapply(seq_along(n1), function(i) getBoundary(n1[[i]],
                                                                     look = i))
        }
        out <- .bfseq_stage_probabilities_from_bounds(
            bounds = bounds, oneCritical = oneCritical, strict = strict,
            direction = regionDirection, dpm = dpm, dpsd = dpsd, dots = dots
        )
        assign(key, out, envir = stageCache)
        out
    }

    function(maxN) {
        n1 <- .bfseq_schedule_n(maxN = maxN, schedule = schedule)
        bounds <- lapply(seq_along(n1), function(i) getBoundary(n1[[i]],
                                                                 look = i))
        n2 <- vapply(bounds, `[[`, numeric(1), "n2")
        searchLimitWarnings <- sum(vapply(bounds, `[[`, integer(1),
                                            "warnings"))
        impossibleWarnings <- sum(vapply(bounds, `[[`, integer(1),
                                          "impossible"))
        if (impossibleWarnings > 0) {
            warning(paste0(
                "No H0 sequential t stopping boundary exists in ",
                impossibleWarnings,
                " boundary search(es); the corresponding H0 stopping ",
                "regions are treated as empty."
            ), call. = FALSE)
        }
        if (searchLimitWarnings > 0) {
            warning(paste0(
                "Adaptive t critical-value search reached the predictive ",
                "tail cutoff in ",
                searchLimitWarnings,
                " sequential boundary search(es); each unresolved boundary ",
                "has marginal tail probability <= ", format(tail.eps),
                ". Pass a wider numeric 'trange' interval to search exact ",
                "bounds."
            ))
        }

        design <- .bfseq_build_t_design(
            k1 = k1, k0 = k0, bounds = bounds, dpm = dpm, dpsd = dpsd,
            plocation = plocation, pscale = pscale, pdf = pdf,
            alternative = alternative, type = type, trange = trange,
            strict = strict, tail.eps = tail.eps, tail.nquad = tail.nquad,
            dots = dots,
            evalStage = function(i, bounds) evalStage(n1[seq_len(i)], bounds)
        )

        list(result = design,
             power = .bfseq_target_probability(design = design,
                                                target = target))
    }
}

.bfseq_schedule_summary <- function(schedule) {
    if (schedule$type == "increase") {
        return(list(type = schedule$type, minN = schedule$minN,
                    by = schedule$by, lookMinN = schedule$lookMinN))
    }
    list(type = schedule$type, looks = schedule$looks, timing = schedule$timing,
         lookMinN = schedule$lookMinN)
}

.bfseq_first_scan <- function(search, schedule) {
    if (identical(schedule$type, "timing") && length(schedule$timing) > 1) {
        return(search)
    }
    "none"
}

.bfseq_ratio_look_min_n <- function(ratio) {
    stopifnot(
        length(ratio) == 1,
        is.numeric(ratio),
        is.finite(ratio),
        ratio > 0
    )
    max(2L, as.integer(floor(1/ratio) + 1L))
}

.bfseq_validate_progress <- function(progress) {
    if (is.null(progress)) {
        return(NULL)
    }
    if (!is.function(progress)) {
        stop("argument 'progress' must be NULL or a function")
    }
    progress
}

.bfseq_extract_progress <- function(dots) {
    dotNames <- names(dots)
    hasProgress <- rep(FALSE, length(dots))
    if (!is.null(dotNames)) {
        hasProgress <- !is.na(dotNames) & dotNames == "progress"
    }
    if (sum(hasProgress) > 1) {
        stop("argument 'progress' matched multiple values", call. = FALSE)
    }
    progress <- if (any(hasProgress)) {
        dots[[which(hasProgress)[1L]]]
    } else {
        NULL
    }
    list(
        progress = .bfseq_validate_progress(progress),
        dots = dots[!hasProgress]
    )
}

.bfseq_call_progress <- function(progress, info) {
    if (is.null(progress)) {
        return(invisible(NULL))
    }

    callbackFormals <- formals(progress)
    if (length(callbackFormals) == 0) {
        progress()
    } else {
        progress(info)
    }
    invisible(NULL)
}

.bfseq_match_vector_arg <- function(arg, choices, name) {
    arg <- as.character(arg)
    if (length(arg) < 1) {
        stop(paste0("argument '", name, "' should be one of ",
                    paste(shQuote(choices), collapse = ", ")))
    }
    bad <- is.na(arg) | !arg %in% choices
    if (any(bad)) {
        stop(paste0("argument '", name, "' should be one of ",
                    paste(shQuote(choices), collapse = ", ")))
    }
    arg
}
