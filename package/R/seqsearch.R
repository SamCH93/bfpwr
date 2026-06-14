## Helper functions for sequential BF sample-size searches
## -----------------------------------------------------------------------------

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
    upper <- as.integer(ceiling(nrange[2]))
    if (lower > upper) {
        stop("the lower sample-size search bound exceeds the upper bound after applying the look schedule")
    }

    c(lower, upper)
}

.bfseq_target_probability <- function(design, target) {
    if (target == "h1") {
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

.bfseq_certification_invalid_candidate <- function(candidate) {
    if (.bfseq_candidate_is_invalid(candidate) && !isTRUE(candidate$terminal)) {
        candidate$error <- paste0(
            candidate$error,
            "; sample-size stability certification encountered a transient ",
            "invalid candidate and cannot certify the returned sample size."
        )
    }
    candidate
}

.bfseq_search <- function(power, target, nrange, schedule, evaluate,
                          nextend = 0, progress = NULL,
                          search = c("adaptive", "exhaustive")) {
    search <- match.arg(search)
    nextend <- .bfseq_normalize_nextend(nextend)
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
            nextend = nextend, search = search,
            getEvaluations = function() evaluations, setPhase = setPhase
        ))
    }

    if (!is.null(increaseCandidates)) {
        return(.bfseq_search_increase(
            power = power, target = target, nrange = bounds,
            schedule = schedule, evalN = evalN,
            candidates = increaseCandidates, nextend = nextend,
            search = search, getEvaluations = function() evaluations,
            setPhase = setPhase
        ))
    }

    phase <- "lower"
    lower <- evalN(lowerN)
    if (.bfseq_candidate_reached(lower)) {
        phase <- "certify"
        certified <- .bfseq_certify_nextend(evalN = evalN, foundN = lower$n,
                                            upperLimit = upperLimit,
                                            nextend = nextend)
        return(.bfseq_solver_result(candidate = certified$candidate,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations,
                                    reached = certified$reached,
                                    nextend = nextend, search = search,
                                    firstCrossingCertified = TRUE))
    }
    if (.bfseq_candidate_is_invalid(lower)) {
        limit <- .bfseq_adaptive_invalid_candidate(lower)
        return(.bfseq_solver_result(candidate = limit,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations,
                                    reached = FALSE, nextend = 0,
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
                                    evaluations = evaluations,
                                    reached = FALSE, nextend = 0,
                                    search = search,
                                    firstCrossingCertified = FALSE))
    }

    phase <- "binary"
    firstScan <- .bfseq_first_scan(search, schedule)
    found <- .bfseq_binary_search(evalN = evalN, lowerN = bracket$lowerN,
                                  upperN = bracket$upper$n,
                                  minimumN = lowerN,
                                  upperLimit = upperLimit,
                                  firstScan = firstScan,
                                  nextend = nextend, setPhase = setPhase)

    .bfseq_solver_result(candidate = found$candidate,
                         target = target, targetPower = power,
                         nrange = bounds, schedule = schedule,
                         evaluations = evaluations, reached = found$reached,
                         nextend = nextend, search = search,
                         firstCrossingCertified = found$firstCrossingCertified)
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

.bfseq_search_increase <- function(power, target, nrange, schedule, evalN,
                                   candidates, nextend, search,
                                   getEvaluations, setPhase) {
    .bfseq_search_full_range(
        power = power, target = target, nrange = nrange, schedule = schedule,
        evalN = evalN, candidates = candidates, nextend = nextend,
        search = search, getEvaluations = getEvaluations,
        setPhase = setPhase
    )
}

.bfseq_search_full_range <- function(power, target, nrange, schedule, evalN,
                                     candidates, nextend, search,
                                     getEvaluations, setPhase) {
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
                                    reached = FALSE, nextend = 0,
                                    search = search,
                                    firstCrossingCertified = FALSE))
    }

    setPhase("certify")
    certified <- .bfseq_certify_increase_candidates(
        evalN = evalN, candidates = candidates, foundN = found$n,
        nextend = nextend
    )
    .bfseq_solver_result(candidate = certified$candidate,
                         target = target, targetPower = power,
                         nrange = nrange, schedule = schedule,
                         evaluations = getEvaluations(),
                         reached = certified$reached, nextend = nextend,
                         search = search,
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

.bfseq_certify_increase_candidates <- function(evalN, candidates, foundN,
                                               nextend) {
    reached <- TRUE
    index <- match(foundN, candidates)
    if (is.na(index)) {
        stop("internal error: increase-search solution is not a candidate")
    }

    if (nextend > 0) {
        repeat {
            end <- index + nextend
            if (end > length(candidates)) {
                warning("Power function may still fall below target power, extend sample size search range")
                reached <- FALSE
                break
            }
            checkIndex <- index:end
            checked <- lapply(candidates[checkIndex], evalN)
            invalid <- vapply(checked, .bfseq_candidate_is_invalid, logical(1))
            if (any(invalid)) {
                return(list(
                    candidate = .bfseq_certification_invalid_candidate(
                        checked[[which(invalid)[1]]]
                    ),
                    reached = FALSE
                ))
            }
            criteria <- vapply(checked, `[[`, numeric(1), "criterion")
            if (all(is.finite(criteria) & criteria >= 0)) {
                break
            }
            below <- checkIndex[!is.finite(criteria) | criteria < 0]
            index <- max(below) + 1L
            if (index > length(candidates)) {
                warning("Power function may still fall below target power, extend sample size search range")
                index <- length(candidates)
                reached <- FALSE
                break
            }
        }
    }

    list(candidate = evalN(candidates[index]), reached = reached)
}

.bfseq_binary_search <- function(evalN, lowerN, upperN, minimumN, upperLimit,
                                 firstScan = c("none", "adaptive", "exhaustive"),
                                 nextend = 0,
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

    if (!is.null(setPhase)) {
        setPhase("certify")
    }
    certified <- .bfseq_certify_nextend(evalN = evalN, foundN = foundN,
                                        upperLimit = upperLimit,
                                        nextend = nextend)
    certified$firstCrossingCertified <- firstScan != "adaptive"
    certified
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

.bfseq_certify_nextend <- function(evalN, foundN, upperLimit, nextend) {
    reached <- TRUE
    if (nextend > 0) {
        repeat {
            if (foundN + nextend > upperLimit) {
                warning("Power function may still fall below target power, extend sample size search range")
                reached <- FALSE
                break
            }
            checkN <- foundN:(foundN + nextend)
            checked <- lapply(checkN, evalN)
            invalid <- vapply(checked, .bfseq_candidate_is_invalid, logical(1))
            if (any(invalid)) {
                return(list(
                    candidate = .bfseq_certification_invalid_candidate(
                        checked[[which(invalid)[1]]]
                    ),
                    reached = FALSE
                ))
            }
            criteria <- vapply(checked, `[[`, numeric(1), "criterion")
            if (all(is.finite(criteria) & criteria >= 0)) {
                break
            }
            below <- checkN[!is.finite(criteria) | criteria < 0]
            foundN <- max(below) + 1L
            if (foundN > upperLimit) {
                warning("Power function may still fall below target power, extend sample size search range")
                foundN <- upperLimit
                reached <- FALSE
                break
            }
        }
    }

    list(candidate = evalN(foundN), reached = reached)
}

.bfseq_solver_result <- function(candidate, target, targetPower, nrange,
                                 schedule, evaluations, reached, nextend,
                                 search, firstCrossingCertified) {
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
            nextend = nextend,
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
        nextend = nextend,
        search = search,
        firstCrossingCertified = isTRUE(reached) && isTRUE(firstCrossingCertified),
        status = status,
        reason = reason,
        terminal = terminal,
        error = candidate$error,
        result = result
    )
}

.bfseq_fixed_solver <- function(n, target, design, schedule, nextend = 0) {
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
        nextend = .bfseq_normalize_nextend(nextend),
        search = NA_character_,
        firstCrossingCertified = NA,
        status = "ok",
        reason = NULL,
        terminal = FALSE,
        error = NULL
    )
}

.bfseq_boundary_data <- function(bounds, oneCritical) {
    se <- vapply(bounds, `[[`, numeric(1), "se")
    if (oneCritical) {
        zk0 <- vapply(bounds, function(x) x$zk0[[1]], numeric(1))
        zk1 <- vapply(bounds, function(x) x$zk1[[1]], numeric(1))
    } else {
        zk0 <- do.call(cbind, lapply(bounds, `[[`, "zk0"))
        zk1 <- do.call(cbind, lapply(bounds, `[[`, "zk1"))
    }
    list(se = se, zk0 = zk0, zk1 = zk1)
}

.bfseq_stage_regions <- function(boundaries, oneCritical, strict,
                                 direction = NULL) {
    if (oneCritical) {
        return(.bfseq_genregions1_stage(zcrit0 = boundaries$zk0,
                                        zcrit1 = boundaries$zk1,
                                        direction = direction))
    }
    .bfseq_genregions2_stage(zcrit0 = boundaries$zk0,
                             zcrit1 = boundaries$zk1,
                             strict = strict)
}

.bfseq_stage_stop_probabilities <- function(regions, se, dpm, dpsd, dots) {
    pars <- predpars(se = se, dpm = dpm, dpsd = dpsd)
    pH1 <- do.call(.bfseq_intstage,
                   c(list(stageregions = regions$H1,
                          mean = pars$mean,
                          sigma = pars$sigma),
                     dots))
    pH0 <- do.call(.bfseq_intstage,
                   c(list(stageregions = regions$H0,
                          mean = pars$mean,
                          sigma = pars$sigma),
                     dots))
    if (!is.numeric(pH1) || length(pH1) != 1 || !is.finite(pH1) ||
        !is.numeric(pH0) || length(pH0) != 1 || !is.finite(pH0)) {
        .bfseq_candidate_invalid(
            "non-finite sequential stage probability",
            reason = "stage_probability",
            terminal = TRUE
        )
    }
    list(pH1 = pH1, pH0 = pH0)
}

.bfseq_sample_size_moments <- function(pH1, pH0, n) {
    stopProb <- pH1 + pH0
    EN <- sum(stopProb*n) + (1 - sum(stopProb))*max(n)
    EN2 <- sum(stopProb*n^2) + (1 - sum(stopProb))*max(n^2)
    list(EN = EN, VarN = EN2 - EN^2)
}

.bfseq_z_schedule_evaluator <- function(k1, k0, usd, null, pm, psd, dpm,
                                         dpsd, type, target, schedule,
                                         strict, dots) {
    relpm <- if (type == "moment") NULL else pm - null
    reldpm <- dpm - null
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

    evalStage <- function(n) {
        ## The whole prefix is the cache key because terminal-stage regions
        ## depend on all previous continuation regions.
        key <- paste(n, collapse = "\r")
        if (exists(key, envir = stageCache, inherits = FALSE)) {
            return(get(key, envir = stageCache, inherits = FALSE))
        }

        bounds <- lapply(n, getBoundary)
        boundaries <- .bfseq_boundary_data(bounds = bounds,
                                           oneCritical = oneCritical)
        regions <- .bfseq_stage_regions(boundaries = boundaries,
                                        oneCritical = oneCritical,
                                        strict = strict)
        out <- .bfseq_stage_stop_probabilities(regions = regions,
                                               se = boundaries$se,
                                               dpm = reldpm,
                                               dpsd = dpsd,
                                               dots = dots)
        assign(key, out, envir = stageCache)
        out
    }

    function(maxN) {
        n <- .bfseq_schedule_n(maxN = maxN, schedule = schedule)
        bounds <- lapply(n, getBoundary)
        boundaries <- .bfseq_boundary_data(bounds = bounds,
                                           oneCritical = oneCritical)
        stages <- lapply(seq_along(n), function(i) evalStage(n[seq_len(i)]))
        pH1 <- vapply(stages, `[[`, numeric(1), "pH1")
        pH0 <- vapply(stages, `[[`, numeric(1), "pH0")
        cumpH1 <- cumsum(pH1)
        cumpH0 <- cumsum(pH0)
        cumpInc <- 1 - cumpH1 - cumpH0
        moments <- .bfseq_sample_size_moments(pH1 = pH1, pH0 = pH0, n = n)

        design <- structure(list(
            k1 = k1, k0 = k0, se = boundaries$se, n = n, pm = relpm,
            psd = psd, dpm = reldpm, dpsd = dpsd, type = type,
            strict = strict, test = "z", zk1 = boundaries$zk1,
            zk0 = boundaries$zk0, EN = moments$EN,
            VarN = moments$VarN, cumpH1 = cumpH1, cumpH0 = cumpH0,
            cumpInc = cumpInc
        ), class = "bfseqdesign")

        list(result = design,
             power = .bfseq_target_probability(design = design,
                                                target = target))
    }
}

.bfseq_t_schedule_evaluator <- function(k1, k0, plocation, pscale, pdf,
                                         dpm, dpsd, type, alternative, target,
                                         ratio, schedule, strict, trange,
                                         tail.eps = 1e-3, dots) {
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
            type = type, trange = trange, search_limit = searchLimit
        )
        zk1Result <- .bfpwr_tcrit_result(
            k = k1, n1 = n1, n2 = n2, plocation = plocation,
            pscale = pscale, pdf = pdf, alternative = alternative,
            type = type, trange = trange, search_limit = searchLimit
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

    evalStage <- function(n1) {
        ## The whole prefix is the cache key because terminal-stage regions
        ## depend on all previous continuation regions.
        key <- paste(n1, collapse = "\r")
        if (exists(key, envir = stageCache, inherits = FALSE)) {
            return(get(key, envir = stageCache, inherits = FALSE))
        }

        bounds <- lapply(seq_along(n1), function(i) getBoundary(n1[[i]],
                                                                 look = i))
        boundaries <- .bfseq_boundary_data(bounds = bounds,
                                           oneCritical = oneCritical)
        regions <- .bfseq_stage_regions(boundaries = boundaries,
                                        oneCritical = oneCritical,
                                        strict = strict,
                                        direction = regionDirection)
        out <- .bfseq_stage_stop_probabilities(regions = regions,
                                               se = boundaries$se,
                                               dpm = dpm,
                                               dpsd = dpsd,
                                               dots = dots)
        assign(key, out, envir = stageCache)
        out
    }

    function(maxN) {
        n1 <- .bfseq_schedule_n(maxN = maxN, schedule = schedule)
        bounds <- lapply(seq_along(n1), function(i) getBoundary(n1[[i]],
                                                                 look = i))
        n2 <- vapply(bounds, `[[`, numeric(1), "n2")
        boundaries <- .bfseq_boundary_data(bounds = bounds,
                                           oneCritical = oneCritical)
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

        stages <- lapply(seq_along(n1), function(i) evalStage(n1[seq_len(i)]))
        pH1 <- vapply(stages, `[[`, numeric(1), "pH1")
        pH0 <- vapply(stages, `[[`, numeric(1), "pH0")
        cumpH1 <- cumsum(pH1)
        cumpH0 <- cumsum(pH0)
        cumpInc <- 1 - cumpH1 - cumpH0

        moments1 <- .bfseq_sample_size_moments(pH1 = pH1, pH0 = pH0,
                                               n = n1)
        moments2 <- .bfseq_sample_size_moments(pH1 = pH1, pH0 = pH0,
                                               n = n2)
        design <- structure(list(
            k1 = k1, k0 = k0, n1 = n1, n2 = n2, dpm = dpm,
            dpsd = dpsd, plocation = plocation, pscale = pscale,
            pdf = pdf, alternative = alternative, type = type,
            trange = trange, strict = strict, test = "t",
            tail.eps = tail.eps,
            zk1 = boundaries$zk1, zk0 = boundaries$zk0,
            EN1 = moments1$EN, EN2 = moments2$EN,
            VarN1 = moments1$VarN, VarN2 = moments2$VarN,
            cumpH1 = cumpH1, cumpH0 = cumpH0, cumpInc = cumpInc
        ), class = "bfseqdesign")

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

.bfseq_normalize_nextend <- function(nextend) {
    stopifnot(
        length(nextend) == 1,
        is.numeric(nextend),
        is.finite(nextend),
        nextend >= 0
    )
    as.integer(ceiling(nextend))
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
