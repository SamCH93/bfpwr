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

.bfseq_search <- function(power, target, nrange, schedule, evaluate,
                          nextend = 0) {
    bounds <- .bfseq_search_bounds(nrange = nrange, schedule = schedule)
    lowerN <- bounds[1]
    upperLimit <- bounds[2]
    cache <- new.env(parent = emptyenv())
    evaluations <- 0L

    evalN <- function(n) {
        n <- as.integer(ceiling(n))
        key <- as.character(n)
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }

        evaluations <<- evaluations + 1L
        value <- try(evaluate(n), silent = TRUE)
        if (inherits(value, "try-error")) {
            out <- list(n = n, criterion = NA_real_, power = NA_real_,
                        error = conditionMessage(attr(value, "condition")),
                        result = NULL)
        } else {
            achieved <- value$power
            out <- list(n = n, criterion = achieved - power, power = achieved,
                        error = NULL, result = value$result,
                        schedule = value$schedule)
            if (!is.finite(out$criterion)) {
                out$error <- "non-finite sequential stopping probability"
            }
        }
        assign(key, out, envir = cache)
        out
    }

    lower <- evalN(lowerN)
    if (is.finite(lower$criterion) && lower$criterion >= 0) {
        return(.bfseq_solver_result(evalN = evalN, candidate = lower,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations,
                                    reached = TRUE, nextend = nextend))
    }

    bracket <- .bfseq_find_bracket(evalN = evalN, lower = lower,
                                   lowerN = lowerN, upperLimit = upperLimit)

    if (is.null(bracket$upper)) {
        limit <- bracket$limit
        if (!is.null(limit$error)) {
            warning("upper bound of sample size search range ('nrange') leads to Power = NaN")
        } else {
            warning("upper bound of sample size search range ('nrange') leads to lower power than specified")
        }
        return(.bfseq_solver_result(evalN = evalN, candidate = limit,
                                    target = target, targetPower = power,
                                    nrange = bounds, schedule = schedule,
                                    evaluations = evaluations,
                                    reached = FALSE, nextend = 0))
    }

    found <- .bfseq_binary_search(evalN = evalN, lowerN = bracket$lowerN,
                                  upperN = bracket$upper$n,
                                  minimumN = lowerN,
                                  upperLimit = upperLimit,
                                  scanFirst = .bfseq_needs_first_scan(schedule),
                                  nextend = nextend)

    .bfseq_solver_result(evalN = evalN, candidate = found$candidate,
                         target = target, targetPower = power,
                         nrange = bounds, schedule = schedule,
                         evaluations = evaluations, reached = found$reached,
                         nextend = nextend)
}

.bfseq_find_bracket <- function(evalN, lower, lowerN, upperLimit) {
    currentN <- lowerN
    lastFinite <- if (is.finite(lower$criterion)) lower else NULL

    while (currentN < upperLimit) {
        candidateN <- .bfseq_next_search_candidate(currentN, upperLimit)
        candidate <- evalN(candidateN)

        if (is.finite(candidate$criterion)) {
            lastFinite <- candidate
            if (candidate$criterion >= 0) {
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
    solved <- if (is.finite(lastFinite$criterion) &&
                  lastFinite$criterion >= 0) lastFinite else NULL

    while ((invalidN - validN) > 1) {
        midpoint <- floor((validN + invalidN)/2)
        current <- evalN(midpoint)

        if (is.finite(current$criterion)) {
            validN <- midpoint
            lastFinite <- current
            if (current$criterion >= 0) {
                solved <- current
            }
        } else {
            invalidN <- midpoint
        }
    }

    if (!is.null(solved)) {
        return(list(upper = solved, limit = NULL))
    }
    if (is.finite(lastFinite$criterion)) {
        return(list(upper = NULL, limit = lastFinite))
    }
    list(upper = NULL, limit = NULL)
}

.bfseq_binary_search <- function(evalN, lowerN, upperN, minimumN, upperLimit,
                                 scanFirst = FALSE, nextend = 0) {
    while ((upperN - lowerN) > 1) {
        midpoint <- floor((lowerN + upperN)/2)
        current <- evalN(midpoint)
        if (is.finite(current$criterion) && current$criterion >= 0) {
            upperN <- midpoint
        } else {
            lowerN <- midpoint
        }
    }

    foundN <- upperN
    while (foundN > minimumN) {
        previous <- evalN(foundN - 1L)
        if (!is.finite(previous$criterion) || previous$criterion < 0) {
            break
        }
        foundN <- foundN - 1L
    }

    if (scanFirst && foundN > minimumN) {
        for (candidateN in minimumN:foundN) {
            current <- evalN(candidateN)
            if (is.finite(current$criterion) && current$criterion >= 0) {
                foundN <- candidateN
                break
            }
        }
    }

    reached <- TRUE
    if (nextend > 0) {
        repeat {
            if (foundN + as.integer(nextend) > upperLimit) {
                warning("Power function may still fall below target power, extend sample size search range")
                reached <- FALSE
                break
            }
            checkN <- foundN:min(upperLimit, foundN + as.integer(nextend))
            checked <- lapply(checkN, evalN)
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

.bfseq_solver_result <- function(evalN, candidate, target, targetPower, nrange,
                                 schedule, evaluations, reached, nextend) {
    n <- if (isTRUE(reached)) candidate$n else NaN
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
        error = candidate$error,
        result = result
    )
}

.bfseq_schedule_summary <- function(schedule) {
    if (schedule$type == "increase") {
        return(list(type = schedule$type, minN = schedule$minN,
                    by = schedule$by, lookMinN = schedule$lookMinN))
    }
    list(type = schedule$type, looks = schedule$looks, timing = schedule$timing,
         lookMinN = schedule$lookMinN)
}

.bfseq_needs_first_scan <- function(schedule) {
    identical(schedule$type, "timing") && length(schedule$timing) > 1
}

.bfseq_ratio_look_min_n <- function(ratio) {
    max(2L, as.integer(floor(1/ratio) + 1L))
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
