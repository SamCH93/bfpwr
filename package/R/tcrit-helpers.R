## Evaluate a root function defensively. Failed, non-scalar, or non-finite
## evaluations are treated as NaN so search code can keep classifying failures.
.bfpwr_root_value <- function(f, x) {
    ans <- try(suppressWarnings(f(x)), silent = TRUE)
    if (inherits(ans, "try-error") || length(ans) != 1 ||
        !is.numeric(ans) || !is.finite(ans)) {
        return(NaN)
    }
    ans
}

## Find a bracketed root and, when a faster scout function was used, certify
## the result against the final function before accepting it.
.bfpwr_certified_root <- function(f, x0, x1, f0 = NaN, f1 = NaN,
                                  final_fun = f, tolerance = 1e-5, ...) {
    if (x0 == x1) {
        return(structure("degenerate root interval", class = "try-error"))
    }

    if (!is.finite(f0)) {
        f0 <- .bfpwr_root_value(f = f, x = x0)
    }
    if (!is.finite(f1)) {
        f1 <- .bfpwr_root_value(f = f, x = x1)
    }
    if (!is.finite(f0) || !is.finite(f1) || f0*f1 > 0) {
        return(structure("root not bracketed", class = "try-error"))
    }
    if (f0 == 0) {
        if (identical(final_fun, f)) {
            return(x0)
        }
        final0 <- .bfpwr_root_value(f = final_fun, x = x0)
        if (is.finite(final0) && abs(final0) <= tolerance) {
            return(x0)
        }
    }
    if (f1 == 0) {
        if (identical(final_fun, f)) {
            return(x1)
        }
        final1 <- .bfpwr_root_value(f = final_fun, x = x1)
        if (is.finite(final1) && abs(final1) <= tolerance) {
            return(x1)
        }
    }

    root <- try(stats::uniroot(f = f, interval = sort(c(x0, x1)),
                               extendInt = "no", ...)$root,
                silent = TRUE)
    if (inherits(root, "try-error")) {
        return(root)
    }
    if (identical(final_fun, f)) {
        return(root)
    }

    residual <- .bfpwr_root_value(f = final_fun, x = root)
    if (is.finite(residual) && abs(residual) <= tolerance) {
        return(root)
    }

    final_f0 <- .bfpwr_root_value(f = final_fun, x = x0)
    final_f1 <- .bfpwr_root_value(f = final_fun, x = x1)
    if (!is.finite(final_f0) || !is.finite(final_f1) ||
        final_f0*final_f1 > 0) {
        return(structure("root not certified by final function",
                         class = "try-error"))
    }
    if (final_f0 == 0) return(x0)
    if (final_f1 == 0) return(x1)

    try(stats::uniroot(f = final_fun, interval = sort(c(x0, x1)),
                       extendInt = "no", ...)$root,
        silent = TRUE)
}

## Keep only integrate() controls from dots so root-search controls are not
## accidentally forwarded to numerical integration.
.bfpwr_integrate_dots <- function(dots, rel.tol.default = NULL) {
    if (length(dots) == 0) {
        out <- list()
    } else {
        dot_names <- names(dots)
        if (is.null(dot_names)) {
            out <- list()
        } else {
            integrate_names <- c("subdivisions", "rel.tol", "abs.tol",
                                 "stop.on.error", "keep.xy")
            keep <- nzchar(dot_names) & dot_names %in% integrate_names
            out <- dots[keep]
        }
    }
    if (!is.null(rel.tol.default) && !("rel.tol" %in% names(out))) {
        out$rel.tol <- rel.tol.default
    }
    out
}

## Keep only uniroot() controls from dots for critical-value root searches.
.bfpwr_uniroot_dots <- function(dots) {
    if (length(dots) == 0) {
        return(list())
    }
    dot_names <- names(dots)
    if (is.null(dot_names)) {
        return(list())
    }
    keep_names <- c("tol", "maxiter", "trace", "check.conv")
    keep <- nzchar(dot_names) & dot_names %in% keep_names
    dots[keep]
}

## Structured tcrit issues preserve whether a warning is handled internally or
## should be surfaced to callers.
.bfpwr_tcrit_issue <- function(code, status, message, handled = TRUE) {
    list(
        code = code,
        status = status,
        message = message,
        handled = handled
    )
}

## Use a custom warning class so tcrit_result() can capture numerical status
## without parsing warning text.
.bfpwr_tcrit_condition <- function(message, code, status,
                                   handled = TRUE) {
    structure(
        list(
            message = message,
            call = NULL,
            code = code,
            status = status,
            handled = handled
        ),
        class = c("bfpwr_tcrit_warning", "warning", "condition")
    )
}

## Emit a structured tcrit warning.
.bfpwr_tcrit_warning <- function(message, code, status,
                                 handled = TRUE) {
    warning(.bfpwr_tcrit_condition(message = message, code = code,
                                   status = status, handled = handled))
}

## Convert any warning from tcrit() into a normalized issue record.
.bfpwr_tcrit_condition_issue <- function(warning) {
    if (inherits(warning, "bfpwr_tcrit_warning")) {
        return(.bfpwr_tcrit_issue(
            code = warning$code,
            status = warning$status,
            message = conditionMessage(warning),
            handled = isTRUE(warning$handled)
        ))
    }

    .bfpwr_tcrit_issue(
        code = "warning",
        status = "search_failed",
        message = conditionMessage(warning),
        handled = FALSE
    )
}

## Collapse possibly multiple issues into one boundary-search status, ordered by
## severity.
.bfpwr_tcrit_status_from_result <- function(value, issues) {
    if (length(issues) > 0) {
        statuses <- vapply(issues, `[[`, character(1), "status")
        for (status in c("error", "search_failed", "impossible",
                         "tail_cutoff")) {
            if (any(statuses == status)) {
                return(status)
            }
        }
    }
    if (length(value) > 0 && all(is.finite(value))) {
        return("ok")
    }
    "search_failed"
}

## Pick the issue code that best explains the final status.
.bfpwr_tcrit_reason <- function(status, issues) {
    if (length(issues) > 0) {
        statuses <- vapply(issues, `[[`, character(1), "status")
        match <- which(statuses == status)
        if (length(match) > 0) {
            return(issues[[match[[1]]]]$code)
        }
        return(issues[[1]]$code)
    }
    if (status == "ok") {
        return(NA_character_)
    }
    "nonfinite_value"
}

## Evaluate tcrit() while collecting structured warnings instead of emitting
## them immediately.
.bfpwr_tcrit_eval <- function(args) {
    tcritIssues <- list()
    value <- withCallingHandlers(
        do.call(tcrit, args),
        warning = function(w) {
            tcritIssues[[length(tcritIssues) + 1L]] <<-
                .bfpwr_tcrit_condition_issue(w)
            invokeRestart("muffleWarning")
        }
    )
    status <- .bfpwr_tcrit_status_from_result(value = value,
                                              issues = tcritIssues)
    warnings <- vapply(tcritIssues, `[[`, character(1), "message")
    list(
        value = value,
        status = status,
        reason = .bfpwr_tcrit_reason(status = status, issues = tcritIssues),
        warnings = warnings,
        issues = tcritIssues
    )
}

## Inspect and update collected tcrit issues.
.bfpwr_tcrit_has_issue <- function(result, code) {
    any(vapply(result$issues, function(issue) identical(issue$code, code),
               logical(1)))
}

## Recompute derived status fields after appending an issue.
.bfpwr_tcrit_append_issue <- function(result, issue) {
    result$issues[[length(result$issues) + 1L]] <- issue
    result$warnings <- vapply(result$issues, `[[`, character(1), "message")
    result$status <- .bfpwr_tcrit_status_from_result(
        value = result$value,
        issues = result$issues
    )
    result$reason <- .bfpwr_tcrit_reason(status = result$status,
                                         issues = result$issues)
    result
}

## Run tcrit() and add an extra diagnostic when a numeric two-sided trange
## misses a boundary that the adaptive range can still find.
.bfpwr_tcrit_result <- function(...) {
    args <- list(...)
    trange <- args$trange
    if (is.null(trange)) {
        trange <- "adaptive"
    }

    result <- .bfpwr_tcrit_eval(args = args)
    if (identical(result$status, "impossible") &&
        identical(args$alternative, "two.sided") &&
        is.numeric(trange) &&
        .bfpwr_tcrit_has_issue(result = result,
                               code = "maximum_bf_below_k")) {
        adaptiveArgs <- args
        adaptiveArgs$trange <- "adaptive"
        adaptiveResult <- .bfpwr_tcrit_eval(args = adaptiveArgs)
        if (identical(adaptiveResult$status, "impossible")) {
            result$status <- "impossible"
            result$reason <- "maximum_bf_below_k"
        } else {
            result <- .bfpwr_tcrit_append_issue(
                result = result,
                issue = .bfpwr_tcrit_issue(
                    code = "numeric_trange_miss",
                    status = "search_failed",
                    message = paste0(
                        "numeric 'trange' does not contain the two-sided ",
                        "BF01 = k boundary"
                    )
                )
            )
        }
    }
    result
}

## Return only unhandled tcrit warnings so sequential callers do not repeat
## warnings they already translated into boundary statuses.
.bfpwr_tcrit_unhandled_warnings <- function(results) {
    issues <- unlist(
        lapply(results, function(result) {
            if (is.null(result$issues)) list() else result$issues
        }),
        recursive = FALSE,
        use.names = FALSE
    )
    if (length(issues) == 0) {
        return(character())
    }
    unhandled <- vapply(issues, function(issue) isFALSE(issue$handled),
                        logical(1))
    unique(vapply(issues[unhandled], `[[`, character(1), "message"))
}

## Extract statuses from a list of t boundary searches.
.bfseq_t_boundary_statuses <- function(results) {
    vapply(results, `[[`, character(1), "status")
}

## Convert t-boundary statuses into a sequential-search error message. Missing
## H0 boundaries can be valid empty stopping regions; missing H1 boundaries are
## terminal failures.
.bfseq_t_boundary_status_message <- function(results, boundary,
                                             looks = seq_along(results)) {
    stopifnot(boundary %in% c("H0", "H1"))
    if (length(looks) != length(results)) {
        stop("internal error: boundary status look labels do not match results",
             call. = FALSE)
    }

    statuses <- .bfseq_t_boundary_statuses(results)
    allowed <- if (boundary == "H0") {
        c("ok", "impossible", "tail_cutoff")
    } else {
        c("ok", "tail_cutoff")
    }
    invalid <- !(statuses %in% allowed)
    if (!any(invalid)) {
        return(NULL)
    }

    look <- which(invalid)[1]
    warningText <- results[[look]]$warnings
    detail <- if (length(warningText) > 0) {
        paste(unique(warningText), collapse = "; ")
    } else {
        paste0("status: ", statuses[[look]])
    }
    paste0(
        "Failed to compute ", boundary,
        " sequential t stopping boundary at look ", looks[[look]], ": ",
        detail,
        ". Widen numeric 'trange' or use adaptive 'trange' with a smaller ",
        "'tail.eps'."
    )
}

## Public-facing fixed-design sequential calls validate t boundaries eagerly.
.bfseq_validate_t_boundary_statuses <- function(results, boundary,
                                               looks = seq_along(results)) {
    msg <- .bfseq_t_boundary_status_message(results = results,
                                            boundary = boundary,
                                            looks = looks)
    if (!is.null(msg)) {
        stop(msg, call. = FALSE)
    }
    invisible(.bfseq_t_boundary_statuses(results))
}

## Search-mode sequential calls warn for tolerated H0 empty regions and tail
## cutoffs, but still propagate unhandled numerical warnings.
.bfseq_warn_t_boundary_statuses <- function(results0, results1, tail.eps) {
    statuses0 <- .bfseq_t_boundary_statuses(results0)
    statuses1 <- .bfseq_t_boundary_statuses(results1)

    unhandled <- .bfpwr_tcrit_unhandled_warnings(c(results0, results1))
    for (msg in unhandled) {
        warning(msg, call. = FALSE)
    }

    impossible <- sum(statuses0 == "impossible")
    if (impossible > 0) {
        warning(paste0(
            "No H0 sequential t stopping boundary exists in ",
            impossible,
            " boundary search(es); the corresponding H0 stopping regions ",
            "are treated as empty."
        ), call. = FALSE)
    }

    tailCutoff <- sum(c(statuses0, statuses1) == "tail_cutoff")
    if (tailCutoff > 0) {
        warning(paste0(
            "Adaptive t critical-value search reached the predictive tail ",
            "cutoff in ",
            tailCutoff,
            " sequential boundary search(es); each unresolved boundary has ",
            "marginal tail probability <= ", format(tail.eps),
            ". Pass a wider numeric 'trange' interval to search exact bounds."
        ), call. = FALSE)
    }
}

## Determine the expected one-sided search direction from the BF value at the
## origin and the requested alternative.
.bfpwr_one_sided_direction <- function(alternative, f_origin) {
    stopifnot(
        alternative %in% c("greater", "less"),
        length(f_origin) == 1,
        is.numeric(f_origin),
        is.finite(f_origin),
        f_origin != 0
    )

    if (alternative == "greater") {
        if (f_origin > 0) 1 else -1
    } else {
        if (f_origin > 0) -1 else 1
    }
}

## Convert a predictive tail-mass cutoff into a finite one-sided t search limit.
.bfpwr_one_sided_tail_limit <- function(direction, origin, step_scale, mean,
                                        sd, tail.eps) {
    stopifnot(
        length(direction) == 1,
        direction %in% c(-1, 1),
        length(origin) == 1,
        is.numeric(origin),
        is.finite(origin),
        length(step_scale) == 1,
        is.numeric(step_scale),
        is.finite(step_scale),
        step_scale > 0,
        length(mean) == 1,
        is.numeric(mean),
        is.finite(mean),
        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        sd > 0,
        length(tail.eps) == 1,
        is.numeric(tail.eps),
        is.finite(tail.eps),
        tail.eps > 0,
        tail.eps < 0.5
    )

    if (direction > 0) {
        limit <- stats::qnorm(p = tail.eps, mean = mean, sd = sd,
                              lower.tail = FALSE)
        if (!is.finite(limit) || limit <= origin) {
            limit <- origin
        }
        tail_probability <- stats::pnorm(q = limit, mean = mean, sd = sd,
                                         lower.tail = FALSE)
    } else {
        limit <- stats::qnorm(p = tail.eps, mean = mean, sd = sd,
                              lower.tail = TRUE)
        if (!is.finite(limit) || limit >= origin) {
            limit <- origin
        }
        tail_probability <- stats::pnorm(q = limit, mean = mean, sd = sd,
                                         lower.tail = TRUE)
    }

    list(
        direction = direction,
        limit = limit,
        search_limit = abs(limit - origin)/step_scale,
        tail_probability = tail_probability,
        tail.eps = tail.eps
    )
}

## Precompute both directional tail limits; the root direction is selected after
## the BF value at the origin is known.
.bfpwr_one_sided_tail_limits <- function(origin, step_scale, mean, sd,
                                         tail.eps) {
    list(
        positive = .bfpwr_one_sided_tail_limit(
            direction = 1, origin = origin, step_scale = step_scale,
            mean = mean, sd = sd, tail.eps = tail.eps
        ),
        negative = .bfpwr_one_sided_tail_limit(
            direction = -1, origin = origin, step_scale = step_scale,
            mean = mean, sd = sd, tail.eps = tail.eps
        ),
        tail.eps = tail.eps
    )
}

## Accept either a scalar search distance or the directional object produced by
## .bfpwr_one_sided_tail_limits().
.bfpwr_select_one_sided_search_limit <- function(search_limit, direction,
                                                 origin, step_scale) {
    stopifnot(
        length(direction) == 1,
        direction %in% c(-1, 1),
        length(origin) == 1,
        is.numeric(origin),
        is.finite(origin),
        length(step_scale) == 1,
        is.numeric(step_scale),
        is.finite(step_scale),
        step_scale > 0
    )

    if (is.null(search_limit)) {
        return(NULL)
    }

    if (is.numeric(search_limit) && length(search_limit) == 1) {
        stopifnot(is.finite(search_limit), search_limit >= 0)
        return(list(
            direction = direction,
            limit = origin + direction*step_scale*search_limit,
            search_limit = search_limit,
            tail_probability = NA_real_,
            tail.eps = NA_real_
        ))
    }

    stopifnot(is.list(search_limit))
    selected <- if (direction > 0) search_limit$positive else search_limit$negative
    stopifnot(
        is.list(selected),
        length(selected$search_limit) == 1,
        is.numeric(selected$search_limit),
        is.finite(selected$search_limit),
        selected$search_limit >= 0
    )
    selected
}

## Standard return object for one-sided adaptive searches.
.bfpwr_one_sided_adaptive_result <- function(root, search_limit_reached,
                                             selected_limit = NULL,
                                             status = NULL) {
    if (is.null(selected_limit)) {
        selected_limit <- list(
            direction = NA_real_,
            limit = NA_real_,
            search_limit = NA_real_,
            tail_probability = NA_real_,
            tail.eps = NA_real_
        )
    }
    if (is.null(status)) {
        status <- if (!inherits(root, "try-error")) {
            "ok"
        } else if (search_limit_reached) {
            "tail_cutoff"
        } else {
            "search_failed"
        }
    }
    list(
        root = root,
        search_limit_reached = search_limit_reached,
        status = status,
        direction = selected_limit$direction,
        limit = selected_limit$limit,
        search_limit = selected_limit$search_limit,
        tail_probability = selected_limit$tail_probability,
        tail.eps = selected_limit$tail.eps
    )
}

## Classify the finite tail cutoff when no sign change was found before it.
.bfpwr_one_sided_limit_status <- function(f_origin, f_limit, certify_fun,
                                          direction, origin, step_scale,
                                          search_limit,
                                          impossible_margin = 0.05,
                                          flat_tolerance = 0.01,
                                          flat_fraction = 0.05) {
    if (!is.finite(f_origin) || !is.finite(f_limit)) {
        return("search_failed")
    }
    if (f_origin*f_limit <= 0) {
        return("ok")
    }

    ## A same-signed finite predictive cutoff usually means only that the root
    ## was not found before the user-visible tail-mass bound. Classify a missing
    ## one-sided H0 root as unattainable only when a farther certified probe is
    ## flat and still materially away from BF01 = k. This avoids labeling roots
    ## just beyond the finite cutoff as mathematically impossible.
    status <- "tail_cutoff"
    if (!is.finite(search_limit) || search_limit < 16 ||
        abs(f_limit) < impossible_margin) {
        return(status)
    }

    probe_step <- search_limit/2
    if (!is.finite(probe_step) || probe_step <= 0 ||
        probe_step >= search_limit) {
        return(status)
    }
    probe_x <- origin + direction*step_scale*probe_step
    f_probe <- .bfpwr_root_value(f = certify_fun, x = probe_x)
    if (!is.finite(f_probe) || f_origin*f_probe <= 0) {
        return(status)
    }

    flat_bound <- max(flat_tolerance, flat_fraction*abs(f_limit))
    if (abs(f_limit - f_probe) <= flat_bound) {
        return("impossible")
    }
    status
}

## Use a fast scout root only if its residual is certified by the stable BF
## evaluation.
.bfpwr_residual_certified_root <- function(scout_fun, certify_fun, x0, x1,
                                           final_fun = certify_fun,
                                           tolerance = 1e-5, ...) {
    root <- try(stats::uniroot(f = scout_fun, interval = sort(c(x0, x1)),
                               extendInt = "no", ...)$root,
                silent = TRUE)
    if (inherits(root, "try-error") || !is.numeric(root) ||
        length(root) != 1 || !is.finite(root)) {
        return(structure("scout root failed", class = "try-error"))
    }

    residual <- .bfpwr_root_value(f = certify_fun, x = root)
    if (!is.finite(residual) || abs(residual) > tolerance) {
        return(structure("scout root not certified", class = "try-error"))
    }

    final_residual <- .bfpwr_root_value(f = final_fun, x = root)
    if (is.finite(final_residual) && abs(final_residual) <= tolerance) {
        return(root)
    }

    structure("scout root not certified by final function",
              class = "try-error")
}

## One-sided adaptive boundary search. The fast scout function is used only to
## locate candidate brackets. search_fun performs stable bracket checks, while
## returned roots and tail cutoffs must be certified by certify_fun.
## The return value is a list with root and search_limit_reached.
.bfpwr_one_sided_adaptive_root <- function(certify_fun, scout_fun, alternative,
                                           search_fun = certify_fun,
                                           origin = 0, step_scale = 1,
                                           try_opposite = FALSE,
                                           search_limit = NULL,
                                           steps = c(0.1, 0.25, 0.5, 1, 1.5,
                                                     2, 2.5, 3, 3.5),
                                           tail_steps = c(4, 8, 16, 32, 64,
                                                          128, 256),
                                           scout_tail_steps = c(4, 5, 6, 7, 8,
                                                                16, 32, 64,
                                                                128, 256),
                                           scout_tolerance = 1e-5,
                                           ...) {
    f_origin <- .bfpwr_root_value(f = certify_fun, x = origin)
    if (!is.finite(f_origin)) {
        return(.bfpwr_one_sided_adaptive_result(
            root = structure("non-finite root start", class = "try-error"),
            search_limit_reached = FALSE
        ))
    }
    if (f_origin == 0) {
        return(.bfpwr_one_sided_adaptive_result(
            root = origin, search_limit_reached = FALSE
        ))
    }
    f_origin_search <- .bfpwr_root_value(f = search_fun, x = origin)
    if (!is.finite(f_origin_search)) {
        return(.bfpwr_one_sided_adaptive_result(
            root = structure("non-finite search start", class = "try-error"),
            search_limit_reached = FALSE
        ))
    }

    expected_direction <- .bfpwr_one_sided_direction(alternative = alternative,
                                                     f_origin = f_origin)
    if (is.null(search_limit)) {
        return(.bfpwr_one_sided_adaptive_result(
            root = structure("missing finite search limit", class = "try-error"),
            search_limit_reached = FALSE
        ))
    }

    directions <- if (try_opposite) {
        c(expected_direction, -expected_direction)
    } else {
        expected_direction
    }
    search_limit_reached <- FALSE
    search_limit_status <- "search_failed"
    selected_limit <- NULL
    for (direction in directions) {
        selected_limit <- .bfpwr_select_one_sided_search_limit(
            search_limit = search_limit, direction = direction,
            origin = origin, step_scale = step_scale
        )
        if (is.null(selected_limit)) {
            return(.bfpwr_one_sided_adaptive_result(
                root = structure("missing finite search limit",
                                 class = "try-error"),
                search_limit_reached = FALSE
            ))
        }
        current_search_limit <- selected_limit$search_limit
        if (current_search_limit <= 0) {
            search_limit_reached <- TRUE
            search_limit_status <- "tail_cutoff"
            next
        }

        scan_steps <- sort(unique(c(steps, scout_tail_steps)))
        scan_steps <- scan_steps[is.finite(scan_steps) & scan_steps > 0 &
                                 scan_steps <= current_search_limit]
        if (!current_search_limit %in% scan_steps) {
            scan_steps <- sort(c(scan_steps, current_search_limit))
        }
        exact_tail_steps <- sort(unique(c(tail_steps, current_search_limit)))
        exact_tail_steps <- exact_tail_steps[
            is.finite(exact_tail_steps) & exact_tail_steps > 0 &
                exact_tail_steps <= current_search_limit
        ]

        last_finite_step <- 0
        last_finite_x <- origin
        finite_x <- numeric(0)
        scout_prev_x <- origin
        scout_prev_f <- f_origin_search

        ## First use the fast direct integral only as a scout. A finite scout
        ## sign change is never returned unless the stable BF path certifies it.
        for (step in scan_steps) {
            x1 <- origin + direction*step_scale*step
            f1 <- .bfpwr_root_value(f = scout_fun, x = x1)
            if (!is.finite(f1)) {
                break
            }
            if (is.finite(scout_prev_f) && scout_prev_f*f1 <= 0) {
                root <- .bfpwr_residual_certified_root(
                    scout_fun = scout_fun, certify_fun = search_fun,
                    final_fun = certify_fun,
                    x0 = scout_prev_x, x1 = x1,
                    tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
                root <- .bfpwr_certified_root(
                    f = search_fun, x0 = scout_prev_x, x1 = x1,
                    final_fun = certify_fun, tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
            }
            scout_prev_x <- x1
            scout_prev_f <- f1
            last_finite_step <- step
            last_finite_x <- x1
            finite_x <- c(finite_x, x1)
        }

        fprev <- f_origin_search
        xprev <- origin
        certified_limit_finite <- FALSE
        if (last_finite_step > 0) {
            f_last <- .bfpwr_root_value(f = search_fun, x = last_finite_x)
            if (is.finite(f_last) && f_origin_search*f_last <= 0) {
                x_bracket0 <- origin
                f_bracket0 <- f_origin_search
                x_bracket1 <- last_finite_x
                f_bracket1 <- f_last
                if (length(finite_x) > 1) {
                    for (x_candidate in rev(finite_x[-length(finite_x)])) {
                        f_candidate <- .bfpwr_root_value(f = search_fun,
                                                         x = x_candidate)
                        if (!is.finite(f_candidate)) {
                            next
                        }
                        if (f_candidate*f_bracket1 <= 0) {
                            x_bracket0 <- x_candidate
                            f_bracket0 <- f_candidate
                            break
                        }
                        x_bracket1 <- x_candidate
                        f_bracket1 <- f_candidate
                    }
                }
                root <- .bfpwr_certified_root(
                    f = search_fun, x0 = x_bracket0, x1 = x_bracket1,
                    f0 = f_bracket0, f1 = f_bracket1,
                    final_fun = certify_fun, tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
            }
            if (is.finite(f_last)) {
                fprev <- f_last
                xprev <- last_finite_x
                if (last_finite_step >= current_search_limit) {
                    f_limit_final <- .bfpwr_root_value(f = certify_fun,
                                                       x = last_finite_x)
                    if (is.finite(f_limit_final) &&
                        f_origin*f_limit_final <= 0) {
                        root <- .bfpwr_certified_root(
                            f = certify_fun, x0 = origin, x1 = last_finite_x,
                            f0 = f_origin, f1 = f_limit_final, ...
                        )
                        if (!inherits(root, "try-error")) {
                            return(.bfpwr_one_sided_adaptive_result(
                                root = root, search_limit_reached = FALSE,
                                selected_limit = selected_limit
                            ))
                        }
                    }
                    certified_limit_finite <- is.finite(f_limit_final) &&
                        f_origin*f_limit_final > 0
                    if (certified_limit_finite) {
                        search_limit_status <- .bfpwr_one_sided_limit_status(
                            f_origin = f_origin,
                            f_limit = f_limit_final,
                            certify_fun = certify_fun,
                            direction = direction,
                            origin = origin,
                            step_scale = step_scale,
                            search_limit = current_search_limit
                        )
                    }
                }
            }
        }
        if (last_finite_step >= current_search_limit) {
            search_limit_reached <- search_limit_reached ||
                certified_limit_finite
            next
        }

        direction_limit_reached <- FALSE
        exact_steps <- exact_tail_steps[exact_tail_steps > last_finite_step]
        f_limit <- NaN
        limit_checked <- FALSE
        limit_ruled_out <- FALSE
        for (step in exact_steps) {
            x1 <- origin + direction*step_scale*step
            f1 <- if (step >= current_search_limit && is.finite(f_limit)) {
                f_limit
            } else {
                .bfpwr_root_value(f = search_fun, x = x1)
            }
            if (is.finite(f1) && fprev*f1 <= 0) {
                root <- .bfpwr_certified_root(
                    f = search_fun, x0 = xprev, x1 = x1, f0 = fprev,
                    f1 = f1, final_fun = certify_fun,
                    tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
            }
            if (is.finite(f1)) {
                xprev <- x1
                fprev <- f1
                if (step >= current_search_limit) {
                    f_limit_final <- .bfpwr_root_value(f = certify_fun,
                                                       x = x1)
                    if (is.finite(f_limit_final) &&
                        f_origin*f_limit_final <= 0) {
                        root <- .bfpwr_certified_root(
                            f = certify_fun, x0 = origin, x1 = x1,
                            f0 = f_origin, f1 = f_limit_final, ...
                        )
                        if (!inherits(root, "try-error")) {
                            return(.bfpwr_one_sided_adaptive_result(
                                root = root, search_limit_reached = FALSE,
                                selected_limit = selected_limit
                            ))
                        }
                    }
                    direction_limit_reached <- is.finite(f_limit_final) &&
                        f_origin*f_limit_final > 0
                    if (direction_limit_reached) {
                        search_limit_status <- .bfpwr_one_sided_limit_status(
                            f_origin = f_origin,
                            f_limit = f_limit_final,
                            certify_fun = certify_fun,
                            direction = direction,
                            origin = origin,
                            step_scale = step_scale,
                            search_limit = current_search_limit
                        )
                    }
                }
            }
            ## If the exact path is still far from the threshold, check the
            ## finite search limit before walking every tail step. Near-root
            ## cases continue locally instead of jumping to the search limit.
            if (!limit_checked && step < current_search_limit &&
                is.finite(fprev) && abs(fprev) > 0.1) {
                x_limit <- origin + direction*step_scale*current_search_limit
                f_limit <- .bfpwr_root_value(f = search_fun, x = x_limit)
                limit_checked <- TRUE
                if (is.finite(f_limit)) {
                    f_limit_final <- .bfpwr_root_value(f = certify_fun,
                                                       x = x_limit)
                    if (is.finite(f_limit_final) &&
                        f_origin*f_limit_final <= 0) {
                        root <- .bfpwr_certified_root(
                            f = certify_fun, x0 = origin, x1 = x_limit,
                            f0 = f_origin, f1 = f_limit_final, ...
                        )
                        if (!inherits(root, "try-error")) {
                            return(.bfpwr_one_sided_adaptive_result(
                                root = root, search_limit_reached = FALSE,
                                selected_limit = selected_limit
                            ))
                        }
                    }
                    direction_limit_reached <- is.finite(f_limit_final) &&
                        f_origin*f_limit_final > 0
                    if (direction_limit_reached) {
                        search_limit_status <- .bfpwr_one_sided_limit_status(
                            f_origin = f_origin,
                            f_limit = f_limit_final,
                            certify_fun = certify_fun,
                            direction = direction,
                            origin = origin,
                            step_scale = step_scale,
                            search_limit = current_search_limit
                        )
                    }
                    if (direction_limit_reached && fprev*f_limit > 0) {
                        limit_ruled_out <- TRUE
                        break
                    }
                }
            }
        }
        if (limit_ruled_out) {
            search_limit_reached <- TRUE
            next
        }

        search_limit_reached <- search_limit_reached || direction_limit_reached
    }

    .bfpwr_one_sided_adaptive_result(
        root = structure("root not bracketed", class = "try-error"),
        search_limit_reached = search_limit_reached,
        selected_limit = selected_limit,
        status = if (search_limit_reached) search_limit_status else
            "search_failed"
    )
}


#' @title Compute Critical T-Values for T-Test Bayes Factors
#'
#' @description Computes critical t-values for T-Test Bayes factors
#'
#' @param k Positive numeric. Bayes factor threshold (BF01 oriented in favor of
#'     H0)
#' @param n1 Sample size in group 1
#' @param n2 Sample size in group 2 (is ignored for one-sample \eqn{t}-tests)
#' @param plocation \eqn{t} prior location
#' @param pscale \eqn{t} prior scale
#' @param pdf \eqn{t} prior degrees of freedom
#' @param type Type of \eqn{t}-test. Can be \code{"two.sample"},
#'     \code{"one.sample"}, or \code{"paired"}
#' @param alternative Direction of the test. Can be either \code{"two.sided"},
#'     \code{"less"}, or \code{"greater"}. The latter two truncate the analysis
#'     prior to negative and positive effects, respectively
#' @param trange Numerical search strategy. Can be either \code{"adaptive"}
#'     (default) or an interval. One-sided adaptive searches require a finite
#'     \code{search_limit} supplied by the caller.
#' @param search_limit Finite one-sided adaptive search limit, either as a
#'     scalar in \eqn{t}-statistic units or as the directional object returned
#'     by \code{.bfpwr_one_sided_tail_limits()}.
#' @param ... Optional numerical controls. For numeric ranges and two-sided
#'     adaptive searches, arguments are passed to \code{stats::uniroot}. In
#'     adaptive one-sided searches, \code{subdivisions}, \code{rel.tol},
#'     \code{abs.tol}, \code{stop.on.error}, and \code{keep.xy} are used for BF
#'     integration, while \code{tol}, \code{maxiter}, \code{trace}, and
#'     \code{check.conv} are passed to \code{stats::uniroot}.
#'
#' @return Numeric vector of critical t-value(s)
#'
#' @examples
#' tseq <- seq(-10, 10, length.out = 100)
#' n1 <- 50
#' n2 <- 60
#' type <- "two.sample"
#' alternative <- "two.sided"
#' plocation <- 0
#' pscale <- 1/sqrt(2)
#' pdf <- 1
#' k <- 3
#' tcrit1 <- tcrit(k = k, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
#'                 pdf = pdf, alternative = alternative, type = type)
#' plot(tseq, tbf01(t = tseq, n1 = n1, n2 = n2, plocation = plocation,
#'                  pscale = pscale, pdf = pdf, alternative = alternative,
#'                  type = type),
#'      type = "l", xlab = "t-statistic", ylab = bquote("BF"["01"]), log = "y")
#' abline(h = k, lty = 2)
#' abline(v = tcrit1, lty = 2)
#'
#' n1 <- n2 <- 100
#' alternative <- "greater"
#' k <- 6
#' search_limit <- .bfpwr_one_sided_tail_limits(origin = 0, step_scale = 1,
#'                                              mean = 0, sd = 2,
#'                                              tail.eps = 1e-3)
#' tcrit2 <- tcrit(k = k, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
#'                 pdf = pdf, alternative = alternative, type = type,
#'                 search_limit = search_limit)
#' plot(tseq, tbf01(t = tseq, n1 = n1, n2 = n2, plocation = plocation,
#'                  pscale = pscale, pdf = pdf, alternative = alternative,
#'                  type = type),
#'      type = "l", xlab = "t-statistic", ylab = bquote("BF"["01"]), log = "y")
#' abline(h = k, lty = 2)
#' abline(v = tcrit2, lty = 2)
#'
#' @noRd
#'
#' @keywords internal
tcrit <- function(k, n1, n2, plocation, pscale, pdf, type, alternative,
                  trange = "adaptive", search_limit = NULL,
                  tail.nquad = .tbf01_tail_nquad_default, ...) {

    ## determine t-statistic for which BF = k
    dots <- list(...)
    searchDots <- .bfpwr_integrate_dots(dots = dots,
                                        rel.tol.default = 1e-2)
    rootDots <- .bfpwr_uniroot_dots(dots = dots)
    ## Final BF evaluation used to accept returned roots.
    rootFun <- function(t) {
        tbf01(t = t, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
              pdf = pdf, type = type, alternative = alternative,
              log = TRUE, tail.nquad = tail.nquad) - log(k)
    }
    ## Search evaluation with integration controls from dots.
    rootFunSearch <- function(t) {
        do.call(tbf01, c(list(
            t = t, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
            pdf = pdf, type = type, alternative = alternative, log = TRUE,
            tail.nquad = tail.nquad
        ), searchDots)) - log(k)
    }
    if (type == "two.sample") {
        pars <- .tbf01_pars(n1 = n1, n2 = n2, type = type)
    } else {
        pars <- .tbf01_pars(n1 = n1, n2 = n1, type = type)
    }
    region <- .tbf01_prior_region(plocation = plocation, pscale = pscale,
                                  pdf = pdf, alternative = alternative)
    ## Fast direct-integral scout used only to locate candidate brackets.
    rootFunFast <- function(t) {
        do.call(.tbf01_log_fast, c(list(
            t = t, df = pars$df, neff = pars$neff,
            plocation = plocation, pscale = pscale, pdf = pdf,
            region = region
        ), searchDots)) - log(k)
    }

    if (alternative == "two.sided") {
        ## guess search range based on search range from z-test BF
        if (!is.numeric(trange) && trange == "adaptive") {
            if (type == "two.sample") {
                neff <- 1/(1/n1 + 1/n2)
            } else {
                neff <- n1
            }
            se <- 1/sqrt(neff)
            X <- (plocation^2/pscale^2 + log(1 + pscale^2/se^2) -
                  2*log(k))*(1 + se^2/pscale^2)
            if (X <= 0) {
                X <- 5/k
            }
            zcrit <- -plocation*se/pscale^2 + c(-1, 1)*sqrt(X)
            meant <- mean(zcrit)
            if (zcrit[1] < zcrit[2]) {
                searchIntLow <- c(zcrit[1] - 2, meant)
                searchIntUp <- c(meant, zcrit[2] + 2)
            } else {
                searchIntLow <- c(zcrit[2] - 2, meant)
                searchIntUp <- c(meant, zcrit[1] + 2)
            }
        } else {
            meant <- mean(trange)
            searchIntLow <- c(trange[1], meant)
            searchIntUp <- c(meant, trange[2])
        }
        if (k > 1) {
            ## Check impossible H0 boundaries only in the interval searched below.
            ## This avoids unconstrained wrong-tail evaluations in tbf01().
            maxInt <- c(searchIntLow[1], searchIntUp[2])
            opt <- try(stats::optimize(f = function(t) {
                                           ans <- suppressWarnings(rootFun(t))
                                           if (is.finite(ans)) ans else -Inf
                                       },
                                       interval = maxInt,
                                       maximum = TRUE),
                       silent = TRUE)
            if (!inherits(opt, "try-error") &&
                is.finite(opt$objective) && opt$objective < 0) {
                .bfpwr_tcrit_warning(
                    "maximum BF is less than k; BF01 = k impossible",
                    code = "maximum_bf_below_k",
                    status = "impossible"
                )
                return(c(NaN, NaN))
            }
        }
        ## search for critical values
        tcrit <- c(NaN, NaN)
        lower <- try(stats::uniroot(f = rootFun, interval = searchIntLow,
                                    extendInt = "upX", ...)$root,
                     silent = TRUE)
        upper <- try(stats::uniroot(f = rootFun, interval = searchIntUp,
                                    extendInt = "downX", ...)$root,
                     silent = TRUE)
        if (inherits(lower, "try-error") || inherits(upper, "try-error")) {
            .bfpwr_tcrit_warning(
                "Numerical problems: Could not find 2 t-roots",
                code = "two_roots_failed",
                status = "search_failed"
            )
        } else {
            tcrit <- c(lower, upper)
        }
    } else { # one-sided cases
        search_limit_reached <- FALSE
        search_limit_missing <- FALSE
        search_status <- "search_failed"
        if (!is.numeric(trange) && trange == "adaptive") {
            if (is.null(search_limit)) {
                search_limit_missing <- TRUE
                .bfpwr_tcrit_warning(
                    paste0(
                        "Adaptive one-sided t critical-value search requires ",
                        "'search_limit'"
                    ),
                    code = "missing_search_limit",
                    status = "error"
                )
                res <- structure("missing finite search limit",
                                 class = "try-error")
            } else {
                search <- do.call(.bfpwr_one_sided_adaptive_root, c(list(
                    certify_fun = rootFun, search_fun = rootFunSearch,
                    scout_fun = rootFunFast, alternative = alternative,
                    origin = 0, step_scale = 1, try_opposite = FALSE,
                    search_limit = search_limit
                ), rootDots))
                res <- search$root
                search_limit_reached <- search$search_limit_reached
                search_status <- search$status
            }
        } else {
            searchint <- trange
            extend <- "no"

            suppressWarnings({
                res <- try(stats::uniroot(f = rootFun, interval = searchint,
                                          extendInt = extend, ...)$root,
                           silent = TRUE)
            })
        }
        if (inherits(res, "try-error")) {
            if (search_limit_missing) {
                tcrit <- NaN
            } else if (identical(search_status, "impossible")) {
                .bfpwr_tcrit_warning(
                    paste0(
                        "BF01 = k appears unattainable for this one-sided ",
                        "t test; no critical value exists in the searched ",
                        "tail."
                    ),
                    code = "one_sided_unattainable",
                    status = "impossible"
                )
            } else if (search_limit_reached) {
                .bfpwr_tcrit_warning(
                    paste0(
                        "Adaptive t critical-value search reached the ",
                        "predictive tail cutoff without bracketing BF01 = k; ",
                        "pass a wider numeric 'trange' interval to search ",
                        "exact bounds."
                    ),
                    code = "tail_cutoff",
                    status = "tail_cutoff"
                )
            } else {
                .bfpwr_tcrit_warning(
                    "Numerical problems finding critical value",
                    code = "critical_value_failed",
                    status = "search_failed"
                )
            }
            tcrit <- NaN
        } else {
            tcrit <- res
        }
    }
    return(tcrit)
}
