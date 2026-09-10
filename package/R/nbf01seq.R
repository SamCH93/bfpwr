nbf01seq. <- function(k1, k0 = 1/k1, power, usd = sqrt(2), null = 0,
                      pm = NULL, psd, dpm = pm, dpsd = psd,
                      type = c("normal", "directional", "moment"),
                      target = c("H1", "H0"), nrange = c(2, 10^5),
                      looks = 1, timing = NULL, minN = NULL, by = NULL,
                      strict = TRUE, integer = TRUE,
                      search = c("adaptive", "exhaustive"),
                      details = FALSE, progress = NULL, ...) {
    type <- match.arg(type)
    target <- match.arg(target)
    search <- match.arg(search)
    if (type == "moment" && is.null(dpm)) {
        stop("argument 'dpm' must be specified when type = \"moment\"")
    }

    ## input checks
    stopifnot(
        length(k1) == 1,
        is.numeric(k1),
        is.finite(k1),
        k1 > 0,
        k1 < 1,

        length(k0) == 1,
        is.numeric(k0),
        is.finite(k0),
        k0 > 1,

        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

        length(usd) == 1,
        is.numeric(usd),
        is.finite(usd),
        0 < usd,

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        psd >= 0,

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        dpsd >= 0,

        length(strict) == 1,
        is.logical(strict),
        !is.na(strict),

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer),

        length(details) == 1,
        is.logical(details),
        !is.na(details)
    )
    if (type != "moment") {
        stopifnot(
            length(pm) == 1,
            is.numeric(pm),
            is.finite(pm)
        )
    }
    if (type != "normal") {
        stopifnot(psd > 0)
    }
    progress <- .bfseq_validate_progress(progress)

    schedule <- .bfseq_schedule_spec(looks = looks, timing = timing,
                                     minN = minN, by = by, nrange = nrange)
    evalDesign <- .bfseq_z_schedule_evaluator(
        k1 = k1, k0 = k0, usd = usd, null = null, pm = pm, psd = psd,
        dpm = dpm, dpsd = dpsd, type = type, target = target,
        schedule = schedule, strict = strict, dots = list(...)
    )

    solver <- .bfseq_search(power = power, target = target, nrange = nrange,
                            schedule = schedule, evaluate = evalDesign,
                            progress = progress, search = search)

    if (details) {
        return(solver)
    }
    if (!is.null(solver$error)) {
        warning(solver$error, call. = FALSE)
    }
    if (integer) {
        return(ceiling(solver$n))
    }
    solver$n
}


#' @title Maximum Sample Size Determination for Sequential z-Test Bayes Factors
#'
#' @description Computes the maximum sample size required for a sequential
#'     z-test Bayes factor design to reach a target probability of stopping for
#'     \eqn{H_1}{H1} or \eqn{H_0}{H0}.
#'
#' @details The function searches over the maximum sample size of the
#'     sequential design. Candidate look schedules are rebuilt for each
#'     maximum sample size according to \code{looks}/\code{timing} or
#'     \code{by}/\code{minN}. Power need not be monotone for every design-prior
#'     and target combination, and rounded interim looks add further
#'     non-monotonicity; \code{search} selects the search rule. For
#'     \code{by}/\code{minN} schedules, scheduled maximum
#'     sample sizes are scanned in increasing order and previous look
#'     calculations are reused; \code{search} is ignored. If the target is not
#'     reached within \code{nrange}, the function returns \code{NaN} and issues
#'     a warning.
#'
#' @inheritParams pbf01seq
#' @param k1 Bayes factor threshold in favor of \eqn{H_1}{H1}. Evidence for
#'     \eqn{H_1}{H1} is obtained when \eqn{\mathrm{BF}_{01} \leq k1}.
#' @param k0 Bayes factor threshold in favor of \eqn{H_0}{H0}. Evidence for
#'     \eqn{H_0}{H0} is obtained when \eqn{\mathrm{BF}_{01} \geq k0}.
#' @param power Target stopping probability.
#' @param usd Unit standard deviation, the standard error of the parameter
#'     estimate at \eqn{n = 1}{n = 1}.
#' @param target Character string. Either \code{"H1"} for the final cumulative
#'     probability of stopping for \eqn{H_1}{H1}, or \code{"H0"} for the final
#'     cumulative probability of stopping for \eqn{H_0}{H0}.
#' @param nrange Maximum sample size search range over which numerical search
#'     is performed. Defaults to \code{c(2, 10^5)}.
#' @param looks Number of sequential looks when \code{timing} is not supplied.
#'     The default \code{looks = 1} gives a fixed design.
#' @param timing Optional cumulative information fractions for the looks. Values
#'     must be positive, strictly increasing, no larger than one, and end at
#'     one. Candidate look sample sizes are
#'     \code{ceiling(maximumN * timing - sqrt(.Machine$double.eps))} with the
#'     final look fixed to \code{maximumN}.
#' @param minN First-look sample size when \code{by} is supplied. Defaults to
#'     the lower search bound.
#' @param by Optional sample-size increase between looks. If supplied, schedules
#'     are generated as \code{seq(minN, maximumN, by = by)} with
#'     \code{maximumN} appended as the final look.
#' @param search Sample-size search rule for schedules generated from
#'     \code{looks}/\code{timing}. \code{"adaptive"} uses bracketing and binary
#'     search with local checks. It is a heuristic that can skip an isolated
#'     crossing on a non-monotone power curve, and stops with a diagnostic if a
#'     transient invalid candidate invalidates the bracket. Use
#'     \code{"exhaustive"} when the first crossing must be guaranteed.
#'     \code{"exhaustive"} scans the
#'     full feasible candidate maximum-sample-size range and returns the first
#'     finite candidate that reaches \code{power}, skipping transient invalid
#'     candidates and stopping at terminal invalid candidates. Ignored when
#'     \code{by} is supplied, where scheduled maximum sample sizes
#'     \code{minN + j * by} are scanned in increasing order. Defaults to
#'     \code{"adaptive"}.
#' @param details Logical indicating whether the full search result should be
#'     returned instead of only the maximum sample size. The detailed result is
#'     scalar; vectorized inputs require \code{details = FALSE}. Defaults to
#'     \code{FALSE}.
#' @param integer Logical indicating whether only integer-valued maximum sample
#'     sizes should be returned. If \code{TRUE}, the required maximum sample
#'     size is rounded to the next larger integer. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{pbf01seq}}.
#'
#' @return The maximum sample size found to achieve the specified power. If
#'     \code{details = TRUE}, returns a list with the sample size, achieved
#'     power, generated design, and search diagnostics.
#'
#' @author František Bartoš
#'
#' @seealso \link{pbf01seq}, \link{powerbf01seq}, \link{nbf01}
#'
#' @examples
#' nbf01seq(k1 = 1/5, k0 = 5, power = 0.8, usd = sqrt(2),
#'          pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
#'          looks = 3, nrange = c(2, 200))
#'
#' @export
nbf01seq <- function(k1, k0 = 1/k1, power, usd = sqrt(2), null = 0,
                     pm = NULL, psd, dpm = pm, dpsd = psd,
                     type = c("normal", "directional", "moment"),
                     target = c("H1", "H0"), nrange = c(2, 10^5),
                     looks = 1, timing = NULL, minN = NULL, by = NULL,
                     strict = TRUE, integer = TRUE,
                     search = c("adaptive", "exhaustive"),
                     details = FALSE, ...) {
    progressInfo <- .bfseq_extract_progress(list(...))
    progress <- progressInfo$progress
    dots <- progressInfo$dots
    type <- if (missing(type)) {
        "normal"
    } else {
        .bfseq_match_vector_arg(type, c("normal", "directional", "moment"),
                                "type")
    }
    target <- if (missing(target)) {
        "H1"
    } else {
        .bfseq_match_vector_arg(target, c("H1", "H0"), "target")
    }

    if (isTRUE(details)) {
        if (length(type) != 1 || length(target) != 1) {
            stop("'details = TRUE' requires scalar 'type' and 'target'")
        }
        return(do.call(nbf01seq., c(list(
            k1 = k1, k0 = k0, power = power, usd = usd, null = null,
            pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, type = type,
            target = target, nrange = nrange, looks = looks,
            timing = timing, minN = minN, by = by, strict = strict,
            integer = integer, search = search, details = TRUE,
            progress = progress
        ), dots)))
    }

    vectorizeArgs <- c("k1", "k0", "power", "usd", "null", "psd", "dpm",
                       "dpsd", "type", "target", "integer")
    if (!is.null(pm)) {
        vectorizeArgs <- c(vectorizeArgs, "pm")
    }
    f <- Vectorize(FUN = nbf01seq., vectorize.args = vectorizeArgs)
    do.call(f, c(list(k1 = k1, k0 = k0, power = power, usd = usd,
      null = null, pm = pm, psd = psd, dpm = dpm, dpsd = dpsd,
      type = type, target = target, nrange = nrange, looks = looks,
      timing = timing, minN = minN, by = by, strict = strict,
      integer = integer, search = search, details = FALSE,
      progress = progress), dots))
}
