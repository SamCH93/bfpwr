nbf01seq. <- function(k1, k0 = 1/k1, power, usd = sqrt(2), null = 0,
                      pm = NULL, psd, dpm = pm, dpsd = psd,
                      type = c("normal", "directional", "moment"),
                      target = c("h1", "h0"), nrange = c(2, 10^5),
                      looks = 1, timing = NULL, minN = NULL, by = NULL,
                      strict = TRUE, integer = TRUE, nextend = 0,
                      details = FALSE, ...) {
    type <- match.arg(type)
    target <- match.arg(target)
    if (type == "moment" && is.null(dpm)) {
        stop("argument 'dpm' must be specified when type = \"moment\"")
    }

    ## input checks
    stopifnot(
        length(k1) == 1,
        is.numeric(k1),
        is.finite(k1),
        k1 > 0,
        k1 <= 1,

        length(k0) == 1,
        is.numeric(k0),
        is.finite(k0),
        k0 >= 1,

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
    nextend <- .bfseq_normalize_nextend(nextend)

    schedule <- .bfseq_schedule_spec(looks = looks, timing = timing,
                                     minN = minN, by = by, nrange = nrange)
    evalDesign <- function(maxN) {
        n <- .bfseq_schedule_n(maxN = maxN, schedule = schedule)
        relpm <- if (type == "moment") NULL else pm - null
        design <- pbf01seq(
            k1 = k1, k0 = k0, se = usd/sqrt(n), n = n,
            pm = relpm, psd = psd, dpm = dpm - null, dpsd = dpsd,
            type = type, strict = strict, ...
        )
        list(result = design,
             power = .bfseq_target_probability(design = design,
                                                target = target))
    }

    solver <- .bfseq_search(power = power, target = target, nrange = nrange,
                            schedule = schedule, evaluate = evalDesign,
                            nextend = nextend)

    if (details) {
        return(solver)
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
#'     \code{by}/\code{minN}. For multi-look timing schedules, the search
#'     verifies the first maximum sample size in \code{nrange} that reaches the
#'     requested stopping probability, because the rounded interim looks can
#'     make the power curve non-monotone. If the target is not reached within
#'     \code{nrange}, the function returns \code{NaN} and issues a warning.
#'
#' @inheritParams pbf01seq
#' @inheritParams nbf01
#' @param k1 Bayes factor threshold in favor of \eqn{H_1}{H1}. Evidence for
#'     \eqn{H_1}{H1} is obtained when \eqn{\mathrm{BF}_{01} \leq k1}.
#' @param k0 Bayes factor threshold in favor of \eqn{H_0}{H0}. Evidence for
#'     \eqn{H_0}{H0} is obtained when \eqn{\mathrm{BF}_{01} \geq k0}.
#' @param null Point null value. The sequential z-test calculation is performed
#'     after centering the analysis and design prior means at \code{null}.
#' @param target Character string. Either \code{"h1"} for the final cumulative
#'     probability of stopping for \eqn{H_1}{H1}, or \code{"h0"} for the final
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
#' @param nextend Number of sample sizes beyond the solution used to check that
#'     the target probability does not drop below \code{power}. Non-integer
#'     values are rounded up. Defaults to \code{0}.
#' @param details Logical indicating whether the full search result should be
#'     returned instead of only the maximum sample size. The detailed result is
#'     scalar; vectorized inputs require \code{details = FALSE}. Defaults to
#'     \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{pbf01seq}}.
#'
#' @return The required maximum sample size to achieve the specified power. If
#'     \code{details = TRUE}, returns a list containing the sample size, achieved
#'     power, generated design, and search diagnostics.
#'
#' @author Samuel Pawel
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
                     target = c("h1", "h0"), nrange = c(2, 10^5),
                     looks = 1, timing = NULL, minN = NULL, by = NULL,
                     strict = TRUE, integer = TRUE, nextend = 0,
                     details = FALSE, ...) {
    type <- if (missing(type)) {
        "normal"
    } else {
        .bfseq_match_vector_arg(type, c("normal", "directional", "moment"),
                                "type")
    }
    target <- if (missing(target)) {
        "h1"
    } else {
        .bfseq_match_vector_arg(target, c("h1", "h0"), "target")
    }

    if (isTRUE(details)) {
        if (length(type) != 1 || length(target) != 1) {
            stop("'details = TRUE' requires scalar 'type' and 'target'")
        }
        return(nbf01seq.(
            k1 = k1, k0 = k0, power = power, usd = usd, null = null,
            pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, type = type,
            target = target, nrange = nrange, looks = looks, timing = timing,
            minN = minN, by = by, strict = strict, integer = integer,
            nextend = nextend, details = TRUE, ...
        ))
    }

    vectorizeArgs <- c("k1", "k0", "power", "usd", "null", "psd", "dpm",
                       "dpsd", "type", "target", "integer")
    if (!is.null(pm)) {
        vectorizeArgs <- c(vectorizeArgs, "pm")
    }
    f <- Vectorize(FUN = nbf01seq., vectorize.args = vectorizeArgs)
    f(k1 = k1, k0 = k0, power = power, usd = usd, null = null,
      pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, type = type,
      target = target, nrange = nrange, looks = looks, timing = timing,
      minN = minN, by = by, strict = strict, integer = integer,
      nextend = nextend, details = FALSE, ...)
}
