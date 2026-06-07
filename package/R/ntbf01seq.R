ntbf01seq. <- function(k1, k0 = 1/k1, power, null = 0,
                       plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                       dpm = plocation, dpsd = pscale,
                       type = c("two.sample", "one.sample", "paired"),
                       alternative = c("two.sided", "less", "greater"),
                       target = c("h1", "h0"), nrange = c(2, 10^4),
                       looks = 1, timing = NULL, minN = NULL, by = NULL,
                       ratio = 1, strict = TRUE, trange = "adaptive",
                       integer = TRUE, nextend = 0,
                       search = c("adaptive", "exhaustive"),
                       details = FALSE,
                       progress = NULL, ...) {
    search <- match.arg(search)
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

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

        length(plocation) == 1,
        is.numeric(plocation),
        is.finite(plocation),

        length(pscale) == 1,
        is.numeric(pscale),
        is.finite(pscale),
        0 < pscale,

        length(pdf) == 1,
        is.numeric(pdf),
        is.finite(pdf),
        0 < pdf,

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        0 <= dpsd,

        length(ratio) == 1,
        is.numeric(ratio),
        is.finite(ratio),
        ratio > 0,

        length(strict) == 1,
        is.logical(strict),
        !is.na(strict),

        (is.numeric(trange) && length(trange) == 2 && all(is.finite(trange)) &&
         trange[2] > trange[1]) || (is.character(trange) && length(trange) == 1 &&
                                    !is.na(trange) && trange == "adaptive"),

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer),

        length(details) == 1,
        is.logical(details),
        !is.na(details)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    target <- match.arg(target)
    nextend <- .bfseq_normalize_nextend(nextend)
    progress <- .bfseq_validate_progress(progress)

    lookMinN <- if (type == "two.sample") .bfseq_ratio_look_min_n(ratio) else 2
    schedule <- .bfseq_schedule_spec(looks = looks, timing = timing,
                                     minN = minN, by = by, nrange = nrange,
                                     lookMinN = lookMinN)
    evalDesign <- .bfseq_t_schedule_evaluator(
        k1 = k1, k0 = k0, plocation = plocation - null,
        pscale = pscale, pdf = pdf, dpm = dpm - null,
        dpsd = dpsd, type = type, alternative = alternative,
        target = target, ratio = ratio, schedule = schedule,
        strict = strict, trange = trange, dots = list(...)
    )

    solver <- .bfseq_search(power = power, target = target, nrange = nrange,
                            schedule = schedule, evaluate = evalDesign,
                            nextend = nextend, progress = progress,
                            search = search)

    if (details) {
        return(solver)
    }
    if (integer) {
        return(ceiling(solver$n))
    }
    solver$n
}


#' @title Maximum Sample Size Determination for Sequential t-Test Bayes Factors
#'
#' @description Computes the maximum sample size required for a sequential
#'     t-test Bayes factor design to reach a target probability of stopping for
#'     \eqn{H_1}{H1} or \eqn{H_0}{H0}.
#'
#' @details The function searches over the maximum group-1 sample size for
#'     two-sample designs, or the maximum sample size for one-sample and paired
#'     designs. Candidate look schedules are rebuilt for each maximum sample
#'     size according to \code{looks}/\code{timing} or \code{by}/\code{minN}.
#'     For multi-look timing schedules, rounded interim looks can make the
#'     power curve non-monotone; \code{search} selects the search rule. For
#'     \code{by}/\code{minN} schedules, scheduled maximum sample sizes are
#'     scanned in increasing order and previous look calculations are reused;
#'     \code{search} is ignored. If the target is not reached within
#'     \code{nrange}, the function returns \code{NaN} and issues a warning.
#'
#' @inheritParams ptbf01seq
#' @inheritParams ntbf01
#' @inheritParams nbf01seq
#' @param k1 Bayes factor threshold in favor of \eqn{H_1}{H1}. Evidence for
#'     \eqn{H_1}{H1} is obtained when \eqn{\mathrm{BF}_{01} \leq k1}.
#' @param k0 Bayes factor threshold in favor of \eqn{H_0}{H0}. Evidence for
#'     \eqn{H_0}{H0} is obtained when \eqn{\mathrm{BF}_{01} \geq k0}.
#' @param ratio Allocation ratio \code{n2 / n1} for two-sample designs.
#'     Candidate group-2 sample sizes are \code{ceiling(n1 * ratio)}. Ignored
#'     for one-sample and paired designs.
#' @param details Logical indicating whether the full search result should be
#'     returned instead of only the maximum sample size. The detailed result is
#'     scalar; vectorized inputs require \code{details = FALSE}. Defaults to
#'     \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{ptbf01seq}}.
#'
#' @return The maximum sample size in group 1 found to achieve the specified
#'     power. If \code{details = TRUE}, returns a list with the sample size,
#'     achieved power, generated design, and search diagnostics.
#'
#' @author František Bartoš
#'
#' @seealso \link{ptbf01seq}, \link{powertbf01seq}, \link{ntbf01}
#'
#' @examples
#' ntbf01seq(k1 = 1/2, k0 = 2, power = 0.4, dpm = 0.5, dpsd = 0,
#'           alternative = "greater", looks = 2, nrange = c(2, 80),
#'           strict = FALSE)
#'
#' @export
ntbf01seq <- function(k1, k0 = 1/k1, power, null = 0,
                      plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                      dpm = plocation, dpsd = pscale,
                      type = c("two.sample", "one.sample", "paired"),
                      alternative = c("two.sided", "less", "greater"),
                      target = c("h1", "h0"), nrange = c(2, 10^4),
                      looks = 1, timing = NULL, minN = NULL, by = NULL,
                      ratio = 1, strict = TRUE, trange = "adaptive",
                      integer = TRUE, nextend = 0,
                      search = c("adaptive", "exhaustive"),
                      details = FALSE,
                      progress = NULL, ...) {
    type <- if (missing(type)) {
        "two.sample"
    } else {
        .bfseq_match_vector_arg(type, c("two.sample", "one.sample", "paired"),
                                "type")
    }
    alternative <- if (missing(alternative)) {
        "two.sided"
    } else {
        .bfseq_match_vector_arg(alternative, c("two.sided", "less",
                                               "greater"), "alternative")
    }
    target <- if (missing(target)) {
        "h1"
    } else {
        .bfseq_match_vector_arg(target, c("h1", "h0"), "target")
    }

    if (isTRUE(details)) {
        if (length(type) != 1 || length(alternative) != 1 ||
            length(target) != 1) {
            stop("'details = TRUE' requires scalar 'type', 'alternative', and 'target'")
        }
        return(ntbf01seq.(
            k1 = k1, k0 = k0, power = power, null = null,
            plocation = plocation, pscale = pscale, pdf = pdf,
            dpm = dpm, dpsd = dpsd, type = type,
            alternative = alternative, target = target, nrange = nrange,
            looks = looks, timing = timing, minN = minN, by = by,
            ratio = ratio, strict = strict, trange = trange,
            integer = integer, nextend = nextend, search = search,
            details = TRUE, progress = progress, ...
        ))
    }

    f <- Vectorize(FUN = ntbf01seq.,
                   vectorize.args = c("k1", "k0", "power", "null",
                                      "plocation", "pscale", "pdf",
                                      "dpm", "dpsd", "type",
                                      "alternative", "target", "ratio",
                                      "integer"))
    f(k1 = k1, k0 = k0, power = power, null = null,
      plocation = plocation, pscale = pscale, pdf = pdf,
      dpm = dpm, dpsd = dpsd, type = type, alternative = alternative,
      target = target, nrange = nrange, looks = looks, timing = timing,
      minN = minN, by = by, ratio = ratio, strict = strict,
      trange = trange, integer = integer, nextend = nextend,
      search = search, details = FALSE, progress = progress, ...)
}
