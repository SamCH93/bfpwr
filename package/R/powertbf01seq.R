#' @title Power and Maximum Sample Size Calculations for Sequential t-Test Bayes Factors
#'
#' @description Computes cumulative stopping probabilities for a sequential
#'     t-test Bayes factor design, or determines the maximum sample size needed
#'     to obtain a target stopping probability.
#'
#' @details This function provides a higher-level interface to
#'     \code{\link{ptbf01seq}} and \code{\link{ntbf01seq}}. The analysis and
#'     design prior locations are centered at \code{null} before calling
#'     \code{\link{ptbf01seq}}. If \code{power} is supplied, the returned
#'     design is evaluated at the searched maximum sample size, and
#'     \code{solver$reached} records whether the target was achieved within
#'     \code{nrange}. If \code{n} is supplied, no search is performed and the
#'     \code{solver} element records the achieved stopping probability for the
#'     fixed schedule.
#'
#' @inheritParams ntbf01seq
#' @inheritParams powertbf01
#' @param n Maximum sample size in group 1 for two-sample tests, or maximum
#'     sample size for one-sample and paired tests. Has to be \code{NULL} if
#'     \code{power} is specified. Defaults to \code{NULL}.
#' @param ... Additional arguments passed to \code{\link{ptbf01seq}}.
#'
#' @return An object of class \code{"bfseqdesign"} containing the sequential
#'     design, augmented with a \code{solver} element. In fixed-\code{n} mode,
#'     \code{solver$targetPower} and \code{solver$reached} are \code{NA}.
#'
#' @author František Bartoš
#'
#' @seealso \link{ptbf01seq}, \link{ntbf01seq}, \link{powertbf01}
#'
#' @examples
#' powertbf01seq(n = 20, k1 = 1/2, k0 = 2, dpm = 0.5, dpsd = 0,
#'               alternative = "greater", looks = 2, strict = FALSE)
#' powertbf01seq(power = 0.4, k1 = 1/2, k0 = 2, dpm = 0.5, dpsd = 0,
#'               alternative = "greater", looks = 2, nrange = c(2, 80),
#'               strict = FALSE)
#'
#' @export
powertbf01seq <- function(n = NULL, power = NULL, k1 = 1/10, k0 = 1/k1,
                          null = 0, plocation = 0, pscale = 1/sqrt(2),
                          pdf = 1,
                          type = c("two.sample", "one.sample", "paired"),
                          alternative = c("two.sided", "less", "greater"),
                          dpm = plocation, dpsd = pscale,
                          target = c("h1", "h0"), nrange = c(2, 10^4),
                          looks = 1, timing = NULL, minN = NULL, by = NULL,
                          ratio = 1, strict = TRUE, trange = "adaptive",
                          nextend = 0,
                          search = c("adaptive", "exhaustive"),
                          progress = NULL, ...) {
    if (is.null(n) == is.null(power)) {
        stop("exactly one of 'n' and 'power' must be NULL")
    }
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    target <- match.arg(target)
    search <- match.arg(search)
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

        length(ratio) == 1,
        is.numeric(ratio),
        is.finite(ratio),
        ratio > 0
    )
    nextend <- .bfseq_normalize_nextend(nextend)
    progress <- .bfseq_validate_progress(progress)

    if (is.null(n)) {
        solver <- ntbf01seq.(
            k1 = k1, k0 = k0, power = power, null = null,
            plocation = plocation, pscale = pscale, pdf = pdf,
            dpm = dpm, dpsd = dpsd, type = type,
            alternative = alternative, target = target, nrange = nrange,
            looks = looks, timing = timing, minN = minN, by = by,
            ratio = ratio, strict = strict, trange = trange, integer = TRUE,
            nextend = nextend, search = search, details = TRUE,
            progress = progress, ...
        )
        design <- solver$result
        if (is.null(design)) {
            stop("no valid sequential design could be computed within 'nrange'")
        }
    } else {
        stopifnot(
            length(n) == 1,
            is.numeric(n),
            is.finite(n),
            n >= 2
        )
        lookMinN <- if (type == "two.sample") .bfseq_ratio_look_min_n(ratio) else 2
        schedule <- .bfseq_schedule_spec(looks = looks, timing = timing,
                                         minN = minN, by = by,
                                         nrange = c(2, max(2, n)),
                                         lookMinN = lookMinN)
        n1 <- .bfseq_schedule_n(maxN = n, schedule = schedule)
        n2 <- if (type == "two.sample") {
            as.integer(ceiling(n1*ratio))
        } else {
            n1
        }
        .bfseq_validate_schedule(n2)
        design <- ptbf01seq(
            k1 = k1, k0 = k0, n1 = n1, n2 = n2,
            plocation = plocation - null, pscale = pscale, pdf = pdf,
            dpm = dpm - null, dpsd = dpsd, type = type,
            alternative = alternative, strict = strict, trange = trange, ...
        )
        solver <- .bfseq_fixed_solver(n = n, target = target, design = design,
                                      schedule = schedule, nextend = nextend)
        design$solver <- solver
    }

    design$null <- null
    design$ratio <- ratio
    design
}
