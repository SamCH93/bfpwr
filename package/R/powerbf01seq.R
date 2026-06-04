#' @title Power and Maximum Sample Size Calculations for Sequential z-Test Bayes Factors
#'
#' @description Computes cumulative stopping probabilities for a sequential
#'     z-test Bayes factor design, or determines the maximum sample size needed
#'     to obtain a target stopping probability.
#'
#' @details This function provides a higher-level interface to
#'     \code{\link{pbf01seq}} and \code{\link{nbf01seq}} for continuous data.
#'     The analysis and design prior means are centered at \code{null} before
#'     calling \code{\link{pbf01seq}}. If \code{power} is supplied, the
#'     returned design is evaluated at the searched maximum sample size, and
#'     \code{solver$reached} records whether the target was achieved within
#'     \code{nrange}. If \code{n} is supplied, no search is performed and the
#'     \code{solver} element records the achieved stopping probability for the
#'     fixed schedule.
#'
#' @inheritParams nbf01seq
#' @inheritParams powerbf01
#' @param n Maximum sample size (per group for two-sample tests). Has to be
#'     \code{NULL} if \code{power} is specified. Defaults to \code{NULL}.
#' @param type Type of sampling design. One of \code{"two.sample"},
#'     \code{"one.sample"}, or \code{"paired"}. Defaults to
#'     \code{"two.sample"}.
#' @param bftype Type of z-test Bayes factor. One of \code{"normal"},
#'     \code{"directional"}, or \code{"moment"}. Defaults to \code{"normal"}.
#' @param pm Analysis prior mean. Not taken into account for \code{bftype =
#'     "moment"}.
#' @param psd Analysis prior standard deviation (\code{bftype = "normal"} and
#'     \code{bftype = "directional"}) or scale (\code{bftype = "moment"}).
#' @param nrange Maximum sample size search range over which numerical search
#'     is performed. Defaults to \code{c(2, 10^5)}.
#' @param ... Additional arguments passed to \code{\link{pbf01seq}}.
#'
#' @return An object of class \code{"bfseqdesign"} containing the sequential
#'     design, augmented with a \code{solver} element. In fixed-\code{n} mode,
#'     \code{solver$targetPower} and \code{solver$reached} are \code{NA}.
#'
#' @author František Bartoš
#'
#' @seealso \link{pbf01seq}, \link{nbf01seq}, \link{powerbf01}
#'
#' @examples
#' powerbf01seq(n = 60, k1 = 1/5, k0 = 5, pm = 0, psd = 1,
#'              dpm = 0.5, dpsd = 0, looks = 3)
#' powerbf01seq(power = 0.8, k1 = 1/5, k0 = 5, pm = 0, psd = 1,
#'              dpm = 0.5, dpsd = 0, looks = 3, nrange = c(2, 200))
#'
#' @export
powerbf01seq <- function(n = NULL, power = NULL, k1 = 1/10, k0 = 1/k1,
                         sd = 1, null = 0, pm, psd,
                         type = c("two.sample", "one.sample", "paired"),
                         bftype = c("normal", "directional", "moment"),
                         dpm = pm, dpsd = psd, target = c("h1", "h0"),
                         nrange = c(2, 10^5), looks = 1, timing = NULL,
                         minN = NULL, by = NULL, strict = TRUE,
                         nextend = 0,
                         search = c("adaptive", "exhaustive"),
                         progress = NULL, ...) {
    pmMissing <- missing(pm)
    dpmMissing <- missing(dpm)

    if (is.null(n) == is.null(power)) {
        stop("exactly one of 'n' and 'power' must be NULL")
    }
    stopifnot(
        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        0 < sd,

        length(null) == 1,
        is.numeric(null),
        is.finite(null)
    )
    type <- match.arg(type)
    bftype <- match.arg(bftype)
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
        k0 >= 1
    )
    nextend <- .bfseq_normalize_nextend(nextend)
    progress <- .bfseq_validate_progress(progress)
    if (bftype == "moment") {
        if (pmMissing) {
            pm <- NULL
        }
        if (dpmMissing) {
            stop("argument 'dpm' must be specified when bftype = \"moment\"")
        }
    } else if (pmMissing) {
        stop("argument 'pm' is missing, with no default")
    }

    usd <- if (type == "two.sample") sqrt(2)*sd else sd

    if (is.null(n)) {
        solver <- nbf01seq.(
            k1 = k1, k0 = k0, power = power, usd = usd, null = null,
            pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, type = bftype,
            target = target, nrange = nrange, looks = looks, timing = timing,
            minN = minN, by = by, strict = strict, integer = TRUE,
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
        schedule <- .bfseq_schedule_spec(looks = looks, timing = timing,
                                         minN = minN, by = by,
                                         nrange = c(2, max(2, n)))
        nseq <- .bfseq_schedule_n(maxN = n, schedule = schedule)
        relpm <- if (bftype == "moment") NULL else pm - null
        design <- pbf01seq(
            k1 = k1, k0 = k0, se = usd/sqrt(nseq), n = nseq,
            pm = relpm, psd = psd, dpm = dpm - null, dpsd = dpsd,
            type = bftype, strict = strict, ...
        )
        solver <- .bfseq_fixed_solver(n = n, target = target, design = design,
                                      schedule = schedule, nextend = nextend)
        design$solver <- solver
    }

    design$null <- null
    design$sd <- sd
    design$sample.type <- type
    design
}
