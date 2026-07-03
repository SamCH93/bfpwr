#' @title Power and Maximum Sample Size Calculations for Sequential z-Test Bayes Factors
#'
#' @description Computes cumulative stopping probabilities for a sequential
#'     z-test Bayes factor design, or determines the maximum sample size needed
#'     to obtain a target stopping probability.
#'
#' @details This function provides a higher-level interface to
#'     \code{\link{pbf01seq}} and \code{\link{nbf01seq}} for continuous data.
#'     If \code{power} is supplied, the returned design is evaluated at the
#'     searched maximum sample size, and \code{solver$reached} records whether
#'     the target was achieved within \code{nrange}. If \code{n} is supplied,
#'     no search is performed and the \code{solver} element records the
#'     achieved stopping probability for the fixed schedule. For fixed-\code{n}
#'     increment schedules with missing \code{minN}, \code{nrange[1]} is used
#'     as the default first look, clamped to \code{n} when necessary.
#'
#' @inheritParams nbf01seq
#' @param n Maximum sample size (per group for two-sample tests). Has to be
#'     \code{NULL} if \code{power} is specified. Defaults to \code{NULL}.
#' @param sd Standard deviation of one observation (for \code{type =
#'     "two.sample"} or \code{type = "one.sample"}) or of one difference within
#'     a pair of observations (\code{type = "paired"}). Is assumed to be known.
#'     Defaults to \code{1}.
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
                         dpm = pm, dpsd = psd, target = c("H1", "H0"),
                         nrange = c(2, 10^5), looks = 1, timing = NULL,
                         minN = NULL, by = NULL, strict = TRUE,
                         search = c("adaptive", "exhaustive"), ...) {
    pmMissing <- missing(pm)
    dpmMissing <- missing(dpm)
    progressInfo <- .bfseq_extract_progress(list(...))
    progress <- progressInfo$progress
    dots <- progressInfo$dots

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
        solver <- do.call(nbf01seq., c(list(
            k1 = k1, k0 = k0, power = power, usd = usd, null = null,
            pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, type = bftype,
            target = target, nrange = nrange, looks = looks,
            timing = timing, minN = minN, by = by, strict = strict,
            integer = TRUE, search = search, details = TRUE,
            progress = progress
        ), dots))
        design <- solver$result
        if (is.null(design)) {
            msg <- "no valid sequential design could be computed within 'nrange'"
            if (!is.null(solver$error)) {
                msg <- paste0(msg, ": ", solver$error)
            }
            stop(msg, call. = FALSE)
        }
    } else {
        stopifnot(
            length(n) == 1,
            is.numeric(n),
            is.finite(n),
            n >= 2
        )
        fixedNrange <- .bfseq_fixed_schedule_range(n = n, nrange = nrange)
        schedule <- .bfseq_schedule_spec(looks = looks, timing = timing,
                                         minN = minN, by = by,
                                         nrange = fixedNrange)
        nseq <- .bfseq_schedule_n(maxN = n, schedule = schedule)
        design <- do.call(pbf01seq, c(list(
            k1 = k1, k0 = k0, se = usd/sqrt(nseq), n = nseq,
            null = null, pm = pm, psd = psd, dpm = dpm, dpsd = dpsd,
            type = bftype, strict = strict
        ), dots))
        solver <- .bfseq_fixed_solver(n = n, target = target, design = design,
                                      schedule = schedule)
        design$solver <- solver
    }

    design$sd <- sd
    design$sample.type <- type
    design
}
