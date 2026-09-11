ntbf01. <- function(k, power, null = 0, plocation = 0, pscale = 1/sqrt(2),
                    pdf = 1, type = c("two.sample", "one.sample", "paired"),
                    alternative = c("two.sided", "less", "greater"),
                    dpm = plocation, dpsd = pscale, lower.tail = TRUE,
                    integer = TRUE, nrange = c(2, 10^4),
                    ratio = 1, drange = "adaptive",
                    tail.eps = .bfpwr_defaults$tail.eps,
                    tail.nquad = .bfpwr_defaults$tail.nquad, ...) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

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

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail),

        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1],
        nrange[1] > 1,

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail),

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer),

        length(ratio) == 1,
        is.numeric(ratio),
        is.finite(ratio),
        ratio > 0,

        .tbf01_valid_drange(drange),

        length(tail.eps) == 1,
        is.numeric(tail.eps),
        is.finite(tail.eps),
        tail.eps > 0,
        tail.eps < 0.5,

        .tbf01_valid_tail_nquad(tail.nquad)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)

    nrangeSearch <- nrange
    if (type == "two.sample") {
        minN <- max(2, floor(1/ratio) + 1)
        nrangeSearch[1] <- max(nrangeSearch[1], minN)
        if (nrangeSearch[2] <= nrangeSearch[1]) {
            stop("sample size search range ('nrange') contains no valid two-sample allocation")
        }
    }

    ## define function for numerical root-finding
    rootFun <- function(n) {
        suppressWarnings({
            n2 <- if (type == "two.sample") ceiling(n*ratio) else n
            ptbf01(k = k, n = n, n1 = n, n2 = n2, null = null,
                   plocation = plocation, pscale = pscale, pdf = pdf,
                   dpm = dpm, dpsd = dpsd, type = type,
                   alternative = alternative, lower.tail = lower.tail,
                   drange = drange, tail.eps = tail.eps,
                   tail.nquad = tail.nquad, ...) - power
        })
    }
    ## Some fixed critical-value ranges are undefined at very small n; start
    ## the sample-size search at the first finite integer endpoint.
    if (is.nan(rootFun(nrangeSearch[1]))) {
        lower <- ceiling(nrangeSearch[1])
        upper <- floor(nrangeSearch[2])
        step <- 1
        previous <- lower - 1
        probe <- lower
        while (probe <= upper && is.nan(rootFun(probe))) {
            previous <- probe
            probe <- min(upper, probe + step)
            step <- step*2
            if (previous == probe) break
        }
        if (probe <= upper && !is.nan(rootFun(probe))) {
            left <- previous + 1
            right <- probe
            while (left < right) {
                mid <- floor((left + right)/2)
                if (is.nan(rootFun(mid))) left <- mid + 1
                else right <- mid
            }
            nrangeSearch[1] <- left
        }
    }

    n <- do.call(searchN, c(list(rootFun = rootFun, nrange = nrangeSearch),
                           .bfpwr_uniroot_dots(list(...))))

    if (integer) return(ceiling(n))
    else return(n)
}


#' @title Sample size calculations for \eqn{t}-test Bayes factor
#'
#' @description This function computes the required sample size to obtain a
#'     \eqn{t}-test Bayes factor (\link{tbf01}) more extreme than a threshold
#'     \code{k} with a specified target power.
#'
#' @inheritParams ptbf01
#' @inheritParams nbf01
#' @param nrange Sample size search range over which numerical search is
#'     performed. Defaults to \code{c(2, 10^4)}
#' @param ratio Allocation ratio \code{n2 / n1} for two-sample designs.
#'     Candidate group-2 sample sizes are \code{ceiling(n1 * ratio)}. Ignored
#'     for one-sample and paired designs.
#' @param tail.eps One-sided adaptive power-boundary searches stop once the
#'     remaining predictive probability in the searched tail is at most this
#'     value. If a fixed-\code{n} power evaluation reaches that cutoff before
#'     finding a boundary, the returned 0/1 tail probability has omitted mass
#'     bounded by \code{tail.eps}. Smaller values search farther. Defaults to
#'     \code{1e-6}
#' @param ... Optional numerical controls passed to \code{\link{ptbf01}} and
#'     to the \code{stats::uniroot} sample-size search. The fixed-design
#'     critical-value range \code{drange} is handled explicitly and is not
#'     passed to the sample-size search.
#'
#' @inherit nbf01 return
#'
#' @author Samuel Pawel, František Bartoš
#'
#' @seealso \link{ptbf01}, \link{powertbf01}, \link{tbf01}
#'
#' @examples
#'  ## example from Schönbrodt and Wagenmakers (2018, p.135)
#'  ntbf01(k = 1/6, power = 0.95, dpm = 0.5, dpsd = 0, alternative = "greater")
#'  ntbf01(k = 1/6, power = 0.95, dpm = 0.5, dpsd = 0.1, alternative = "greater")
#'  ntbf01(k = 6, power = 0.95, dpm = 0, dpsd = 0, alternative = "greater",
#'         lower.tail = FALSE, nrange = c(2, 10000))
#'
#' @export
ntbf01 <- Vectorize(FUN = ntbf01.,
                    vectorize.args = c("k", "power", "null", "plocation",
                                       "pscale", "pdf", "type", "alternative",
                                       "dpm", "dpsd", "lower.tail", "integer",
                                       "ratio"))
