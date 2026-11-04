nnmbf01. <- function(k, power, usd, null = 0, psd, dpm, dpsd,
                     nrange = c(1, 10^5), lower.tail = TRUE, integer = TRUE,
                     ...,
                     alternative = c("two.sided", "less", "greater")) {
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
        0 < psd,

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        0 <= dpsd,

        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1],

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail),

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer)
    )


    alternative <- match.arg(alternative)

    ## define function for numerical root-finding
    rootFun <- function(n) {
        pnmbf01(k = k, n = n, usd = usd, null = null, psd = psd, dpm = dpm,
                dpsd = dpsd, lower.tail = lower.tail, alternative = alternative) - power
    }

    ## determine sample size numerically
    n <- searchN(rootFun = rootFun, nrange = nrange, peak = TRUE, ...)
    if (integer) return(ceiling(n))
    else return(n)
}


#' @title Sample size determination for normal moment prior Bayes factor
#'
#' @description This function computes the required sample size to obtain a
#'     normal moment prior Bayes factor (\link{nbf01}) more extreme than a
#'     threshold \code{k} with a specified target power.
#'
#' @inheritParams pnmbf01
#' @inheritParams nbf01
#' @param alternative Direction of the alternative hypothesis, one of
#'     \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#'     For a one-sided alternative, the normal moment analysis prior is
#'     restricted to values below \code{null} (\code{"less"}) or above
#'     \code{null} (\code{"greater"}) and renormalized. The point null
#'     hypothesis is unchanged; \code{psd} is the scale before truncation.
#'     Must be specified by name: unnamed arguments after \code{integer} are
#'     passed to \code{stats::uniroot}.
#'
#' @inherit pbf01 details
#'
#' @inherit nbf01 return
#'
#' @note The function uses a numerical root search within \code{nrange} and
#'     warns when it cannot bracket the target. If both endpoints are below
#'     the target, the search looks for an interior peak. A nonmonotone power
#'     curve can cross the target more than once, and the search may return a
#'     later crossing. To find the first integer crossing, use \link{nbf01seq}
#'     with \code{type = "moment"}, \code{looks = 1}, and \code{search = "exhaustive"}.
#'     For H1 (\code{lower.tail = TRUE}) with \code{k < 1}, set
#'     \code{k1 = k, k0 = 1/k, target = "H1"}. For H0 (\code{lower.tail = FALSE})
#'     with \code{k > 1}, set \code{k1 = 1/k, k0 = k, target = "H0"}.
#'     Specify the intended \code{nrange} explicitly because the defaults differ.
#'
#' @author Samuel Pawel
#'
#' @seealso \link{nmbf01}, \link{pnmbf01}, \link{powernmbf01}
#'
#' @examples
#' nnmbf01(k = 1/10, power = 0.9, usd = 1, null = 0, psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0)
#' @export
nnmbf01 <- Vectorize(FUN = nnmbf01.,
                     vectorize.args = c("k", "power", "usd", "null", "psd",
                                        "dpm", "dpsd", "integer", "alternative"))
