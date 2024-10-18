nnmbf01. <- function(k, power, usd, null = 0, psd, dpm, dpsd,
                     nrange = c(1, 10^5), lower.tail = TRUE, integer = TRUE,
                     ...) {
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


    ## define function for numerical root-finding
    rootFun <- function(n) {
        pnmbf01(k = k, n = n, usd = usd, null = null, psd = psd, dpm = dpm,
                dpsd = dpsd, lower.tail = lower.tail) - power
    }

    ## determine sample size numerically
    n <- searchN(rootFun = rootFun, nrange = nrange, ... = ...)
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
#'
#' @inherit pbf01 details
#'
#' @inherit nbf01 return
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
                                        "dpm", "dpsd", "integer"))
