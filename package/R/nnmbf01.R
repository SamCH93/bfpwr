nnmbf01. <- function(k, power, sd, null = 0, psd, dpm, dpsd,
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

        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        0 < sd,

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
        pnmbf01(k = k, n = n, sd = sd, null = null, psd = psd, dpm = dpm,
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
#'     Bayes factor (\link{bf01}) less or greater than a threshold \code{k} with
#'     a specified target power.
#'
#' @param k Bayes factor threshold
#' @param power Target power
#' @param sd Standard deviation of one unit
#' @param null Parameter value under the point null hypothesis. Defaults to 0
#' @param psd Spread of the normal moment prior assigned to the parameter under
#'     the alternative in the analysis. The modes of the prior are located at
#'     \eqn{\pm\sqrt{2}\,\code{psd}}{+-sqrt(2)*\code{psd}}
#' @param dpm Mean of the normal design prior assigned to the parameter
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     parameter. Set to 0 to obtain a point prior at the prior mean
#' @param nrange Sample size search range over which numerical search is
#'     performed. Defaults to \code{c(1, 10^5)}
#' @param lower.tail Logical indicating whether Pr(BF \eqn{\leq} \code{k})
#'     (\code{TRUE}) or Pr(BF \eqn{>} \code{k}) (\code{FALSE}) should be
#'     computed. Defaults to \code{TRUE}
#' @param integer Logical indicating whether only integer valued sample sizes
#'     should be returned. If \code{TRUE} the required sample size is rounded to
#'     the next larger integer. Defaults to \code{TRUE}
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @seealso \link{nmbf01}, \link{pnmbf01}, \link{powernmbf01}
#'
#' @examples
#' nnmbf01(k = 1/10, power = 0.9, sd = 1, null = 0, psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0)
#' @export
nnmbf01 <- Vectorize(FUN = nnmbf01.,
                     vectorize.args = c("k", "power", "sd", "null", "psd",
                                        "dpm", "dpsd", "integer"))
