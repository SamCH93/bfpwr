ntbf01. <- function(k, power, null = 0, plocation = 0, pscale = 1/sqrt(2),
                    pdf = 1, type = c("two.sample", "one.sample", "paired"),
                    alternative = c("two.sided", "less", "greater"),
                    dpm = plocation, dpsd = pscale, lower.tail = TRUE,
                    integer = TRUE, nrange = c(2, 10^4), ...) {
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
        !is.na(integer)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)


    ## define function for numerical root-finding
    rootFun <- function(n) {
        suppressWarnings({
            ptbf01(k = k, n = n, null = null, plocation = plocation,
                   pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd,
                   type = type, alternative = alternative,
                   lower.tail = lower.tail, ... = ...) - power
        })
    }

    n <- searchN(rootFun = rootFun, nrange = nrange, ... = ...)

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
#'
#' @inherit nbf01 return
#'
#' @author Samuel Pawel
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
                                       "dpm", "dpsd", "lower.tail", "integer"))
