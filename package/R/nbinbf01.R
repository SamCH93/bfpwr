nbinbf01. <- function(k, power, p0 = 0.5, type = c("point", "direction"), a = 1,
                      b = 1, dp = NA, da = a, db = b, dl = 0, du = 1,
                      lower.tail = TRUE, integer = TRUE, nrange = c(1, 10^3),
                      ...) {
    ## input checks
    stopifnot(
        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1],
        nrange[1] > 0,

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer)
    )


    ## define function for numerical root-finding
    rootFun <- function(n) {
        suppressWarnings({
            pbinbf01(k = k, n = n, p0 = p0, type = type, a = a, b = b, dp = dp,
                     da = da, db = db, dl = dl, du = du,
                     lower.tail = lower.tail) - power
        })
    }

    n <- searchN(rootFun = rootFun, nrange = nrange, ... = ...)

    if (integer) return(ceiling(n))
    else return(n)
}


#' @title Sample size determination for binomial Bayes factor
#'
#' @description This function computes the required sample size to obtain a
#'     Bayes factor (\link{binbf01}) more extreme than a threshold \code{k} with
#'     a specified target power.
#'
#' @inheritParams pbinbf01
#' @param power Target power
#' @param nrange Sample size search range over which numerical search is
#'     performed. Defaults to \code{c(1, 10^3)}
#' @param integer Logical indicating whether only integer valued sample sizes
#'     should be returned. If \code{TRUE} the required sample size is rounded to
#'     the next larger integer. Defaults to \code{TRUE}
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbinbf01}, \link{binbf01}
#'
#' @examples
#' ## directional null testing
#' nbinbf01(k = 1/10, power = 0.9, p0 = 3/4, type = "direction", a = 1, b = 1,
#'          da = 1, db = 1, dl = 3/4, du = 1)
#' pbinbf01(k = 1/10, n = 280, p0 = 3/4, type = "direction", a = 1, b = 1,
#'          da = 1, db = 1, dl = 3/4, du = 1)
#'
#' ## point null testing
#' nbinbf01(k = 1/10, power = 0.9, p0 = 3/4, type = "point", a = 1, b = 1)
#' nbinbf01(k = 3, power = 0.8, p0 = 3/4, type = "point", dp = 3/4,
#'          lower.tail = FALSE, nrange = c(1, 10^5))
#' ## FIXME doesn't work, probably uniroot has issues with oscillations!
#' pbinbf01(k = 3, n = 9, p0 = 3/4, type = "point", dp = 3/4,
#'          lower.tail = FALSE)
#'
#' @export
nbinbf01 <- Vectorize(FUN = nbinbf01.,
                      vectorize.args = c("k", "power", "p0", "type", "a", "b",
                                         "dp", "da", "db", "dl", "du",
                                         "lower.tail", "integer"))
