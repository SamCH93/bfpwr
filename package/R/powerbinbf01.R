#' @title Power and sample size calculations for binomial Bayes factor
#'
#' @description Compute probability that binomial Bayes factor (\link{binbf01})
#'     is smaller than a specified threshold (the power), or determine sample
#'     size to obtain a target power.
#'
#' @details This function provides a similar interface as
#'     \code{stats::power.prop.test}. For some users, the low-level functions
#'     \link{nbinbf01} (to directly compute the sample size for a fixed power)
#'     and \link{pbinbf01} (to directly compute the power for a fixed sample
#'     size) may also be useful.
#'
#' @inheritParams pbinbf01
#' @inheritParams powerbf01
#' @param n Sample size. Has to be \code{NULL} if \code{power} is specified.
#'     Defaults to \code{NULL}
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param nrange Sample size search range over which numerical search is
#'     performed (only taken into account when \code{n} is \code{NULL}).
#'     Defaults to \code{c(1, 10^3)}
#'
#' @inherit powerbf01 return
#'
#' @author Samuel Pawel
#'
#' @seealso \link{plot.power.bftest}, \link{pbinbf01}, \link{nbinbf01},
#'     \link{binbf01}
#'
#' @examples
#' ## determine sample size
#' (nres <- powerbinbf01(power = 0.8, p0 = 0.2, type = "direction", dl = 0.2))
#' \dontrun{
#' plot(nres, nlim = c(1, 250), ngrid = 250, type = "s")
#' }
#'
#' ## determine power
#' (powres <- powerbinbf01(n = 100, type = "point"))
#' \dontrun{
#' plot(powres)
#' }
#'
#' @export
powerbinbf01 <- function(n = NULL, power = NULL, k = 1/10, p0 = 0.5,
                         type = c("point", "direction"), a = 1, b = 1, dp = NA,
                         da = a, db = b, dl = 0, du = 1, nrange = c(1, 10^4)) {
    ## input checks
    if (is.null(n) && is.null(power)) {
        stop("exactly one of 'n' and 'power' must be NULL")
    }
    if (is.null(n)) {
        stopifnot(
            length(power) == 1,
            is.numeric(power),
            is.finite(power),
            0 < power, power < 1
        )
    } else {
        stopifnot(
            length(n) == 1,
            is.numeric(n),
            is.finite(n),
            0 < n
        )
    }
    ## determine sample size
    if (is.null(n)) {
        n <- nbinbf01(k = k, power = power, p0 = p0, type = type, a = a, b = b,
                      dp = dp, da = da, db = db, dl = dl, du = du,
                      nrange = nrange)
    } else {
        ## determine power
        power <- pbinbf01(k = k, n = n, p0 = p0, type = type, a = a, b = b,
                          dp = dp, da = da, db = db, dl = dl, du = du,
                          lower.tail = TRUE)
    }

    ## return object
    structure(list(n = n, power = power, p0 = p0, type = type, a = a, b = b,
                   dp = dp, da = da, db = db, dl = dl, du = du, k = k,
                   nrange = nrange, type = type, test = "binomial"),
              class = "power.bftest")

}
