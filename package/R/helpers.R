#' @title Numerical search for sample size
#'
#' @description Determine the sample size such that a given root function is
#'     zero
#'
#' @param rootFun Function that returns power - target power as a function of
#'     the sample size n
#' @param nrange Sample size search range over which numerical search is
#'     performed
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @noRd
#'
#' @keywords internal

searchN <- function(rootFun, nrange, ...) {
    ## check boundaries of sample size search range
    lower <- rootFun(nrange[1])
    upper <- rootFun(nrange[2])
    if (is.nan(lower)) {
        warning("lower bound of sample size search range ('nrange') leads to Power = NaN")
        n <- NaN
    } else if (is.nan(upper)) {
        warning("upper bound of sample size search range ('nrange') leads to Power = NaN")
        n <- NaN
    } else if (lower > 0) {
        warning("lower bound of sample size search range ('nrange') leads to higher power than specified")
        n <- NaN
    } else if (upper < 0) {
        warning("upper bound of sample size search range ('nrange') leads to lower power than specified")
        n <- NaN
    } else {
        ## perform root-finding
        res <- try(stats::uniroot(f = rootFun, interval = nrange, ... = ...)$root)
        if (inherits(res, "try-error")) {
            warning("problems while running uniroot")
            n <- NaN
        } else {
            n <- res
        }
    }
    return(n)
}
