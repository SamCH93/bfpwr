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


#' @title Numerical search for sample size with oscillating power function
#'
#' @description Determine the sample size such that a given root function is
#'     zero based on an ascending but oscillating power function. After a
#'     solution is found the search terminates only if the power function does
#'     not drop below the target power for a specified range of sample sizes
#'     above the solution.
#'
#' @param rootFun Function that returns power - target power as a function of
#'     the sample size n
#' @param nrange Sample size search range over which numerical search is
#'     performed
#' @param nextend Number of samples sizes for which solution
#' @param maxcycles Maximum number of cycles. Defaults to \code{5}
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @noRd
#'
#' @keywords internal

searchNoscil <- function(rootFun, nrange, nextend = 10, maxcycles = 5, ...) {
    cycles <- 0
    nrangei <- nrange
    while (cycles <= maxcycles) {
        ni <- searchN(rootFun = rootFun, nrange = nrangei, ... = ...)
        nextendi <- seq(ni + 1, ni + 10, 1)
        rootextendi <- rootFun(nextendi)
        if (all(rootextendi >= 0)) {
            return(ni)
        } else {
            cycles <- cycles + 1
            nrangei <- c(nextendi[which(rootextendi < 0)][1], nrange[2])
        }
    }
    warning(paste("Power function may still fall below target power, extend root-fnding cycles"))
    return(ni)
}
