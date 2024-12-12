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
#' @param nextend Number of samples sizes beyond solution for which it should be
#'     verified that the power does not drop below the target
#' @param maxcycles Maximum number of cycles to check that power does not drop
#'     beyond target. Defaults to \code{5}
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
        nextendi <- seq(ni, ni + nextend, 1)
        rootextendi <- rootFun(nextendi)
        if (all(rootextendi >= 0)) {
            ## power doesn't drop below target power
            return(ni)
        } else {
            if (all(rootextendi[-1] >= 0)) {
                ## may happen that the jump is just between ni and ni + 1, but
                ## afterwards power doesn't drop anymore below target, in this
                ## case return ni + 1
                return(ni + 1)
            } else {
                ## otherwise search again
                cycles <- cycles + 1
                nbelow <- nextendi[which(rootextendi < 0)]
                nrangei <- c(nbelow[length(nbelow)], nrange[2])
            }
        }
    }
    warning(paste("Power function may still fall below target power, extend root-finding cycles"))
    return(ni)
}
