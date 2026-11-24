#' @title Numerical search for sample size
#'
#' @description Determine the sample size such that a given root function is
#'     zero
#'
#' @param rootFun Function that returns power - target power as a function of
#'     the sample size n
#' @param nrange Sample size search range over which numerical search is
#'     performed
#' @param peak Logical indicating whether to check for an interior power
#'     maximum if both endpoints fall below the target. Power need not be
#'     monotone in n; for example, the probability of misleading evidence
#'     first rises and then falls. The crossing on the rising side is returned.
#' @param ... Other arguments passed to \code{stats::uniroot}. The sample-size
#'     root tolerance defaults to \code{tol = 1e-8} and can be overridden.
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @noRd
#'
#' @keywords internal

searchN <- function(rootFun, nrange, tol = getOption("bfpwr.tol", 1e-8),
                    peak = FALSE, ...) {
    ## check boundaries of sample size search range
    lower <- rootFun(nrange[1])
    upper <- rootFun(nrange[2])
    if (peak && is.finite(lower) && lower <= 0 &&
        is.finite(upper) && upper < 0) {
        ## Power need not be increasing in N. Locate its interior peak on
        ## log(N), then bracket the rising crossing against the lower end.
        maximum <- stats::optimize(function(logn) rootFun(exp(logn)),
                                    interval = log(nrange), maximum = TRUE,
                                    tol = tol)
        if (maximum$objective >= 0) {
            nrange[2] <- exp(maximum$maximum)
            upper <- maximum$objective
        }
    }
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
        ## uniroot otherwise evaluates these endpoints again. We already need
        ## them for the range checks, so pass them through directly.
        res <- try(stats::uniroot(f = rootFun, interval = nrange,
                                  f.lower = lower, f.upper = upper,
                                  tol = tol, ...)$root)
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
#' @details This rule is intended for power functions of discrete data, such
#'     as the binomial power in \code{nbinbf01()}. These increase with n
#'     overall but oscillate in a saw-tooth pattern, so power can reach the
#'     target at n and fall below it at n + 1 (Chernick and Liu, 2002,
#'     \doi{10.1198/000313002317572835}). A design based on the first crossing
#'     then loses the target power if one more observation is collected, so
#'     the function deliberately does not return the first crossing. This is
#'     different from a genuinely non-monotone power function of continuous
#'     data (for example, the probability of misleading evidence), for which
#'     \code{searchN(peak = TRUE)} returns the rising crossing.
#'
#'     \code{searchN()} first locates a crossing. If power falls below the
#'     target within the next \code{nextend} sample sizes, the search is
#'     repeated from the largest of these sample sizes with power below the
#'     target, at most \code{maxcycles} times. If \code{searchN()} cannot
#'     bracket the target, or power falls below the target beyond the upper
#'     bound of \code{nrange}, \code{NaN} is returned with a warning, as in the
#'     other sample-size functions.
#'
#' @param rootFun Function that returns power - target power as a function of
#'     the sample size n
#' @param nrange Sample size search range over which numerical search is
#'     performed
#' @param nextend Number of samples sizes beyond solution for which it should be
#'     verified that the power does not drop below the target
#' @param maxcycles Maximum number of cycles to check that power does not drop
#'     beyond target. Defaults to \code{5}
#' @param ... Other arguments passed to \code{stats::uniroot}. The sample-size
#'     root tolerance defaults to \code{tol = 1e-8} and can be overridden.
#'
#' @return The required sample size to achieve the specified power, or
#'     \code{NaN} with a warning if no such sample size is found within
#'     \code{nrange}
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
        ni <- searchN(rootFun = rootFun, nrange = nrangei, ...)
        if (is.nan(ni)) {
            ## searchN() has already warned why nrange does not bracket the
            ## target (e.g., power already above target at nrange[1])
            return(ni)
        }
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
                ## otherwise search again above the last sample size with
                ## power below target
                cycles <- cycles + 1
                nbelow <- nextendi[which(rootextendi < 0)]
                nlast <- nbelow[length(nbelow)]
                if (nlast >= nrange[2]) {
                    warning("power falls below the target power above the upper bound of the sample size search range ('nrange')")
                    return(NaN)
                }
                nrangei <- c(nlast, nrange[2])
            }
        }
    }
    warning(paste("Power function may still fall below target power, extend root-finding cycles"))
    return(ni)
}
