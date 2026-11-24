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
#'     as the binomial power in \code{nbinbf01()}, and must be kept. Power,
#'     coverage, and type I error rates of discrete-data models can oscillate
#'     in a saw-tooth (zig-zag) pattern, so power can reach the target at n and
#'     fall below it at a slightly larger n (Chernick and Liu, 2002,
#'     \doi{10.1198/000313002317572835}). Choosing the sample size so that
#'     power does not drop below the target again is common practice for
#'     discrete data (e.g., Julious, 2023, Chapter 17,
#'     \doi{10.1201/9780429503658}), so the function deliberately does not
#'     simply return the first crossing. This is qualitatively different from
#'     a non-monotone power function of continuous data (for example, the
#'     probability of misleading evidence), for which
#'     \code{searchN(peak = TRUE)} can locate an interior maximum to bracket
#'     the rising crossing.
#'
#'     \code{searchN()} first locates a crossing, and power is checked at the
#'     crossing and the next \code{nextend} sample sizes. If every checked
#'     point meets the target, the crossing is returned. If only the crossing
#'     falls below the target, the crossing + 1 is returned (so only
#'     \code{nextend - 1} successors are checked). Otherwise, the search is
#'     restarted from the last checked point below the target, up to
#'     \code{maxcycles} times after the initial search; if all restarts are
#'     used, a warning is issued and the last crossing is returned. A
#'     \code{NaN} from \code{searchN()} is returned with its warning.
#'
#' @param rootFun Function that returns power - target power as a function of
#'     the sample size n
#' @param nrange Sample size search range over which numerical search is
#'     performed
#' @param nextend Number of samples sizes beyond solution for which it should be
#'     verified that the power does not drop below the target
#' @param maxcycles Maximum number of restarts after the initial search.
#'     Defaults to \code{5}
#' @param ... Other arguments passed to \code{stats::uniroot}. The sample-size
#'     root tolerance defaults to \code{tol = 1e-8} and can be overridden.
#'
#' @return The selected numerical candidate, or \code{NaN} with the warning
#'     from \code{searchN()} when it returns \code{NaN}
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
            ## searchN() has already warned about the failure
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
