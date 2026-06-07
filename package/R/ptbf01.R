ptbf01. <- function(k, n, n1 = n, n2 = n, null = 0, plocation = 0,
                    pscale = 1/sqrt(2), pdf = 1, dpm = plocation, dpsd = pscale,
                    type = c("two.sample", "one.sample", "paired"),
                    alternative = c("two.sided", "less", "greater"),
                    lower.tail = TRUE, drange = "adaptive", ...) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(n1) == 1,
        is.numeric(n1),
        is.finite(n1),
        1 < n1,

        length(n2) == 1,
        is.numeric(n2),
        is.finite(n2),
        1 < n2,

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

        (is.numeric(drange) && length(drange) == 2 && all(is.finite(drange)) &&
         drange[2] > drange[1]) || (is.character(drange) && length(drange) == 1 &&
                                    !is.na(drange) && drange == "adaptive")
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    if (type != "two.sample") {
        if (n1 != n2) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
            n2 <- n1
        }
    }

    ## determine df and effective sample size
    if (type == "two.sample") {
        df <- n1 + n2 - 2
        neff <- 1/(1/n1 + 1/n2)
    } else {
        df <- n1 - 1
        neff <- n1
    }

    ## determine effect estimate region where BF < k for specified sample size
    dots <- list(...)
    searchDots <- .bfpwr_integrate_dots(dots = dots,
                                        rel.tol.default = 1e-2)
    rootDots <- .bfpwr_uniroot_dots(dots = dots)
    se <- 1/sqrt(neff) # standard error of SMD assuming variance is known
    estsd <- sqrt(se^2 + dpsd^2) # standard deviation of SMD under design prior
    rootFun <- function(est) {
        ## tbf01() tests against zero, so shift the analysis prior by null.
        tbf01(t = (est - null)/se, n1 = n1, n2 = n2,
              plocation = plocation - null, pscale = pscale, pdf = pdf,
              type = type, alternative = alternative, log = TRUE) - log(k)
    }
    rootFunSearch <- function(est) {
        do.call(tbf01, c(list(
            t = (est - null)/se, n1 = n1, n2 = n2,
            plocation = plocation - null, pscale = pscale, pdf = pdf,
            type = type, alternative = alternative, log = TRUE
        ), searchDots)) - log(k)
    }
    region <- .tbf01_prior_region(plocation = plocation - null,
                                  pscale = pscale, pdf = pdf,
                                  alternative = alternative)
    ## For boundary bracketing, use the fast direct integral as a scout. Any
    ## candidate root is certified against the stable BF path before use.
    rootFunFast <- function(est) {
        do.call(.tbf01_log_fast, c(list(
            t = (est - null)/se, df = df, neff = neff,
            plocation = plocation - null, pscale = pscale,
            pdf = pdf, region = region
        ), searchDots)) - log(k)
    }

    if (alternative == "two.sided") {
        if (k > 1) {
            ## check whether BF > k > 1 is achievable for given sample size
            ## maximum BF is obtained when est = null
            maxBF <- exp(rootFun(null) + log(k))
            if (is.nan(maxBF)) return(NaN)
            if (maxBF < k) {
                if (lower.tail == FALSE) {
                    return(0)
                } else {
                    return(1)
                }
            }
        }

        ## guess search range based on search range from z-test BF
        ## TODO improve robustness of adaptive strategy
        if (!is.numeric(drange) && drange == "adaptive") {
            suppressWarnings({
                X <- (log(1 + neff*pscale^2) + (null - plocation)^2/pscale^2 - 2*log(k))*
                    (1 + 1/neff/pscale^2)/neff
                sqrtX <- sqrt(X)
                if (is.nan(sqrtX)) sqrtX <- 0.3
                M <- (-null - (null - plocation)/neff/pscale^2)
                critz1 <- -sqrtX - M # limit 1 z-test BF
                critz2 <- sqrtX - M # limit 2 z-test BF
            })
            meancritz <- (critz1 + critz2)/2
            if (critz1 < critz2) {
                searchIntLow <- c(critz1 - 0.01, meancritz)
                searchIntUp <- c(meancritz, critz2 + 0.01)
            } else {
                searchIntLow <- c(critz2 - 0.01, meancritz)
                searchIntUp <- c(meancritz, critz1 + 0.01)
            }
        } else {
            meant <- mean(drange)
            searchIntLow <- c(drange[1], meant)
            searchIntUp <- c(meant, drange[2])
        }
        ## search for critical values. The lower and upper roots have opposite
        ## crossing directions, so use directional interval extension.
        upper <- try(stats::uniroot(f = rootFun, interval = searchIntUp,
                                    extendInt = "downX", ...)$root,
                     silent = TRUE)
        lower <- try(stats::uniroot(f = rootFun, interval = searchIntLow,
                                    extendInt = "upX", ...)$root,
                     silent = TRUE)

        ## compute power
        if (inherits(upper, "try-error")) {
            logpowup <- -Inf
            uperr <- TRUE
        } else {
            logpowup <- stats::pnorm(q = upper, mean = dpm, sd = estsd,
                                      lower.tail = FALSE, log.p = TRUE)
            uperr <- FALSE
        }
        if (inherits(lower, "try-error")) {
            logpowlow <- -Inf
            lowerr <- TRUE
        } else {
            logpowlow <- stats::pnorm(q = lower, mean = dpm, sd = estsd,
                                       lower.tail = TRUE, log.p = TRUE)
            lowerr <- FALSE
        }
        if ((uperr == TRUE) && (lowerr == TRUE)) {
            warning("Numerical problems finding critical value")
            logpow <- NaN
            logcomp <- NaN
        } else {
            ## BF01 <= k is the union of the two normal tails; the complement
            ## is the interval between any roots that were found.
            logpow <- min(0, .bfpwr_logspace_sum(c(logpowup, logpowlow)))
            if (!uperr && !lowerr) {
                logcomp <- .bfpwr_lpnorm_interval(lower = lower, upper = upper,
                                                  mean = dpm, sd = estsd)
            } else if (!uperr) {
                logcomp <- stats::pnorm(q = upper, mean = dpm, sd = estsd,
                                         lower.tail = TRUE, log.p = TRUE)
            } else {
                logcomp <- stats::pnorm(q = lower, mean = dpm, sd = estsd,
                                         lower.tail = FALSE, log.p = TRUE)
            }
        }
    } else {
        ## one-sided alternatives
        if (!is.numeric(drange) && drange == "adaptive") {
            searchLimit <- 256
            f0 <- .bfpwr_root_value(f = rootFunSearch, x = null)
            search <- do.call(.bfpwr_one_sided_adaptive_root, c(list(
                certify_fun = rootFunSearch, scout_fun = rootFunFast,
                alternative = alternative, origin = null, step_scale = se,
                try_opposite = FALSE, search_limit = searchLimit
            ), rootDots))
            crit <- search$root
            if (inherits(crit, "try-error") && search$search_limit_reached) {
                crit <- structure("adaptive search limit reached",
                                  class = c("bfpwr_ptbf01_search_limit",
                                            "try-error"))
            }
        } else {
            crit <- try(stats::uniroot(f = rootFun,
                                       interval = drange,
                                       extendInt = "no", ...)$root,
                        silent = TRUE)
        }
        if (inherits(crit, "try-error")) {
            if (!is.numeric(drange) && drange == "adaptive" &&
                exists("f0", inherits = FALSE) && is.finite(f0) &&
                inherits(crit, "bfpwr_ptbf01_search_limit")) {
                warning(paste0(
                    "Adaptive t power-boundary search reached |t| <= ",
                    searchLimit,
                    " without bracketing BF01 = k; returning the ",
                    "boundary-free probability implied by the search. Pass ",
                    "a wider numeric 'drange' interval to search for exact ",
                    "bounds beyond this limit."
                ))
                ## No crossing was found within the finite scan. The sign at
                ## the null determines whether all searched values are successes.
                if (f0 < 0) {
                    logpow <- 0
                    logcomp <- -Inf
                } else {
                    logpow <- -Inf
                    logcomp <- 0
                }
            } else {
                warning("Numerical problems finding critical value")
                logpow <- NaN
                logcomp <- NaN
            }
        } else {
            if (alternative == "greater") {
                logpow <- stats::pnorm(q = crit, mean = dpm, sd = estsd,
                                        lower.tail = FALSE, log.p = TRUE)
                logcomp <- stats::pnorm(q = crit, mean = dpm, sd = estsd,
                                         lower.tail = TRUE, log.p = TRUE)
            } else {
                logpow <- stats::pnorm(q = crit, mean = dpm, sd = estsd,
                                        lower.tail = TRUE, log.p = TRUE)
                logcomp <- stats::pnorm(q = crit, mean = dpm, sd = estsd,
                                         lower.tail = FALSE, log.p = TRUE)
            }
        }

    }

    if (lower.tail == TRUE) return(exp(logpow))
    else return(exp(logcomp))
}


#' @title Cumulative distribution function of the t-test Bayes factor
#'
#' @description This function computes the probability of obtaining a
#'     \eqn{t}-test Bayes factor (\link{tbf01}) more extreme than a threshold
#'     \code{k} with a specified sample size.
#'
#' @inheritParams tbf01
#' @inheritParams pbf01
#' @param null Standardized mean difference under the point null hypothesis.
#'     Defaults to \code{0}
#' @param alternative Direction of the test. Can be either \code{"two.sided"}
#'     (default), \code{"less"}, or \code{"greater"}. The latter two truncate
#'     the analysis prior to negative and positive effects, respectively. If set
#'     to \code{"less"} or \code{"greater"}, the power is only computed based on
#'     data with effect estimates in the direction of the alternative
#' @param dpm Mean of the normal design prior assigned to the standardized mean
#'     difference. Defaults to the analysis prior location
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     standardized mean difference. Set to \code{0} to obtain a point prior at
#'     the design prior mean. Defaults to the analysis prior scale
#' @param drange Standardized mean difference search range over which the
#'     critical values are searched for. Can be either set to a numerical range
#'     or to \code{"adaptive"} (default) which determines the range in an
#'     adaptive way from the other input parameters
#' @param ... Optional numerical controls. For numeric ranges and two-sided
#'     adaptive searches, arguments are passed to \code{stats::uniroot}. In
#'     adaptive one-sided searches, \code{subdivisions}, \code{rel.tol},
#'     \code{abs.tol}, \code{stop.on.error}, and \code{keep.xy} are used for BF
#'     integration, while \code{tol}, \code{maxiter}, \code{trace}, and
#'     \code{check.conv} are passed to \code{stats::uniroot}.
#'
#' @inherit pbf01 return
#'
#' @author Samuel Pawel
#'
#' @seealso \link{tbf01}, \link{ntbf01}, \link{powertbf01}
#'
#' @examples
#' ## example from Schönbrodt and Wagenmakers (2018, p. 135)
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0, alternative = "greater")
#' ptbf01(k = 6, n = 146, dpm = 0, dpsd = 0, alternative = "greater",
#'        lower.tail = FALSE)
#'
#' ## two-sided
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0)
#' ptbf01(k = 6, n = 146, dpm = 0, dpsd = 0, lower.tail = FALSE)
#'
#' ## one-sample test
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0, alternative = "greater",
#'        type = "one.sample")
#'
#' @export
ptbf01 <- Vectorize(FUN = ptbf01.,
                    vectorize.args = c("k", "n", "n1", "n2", "null",
                                       "plocation", "pscale", "pdf", "type",
                                       "alternative", "dpm", "dpsd", "type",
                                       "lower.tail"))

## ## verify with simulation
## set.seed(10000)
## nsim <- 10000
## dpm <- 0.5
## dpsd <- 0.1
## n1 <- 30
## n2 <- 120
## r <- 1/sqrt(2)
## plocation <- 0.3
## bfs <- replicate(n = nsim, expr = {
##     if (dpsd == 0) m <- dpm
##     else m <- rnorm(n = 1, mean = dpm, sd = dpsd)
##     x <- rnorm(n = n1, mean = 0, sd = 1)
##     y <- rnorm(n = n2, mean = m, sd = 1)
##     t <- t.test(y, x)$statistic
##     ## se <- sqrt(1/n1 + 1/n2)
##     ## est <- rnorm(n = 1, mean = dpm, sd = sqrt(se^2 + dpsd^2))
##     ## t <- est/se
##     tbf01(t = t, n1 = n1, n2 = n2, pscale = r, plocation = plocation,
##           type ="two.sample", alternative = "two.sided")
## })
## mean(bfs < 1/10)
## ptbf01(k = 1/10, n1 = n1, n2 = n2, pscale = r, plocation = plocation, dpm = dpm,
##        dpsd = dpsd, type = "two.sample", alternative = "two.sided", drange = "adaptive")
