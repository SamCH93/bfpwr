## function to compute Pr(binbf01 <= k)
pbinbf01. <- function(k, n, p0 = 0.5, type = c("point", "direction"), a = 1,
                      b = 1, dp = NA, da = a, db = b, dl = 0, du = 1,
                      lower.tail = TRUE) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(n) == 1,
        is.numeric(n),
        is.finite(n),
        0 < n,

        length(p0) == 1,
        is.numeric(p0),
        is.finite(p0),
        0 < p0, p0 < 1,

        length(a) == 1,
        is.numeric(a),
        is.finite(a),
        0 < a,

        length(b) == 1,
        is.numeric(b),
        is.finite(b),
        0 < b,

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail)
    )
    type <- match.arg(arg = type)

    ## n has to be an integer
    n <- ceiling(n)

    pointDesign <- !is.na(dp)
    if (pointDesign) {
        ## point design prior
        stopifnot(
            length(dp) == 1,
            is.numeric(dp),
            is.finite(dp),
            0 < dp, dp < 1
            )
        ## predictive PMF under the point design prior
        predlogpmf <- function(x) {
            stats::dbinom(x = x, size = n, prob = dp, log = TRUE)
        }
        predlogcdf <- function(q) {
            if (q < 0) {
                return(-Inf)
            }
            if (q >= n) {
                return(0)
            }
            stats::pbinom(q = q, size = n, prob = dp, log.p = TRUE)
        }
        predlogsf <- function(q) {
            ## P(X >= q), on the log scale.
            if (q <= 0) {
                return(0)
            }
            if (q > n) {
                return(-Inf)
            }
            stats::pbinom(q = q - 1, size = n, prob = dp,
                          lower.tail = FALSE, log.p = TRUE)
        }
    } else {
        ## Beta design prior
        stopifnot(
            length(da) == 1,
            is.numeric(da),
            is.finite(da),
            0 < da,

            length(db) == 1,
            is.numeric(db),
            is.finite(db),
            0 < db,

            length(dl) == 1,
            is.numeric(dl),
            is.finite(dl),
            0 <= dl, dl < 1,

            length(du) == 1,
            is.numeric(du),
            is.finite(du),
            dl < du, du <= 1
        )
        ## predictive PMF under the truncated Beta design prior. Work on the
        ## log scale because extreme truncation intervals can make the beta
        ## normalizing constants very small. The common full and one-sided
        ## truncation intervals can be evaluated vectorwise; only interior
        ## intervals need the scalar stable interval helper.
        log_norm_const <- .bfpwr_lpbeta_interval(lower = dl, upper = du,
                                                 shape1 = da, shape2 = db)
        beta_bin_logpmf <- function(x) {
            lchoose(n, x) + lbeta(da + x, db + n - x) - lbeta(da, db)
        }
        if (dl <= 0 && du >= 1) {
            predlogpmf <- beta_bin_logpmf
        } else {
            log_interval <- if (dl <= 0) {
                function(x) {
                    stats::pbeta(q = du, shape1 = da + x,
                                 shape2 = db + n - x, log.p = TRUE)
                }
            } else if (du >= 1) {
                function(x) {
                    stats::pbeta(q = dl, shape1 = da + x,
                                 shape2 = db + n - x, lower.tail = FALSE,
                                 log.p = TRUE)
                }
            } else {
                function(x) {
                    vapply(
                        X = x,
                        FUN.VALUE = numeric(1),
                        FUN = function(xi) {
                            .bfpwr_lpbeta_interval(lower = dl, upper = du,
                                                    shape1 = da + xi,
                                                    shape2 = db + n - xi)
                        }
                    )
                }
            }
            predlogpmf <- function(x) {
                beta_bin_logpmf(x) + log_interval(x) - log_norm_const
            }
        }
    }

    ## BF as a function of the data
    logbf <- function(x) {
        binbf01(x = x, n = n, p0 = p0, type = type, a = a, b = b, log = TRUE)
    }

    ## The data are integer counts, so determine the successful counts on that
    ## grid. Continuous optimization followed by floor/ceiling can lose a
    ## boundary count when the root is numerically just to one side of an
    ## integer.
    logk <- log(k)
    belowThreshold <- function(x) {
        logbf(x) <= logk
    }
    firstBelow <- function(lower, upper) {
        ## The predicate is FALSE then TRUE on this interval.
        while (lower < upper) {
            midpoint <- floor((lower + upper)/2)
            if (belowThreshold(midpoint)) {
                upper <- midpoint
            } else {
                lower <- midpoint + 1
            }
        }
        lower
    }
    lastBelow <- function(lower, upper) {
        ## The predicate is TRUE then FALSE on this interval.
        while (lower < upper) {
            midpoint <- ceiling((lower + upper)/2)
            if (belowThreshold(midpoint)) {
                lower <- midpoint
            } else {
                upper <- midpoint - 1
            }
        }
        lower
    }

    logpow <- NULL
    xsuccess <- integer(0)
    if (type == "direction") {
        ## The directional BF01 decreases with the number of successes.
        if (belowThreshold(0)) {
            if (pointDesign) {
                logpow <- 0
            } else {
                xsuccess <- 0:n
            }
        } else if (belowThreshold(n)) {
            xcrit <- firstBelow(0, n)
            if (pointDesign) {
                logpow <- predlogsf(xcrit)
            } else {
                xsuccess <- xcrit:n
            }
        }
    } else { ## type == "point"
        ## The point-null log BF is concave on the integer grid. Its increment
        ## from x to x + 1 changes sign at the value below, so the maximum is
        ## one of the adjacent integer counts. Checking a small neighborhood
        ## also protects the split against floating-point rounding at a flat
        ## maximum.
        turningPoint <- p0*(b + n - 1) - (1 - p0)*a
        maximumCandidates <- unique(pmax(
            0, pmin(n, floor(turningPoint) + (-1:2))
        ))
        xmax <- maximumCandidates[which.max(logbf(maximumCandidates))]

        if (belowThreshold(xmax)) {
            if (pointDesign) {
                logpow <- 0
            } else {
                xsuccess <- 0:n
            }
        } else {
            leftCrit <- NULL
            rightCrit <- NULL
            if (belowThreshold(0)) {
                leftCrit <- lastBelow(0, xmax)
            }
            if (belowThreshold(n)) {
                rightCrit <- firstBelow(xmax, n)
            }

            if (pointDesign) {
                leftProbability <- if (is.null(leftCrit)) {
                    -Inf
                } else {
                    predlogcdf(leftCrit)
                }
                rightProbability <- if (is.null(rightCrit)) {
                    -Inf
                } else {
                    predlogsf(rightCrit)
                }
                logpow <- .bfpwr_logspace_sum(c(leftProbability,
                                                 rightProbability))
            } else {
                if (!is.null(leftCrit)) {
                    xsuccess <- c(xsuccess, 0:leftCrit)
                }
                if (!is.null(rightCrit)) {
                    xsuccess <- c(xsuccess, rightCrit:n)
                }
            }
        }
    }

    ## compute probability of BF01 <= k under the design prior
    ## Sum the selected predictive probabilities on the log scale; xsuccess
    ## may be a far tail set for stringent thresholds.
    if (is.null(logpow)) {
        logpow <- if (length(xsuccess) == 0) {
            -Inf
        } else {
            .bfpwr_logspace_sum(predlogpmf(xsuccess))
        }
    }
    logpow <- min(0, logpow)
    if (lower.tail == TRUE) {
        return(exp(logpow))
    } else {
        return(exp(.bfpwr_logspace_sub(0, logpow)))
    }
}


#' @title Cumulative distribution function of the binomial Bayes factor
#'
#' @description This function computes the probability of obtaining a binomial
#'     Bayes factor (\link{binbf01}) more extreme than a threshold \code{k} with
#'     a specified sample size.
#'
#' @inheritParams binbf01
#' @param k Bayes factor threshold
#' @param a Number of successes parameter of the beta analysis prior
#'     distribution. Defaults to \code{1}
#' @param b Number of failures parameter of the beta analysis prior
#'     distribution. Defaults to \code{1}
#' @param dp Fixed binomial proportion assumed for the power calculation. Set to
#'     \code{NA} to use a truncated beta design prior instead (specified via the
#'     \code{da}, \code{db}, \code{dl}, and \code{du} arguments). Defaults to
#'     \code{NA}
#' @param da Number of successes parameter of the truncated beta design prior
#'     distribution. Is only taken into account if \code{dp = NA}. Defaults to
#'     the same value \code{a} as specified for the analysis prior
#' @param db Number of failures parameter of the truncated beta design prior
#'     distribution. Is only taken into account if \code{dp = NA}. Defaults to
#'     the same value \code{b} as specified for the analysis prior
#' @param dl Lower truncation limit of of the truncated beta design prior
#'     distribution. Is only taken into account if \code{dp = NA}. Defaults to
#'     \code{0}
#' @param du Upper truncation limit of of the truncated beta design prior
#'     distribution. Is only taken into account if \code{dp = NA}. Defaults to
#'     \code{1}
#' @param lower.tail Logical indicating whether Pr(\eqn{\mathrm{BF}_{01}}{BF01}
#'     \eqn{\leq}{<=} \code{k}) (\code{TRUE}) or Pr(\eqn{\mathrm{BF}_{01}}{BF01}
#'     \eqn{>} \code{k}) (\code{FALSE}) should be computed. Defaults to
#'     \code{TRUE}
#'
#' @return The probability that the Bayes factor is less or greater (depending
#'     on the specified \code{lower.tail}) than the specified threshold \code{k}
#'
#' @author Samuel Pawel
#'
#' @seealso \link{binbf01}, \link{nbinbf01}
#'
#' @examples
#' ## compute probability that BF > 10 under the point null
#' a <- 1
#' b <- 1
#' p0 <- 3/4
#' k <- 10
#' nseq <- seq(1, 1000, length.out = 50)
#' powH0 <- pbinbf01(k = k, n = nseq, p0 = p0, type = "point", a = a, b = b,
#'                   dp = p0, lower.tail = FALSE)
#' plot(nseq, powH0, type = "s", xlab = "n", ylab = "Power")
#'
#' ## compare to normal approximation
#' pm <- a/(a + b) # prior mean under H1
#' psd <- sqrt(a*b/(a + b)^2/(a + b + 1)) # prior standard deviation under H1
#' pownormH0 <- pbf01(k = k, n = nseq, usd = sqrt(p0*(1 - p0)), null = p0,
#'                    pm = pm, psd = psd, dpm = p0, dpsd = 0, lower.tail = FALSE)
#' lines(nseq, pownormH0, type = "s", col = 2)
#' legend("right", legend = c("Exact", "Normal approximation"), lty = 1,
#'        col = c(1, 2))
#'
#' ## compute probability that BF < 1/10 under the p|H1 ~ Beta(a, b) alternative
#' a <- 10
#' b <- 5
#' p0 <- 3/4
#' k <- 1/10
#' powH1 <- pbinbf01(k = k, n = nseq, p0 = p0, type = "point", a = a, b = b,
#'                   da = a, db = b, dl = 0, du = 1)
#' plot(nseq, powH1, type = "s", xlab = "n", ylab = "Power")
#'
#' ## compare to normal approximation
#' pm <- a/(a + b) # prior mean under H1
#' psd <- sqrt(a*b/(a + b)^2/(a + b + 1)) # prior standard deviation under H1
#' pownormH1 <- pbf01(k = k, n = nseq, usd = sqrt(pm*(1 - pm)), null = p0,
#'                    pm = pm, psd = psd, dpm = pm, dpsd = psd)
#' lines(nseq, pownormH1, type = "s", col = 2)
#' legend("right", legend = c("Exact", "Normal approximation"), lty = 1,
#'        col = c(1, 2))
#'
#' ## probability that directional BF <= 1/10 under uniform [3/4, 1] design prior
#' pow <- pbinbf01(k = 1/10, n = nseq, p0 = 3/4, type = "direction", a = 1, b = 1,
#'                 da = 1, db = 1, dl = 3/4, du = 1)
#' plot(nseq, pow, type = "s", xlab = "n", ylab = "Power")
#' @export
pbinbf01 <- Vectorize(FUN = pbinbf01.)
