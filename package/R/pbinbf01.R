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

    if (!is.na(dp)) {
        ## point design prior
        stopifnot(
            length(dp) == 1,
            is.numeric(dp),
            is.finite(dp),
            0 < dp, dp < 1
            )
        ## predictive PMF under the point design prior
        predpmf <- function(x) stats::dbinom(x = x, size = n, prob = dp)
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
        ## predictive PMF under the truncated Beta design prior
        predpmf. <- function(x) {
            exp(lchoose(n, x) + lbeta(da + x, db + n - x) - lbeta(da, db)) *
                diff(stats::pbeta(q = c(dl, du), shape1 = da + x,
                                  shape2 = db + n - x)) /
                diff(stats::pbeta(q = c(dl, du), shape1 = da, shape2 = db))
        }
        predpmf <- Vectorize(FUN = predpmf.)
    }

    ## BF as a function of the data
    logbf <- function(x) {
        binbf01(x = x, n = n, p0 = p0, type = type, a = a, b = b, log = TRUE)
    }

    ## find BF maximum
    xmax <- stats::optim(par = p0*n, fn = logbf, control = list(fnscale = -1),
                         lower = 0, upper = n, method = "Brent")
    ## is BF below threshold k for any possible data value?
    if (xmax$value <= log(k)) {
        if (lower.tail == TRUE) {
            return(1)
        } else {
            return(0)
        }
    }

    ## find BF minimum
    xmin <- stats::optim(par = p0*n, fn = logbf, lower = 0, upper = n,
                         method = "Brent")
    ## is BF above threshold k for any possible data value?
    if (xmin$value >= log(k)) {
        if (lower.tail == TRUE) {
            return(0)
        } else {
            return(1)
        }
    }

    ## find the critical value(s) with root-finding
    rootFun <- function(x) logbf(x) - log(k)

    ## point null test
    if (type == "point") {
        xcrit1 <- try(stats::uniroot(f = rootFun, lower = 0, upper = xmax$par)$root,
                      silent = TRUE)
        xcrit2 <- try(stats::uniroot(f = rootFun, lower = xmax$par, upper = n)$root,
                      silent = TRUE)

        ## ## plot critical values
        ## xseq <- seq(0, n)
        ## plot(xseq, exp(logbf(xseq)), type = "b", log = "y", xlab = "x", ylab = "BF")
        ## abline(h = k, lty = 2)
        ## abline(v = xcrit1, lty = 2)
        ## abline(v = xcrit2, lty = 2)

        ## data values for which bf01 <= k
        if (inherits(xcrit1, "try-error")) {
            xsuccess <- seq(ceiling(xcrit2), n)
        } else if (inherits(xcrit2, "try-error")) {
            xsuccess <- seq(0, floor(xcrit1))
        } else {
            xsuccess <- c(seq(0, floor(xcrit1)), seq(ceiling(xcrit2), n))
        }
    } else { ## type == "direction"
        xcrit <- stats::uniroot(f = rootFun, lower = 0, upper = n)$root

        ## ## plot critical values
        ## xseq <- seq(0, n)
        ## plot(xseq, exp(logbf(xseq)), type = "b", log = "y", xlab = "x", ylab = "BF")
        ## abline(h = k, lty = 2)
        ## abline(v = xcrit, lty = 2)

        ## data values for which bf01 <= k
        xsuccess <- seq(ceiling(xcrit), n)
    }

    ## compute probability of BF01 <= k under the design prior
    pow <- sum(predpmf(xsuccess))
    if (lower.tail == TRUE) return(pow)
    else return(1 - pow)
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
#' nseq <- seq(1, 1000, length.out = 100)
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
