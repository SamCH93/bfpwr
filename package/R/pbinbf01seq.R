## function to compute critical value for which we get a directional binomial
## BF01(x) <= k for x = xcrit, xcrit + 1, ..., n
xcritbinbf01. <- function(k, n, p0 = 0.5, a = 1, b = 1) {

    ## BF as a function of the data
    logbf <- function(x) {
        binbf01(x = x, n = n, p0 = p0, type = "direction", a = a, b = b,
                log = TRUE)
    }

    ## compute minBF01 (x = n) and maxBF01 (x = 0)
    minlogBF <- logbf(x = n)
    maxlogBF <- logbf(x = 0)

    if (maxlogBF <= log(k)) {
        ## is BF below threshold k for any possible data value?
        xcrit <- 0
    } else if (minlogBF > log(k)) {
        ## is BF above threshold k for any possible data value?
        xcrit <- NaN
    } else {
        ## find the critical value with root-finding
        rootFun <- function(x) logbf(x) - log(k)
        res <- stats::uniroot(f = rootFun, lower = 0, upper = n)$root
        xcrit <- ceiling(res)
    }

    ## return critical value
    return(xcrit)
}
xcritbinbf01 <- Vectorize(FUN = xcritbinbf01.)

#' @title Power calculations for binomial Bayes factor with two-stage designs
#'
#' @description Compute probability that binomial Bayes factor (\link{binbf01})
#'     is equal or below a threshold at final stage, taking into account that
#'     the study can stop if it is above a threshold at interim.
#'
#'
#' @inheritParams pbinbf01
#' @inheritParams powerbf01
#' @param n1 Sample size at interim stage
#' @param n2 Sample size at final stage
#' @param k Bayes factor threshold for evidence for efficacy . Defaults to
#'     \code{1/10}, Jeffreys' threshold for 'strong evidence' against the null
#'     hypothesis
#' @param kf Bayes factor threshold for evidence for futility. Defaults to
#'     \code{1/k}.
#' @param type Typo of test. Currently only \code{"directional"} is implemented
#'
#' @return The probability of obtaining a Bayes factor equal or below k at the
#'     final stage taking into account potential futility stopping with BF01 >
#'     kf at interim
#'
#' @author Samuel Pawel
#'
#' @seealso \link{binbf01}
#'
#' @examples
#' ## probability of BF<1/10 under truncated beta prior
#' pbinbf01seq(n1 = 10, n2 = 20, k = 1/10, kf = 3, p0 = 0.2, a = 1, b = 1,
#'             dl = 0.2, du = 1, type = "direction") # with interim analysis
#' pbinbf01(n = 20, k = 1/10, p0 = 0.2, a = 1, b = 1, dl = 0.2, du = 1,
#'          type = "direction") # without interim analysis
#'
#' ## probability to stop at interim under Bayesian null
#' pbinbf01(n = 10, k = 3, p0 = 0.2, a = 1, b = 1, dl = 0, du = 0.2,
#'          type = "direction", lower.tail = FALSE)
#'
#' ## probability of BF<1/10 under point null
#' pbinbf01seq(n1 = 10, n2 = 20, k = 1/10, kf = 3, p0 = 0.2, a = 1, b = 1,
#'             dp = 0.2, type = "direction") # with interim analysis
#' pbinbf01(n = 20, k = 1/10, p0 = 0.2, a = 1, b = 1, dp = 0.2,
#'          type = "direction") # without interim analysis
#'
#' ## probability to stop at interim under point null
#' pbinbf01(n = 10, k = 3, p0 = 0.2, a = 1, b = 1, dp = 0.2,
#'          type = "direction", lower.tail = FALSE)
#'
#' @export
pbinbf01seq <- function(n1, n2, k = 1/10, kf = 1/k, p0 = 0.5,
                        type = c("direction"), a = 1, b = 1, dp = NA, da = a,
                        db = b, dl = 0, du = 1) {
    ## input checks
    stopifnot(
        length(n1) == 1,
        is.numeric(n1),
        is.finite(n1),
        0 < n1,

        length(n2) == 1,
        is.numeric(n2),
        is.finite(n2),
        0 < n2, n1 < n2
    )
    type <- match.arg(arg = type)

    ## n1 and n2 have to be integers
    n1 <- ceiling(n1)
    n2 <- ceiling(n2)

    if (!is.na(dp)) {
        ## point design prior
        stopifnot(
            length(dp) == 1,
            is.numeric(dp),
            is.finite(dp),
            0 < dp, dp < 1
            )
        ## predictive PMF under the point design prior
        predpmf <- function(x1, x2, n1, n2) {
            stats::dbinom(x = x1, size = n1, prob = dp) *
                stats::dbinom(x = x2, size = n2, prob = dp)
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
        ## predictive PMF under the truncated Beta design prior
        predpmf. <- function(x1, x2, n1, n2) {
            exp(lchoose(n1, x1) + lchoose(n2, x2) +
                lbeta(da + x1 + x2, db + n1 + n2 - x1 - x2) -
                lbeta(da, db)) *
                diff(stats::pbeta(q = c(dl, du), shape1 = da + x1 + x2,
                                  shape2 = db + n1 + n2 - x1 - x2)) /
                diff(stats::pbeta(q = c(dl, du), shape1 = da, shape2 = db))
        }
        predpmf <- Vectorize(FUN = predpmf.)
    }


    ## probability to stop at interim Pr(interim stage BF01 >= kf)
    pinterim <- pbinbf01(k = kf, n = n1, p0 = p0, type = type, a = a, b = b,
                         dp = dp, da = da, db = db, dl = dl, du = du,
                         lower.tail = FALSE)

    ## probability of Pr(final stage BF01 <= k) when no interim analysis
    pfinal <- pbinbf01(k = k, n = n2, p0 = p0, type = type, a = a, b = b,
                       dp = dp, da = da, db = db, dl = dl, du = du,
                       lower.tail = TRUE)

    if (type == "direction") {
        ## determine critical values where BF01 <= k at interim/final stages for
        ## efficacy/futility thresholds
        xcrit1f <- xcritbinbf01(k = kf, n = n1, p0 = p0, a = a, b = b)
        xcrit2 <- xcritbinbf01(k = k, n = n2, p0 = p0, a = a, b = b)

        ## ## should be the same as pfinal
        ## sum(predpmf(x = seq(xcrit2, n2), n = n2))

        if (xcrit1f == 0) {
            ## BF01 is always below futility threshold kf at interim
            pcorrect <- 0
        } else {
            ## probability to stop for futility at interim and for efficacy at final
            pcorrect <- sum(sapply(X = seq(0, xcrit1f - 1), FUN = function(x1i) {
                if (xcrit2 - x1i > n2 - n1) {
                    ## final stage BF01 <= k is not possible anymore when x1 = x1i
                    return(0)
                } else {
                    x2j <- seq(xcrit2 - x1i, n2 - n1)
                    sum(predpmf(x1 = x1i, x2 = x2j, n1 = n1, n2 = n2 - n1))
                }
            }))
        }
    } else {
        ## TODO implement for point null hypothesis BF
    }

    ## adjust unconditional probability
    padj <- pfinal - pcorrect

    return(padj)

    ## ## return object
    ## structure(list(power = padj, n1 = n1, n2 = n2, p0 = p0, type = type, a = a,
    ##                b = b, dp = dp, da = da, db = db, dl = dl, du = du, k = k,
    ##                kf = kf, type = type, test = "binomial"),
    ##           class = "power.bftestseq")

}

## plotcritvals <- function(n1, n2, k, kf, a, b, p0, type = "direction") {
##     xseq1 <- seq(0, n1)
##     bf1 <- binbf01(x = xseq1, n = n1, type = type, a = a, b = b, p0 = p0)
##     xseq2 <- seq(0, n2)
##     bf2 <- binbf01(x = xseq2, n = n2, type = type, a = a, b = b, p0 = p0)
##     ## determine critical values
##     ## efficacy
##     xcrit1e <- xcritbinbf01(k = k, n = n1, p0 = p0, a = a, b = b)
##     xcrit2e <- xcritbinbf01(k = k, n = n2, p0 = p0, a = a, b = b)
##     ## futility (need to subtract one because BF needs to be greater than k)
##     xcrit1f <- xcritbinbf01(k = kf, n = n1, p0 = p0, a = a, b = b) - 1
##     if (xcrit1f < 0) xcrit1f <- NaN
##     xcrit2f <- xcritbinbf01(k = kf, n = n2, p0 = p0, a = a, b = b) - 1
##     if (xcrit2f < 0) xcrit2f <- NaN
##     ## plot BFs and critical values
##     par(mfrow = c(1, 2))
##     ## at interim stage
##     plot(xseq1, bf1, type = "b", xlab = "x", ylab = bquote("BF"["01"]),
##          main = "Interim stage", ylim = range(bf2), log = "y")
##     abline(h = 1, col = adjustcolor(col = 1, alpha.f = 0.2))
##     abline(h = c(k, kf), lty = 2)
##     axis(side = 4, at = c(kf, 1, k), las = 1)
##     axis(side = 3, at = c(xcrit1f, xcrit1e))
##     abline(v = c(xcrit1e, xcrit1f), lty = 3)
##     ## at final stage
##     plot(xseq2, bf2, type = "b",  xlab = "x", ylab = bquote("BF"["01"]),
##          main = "Final stage", ylim = range(bf2), log = "y")
##     abline(h = 1, col = adjustcolor(col = 1, alpha.f = 0.2))
##     abline(h = c(k, kf), lty = 2)
##     axis(side = 4, at = c(k, 1, kf), las = 1)
##     axis(side = 3, at = c(xcrit2f, xcrit2e))
##     abline(v = c(xcrit2e, xcrit2f), lty = 3)
## }

## ## check that calculations correct
## n1 <- 20
## n2 <- 40
## k <- 1/10
## kf <- 1 # set a low futility threshold to better see differences
## a <- 1
## b <- 1
## p0 <- 0.2
## p1 <- 0.4
## da <- a
## db <- b
## dl <- p0
## du <- 1
## plotcritvals(n1, n2, k, kf, a, b, p0, "direction")

## ## simulate BFs to compute probabilities empirically
## set.seed(4242)
## nsim <- 10^5

## ## point prior H0 sequential
## x1simpH0 <- rbinom(n = nsim, size = n1, prob = p0)
## x2simpH0 <- rbinom(n = nsim, size = n2 - n1, prob = p0)
## xsimpH0 <- x1simpH0 + x2simpH0
## bf1simpH0 <- binbf01(x = x1simpH0, n = n1, p0 = p0, type = "direction", a = a, b = b)
## bf2simpH0 <- binbf01(x = xsimpH0, n = n2, p0 = p0, type = "direction", a = a, b = b)
## mean(bf1simpH0 < kf & bf2simpH0 < k)
## pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0, a = a, b = b,
##             dp = p0, type = "direction")
## ## point prior H0 final only
## mean(bf2simpH0 < k)
## pbinbf01(n = n2, k = k, p0 = p0, a = a, b = b, dp = p0, type = "direction")

## ## point prior H1 sequential
## x1simpH1 <- rbinom(n = nsim, size = n1, prob = p1)
## x2simpH1 <- rbinom(n = nsim, size = n2 - n1, prob = p1)
## xsimpH1 <- x1simpH1 + x2simpH1
## bf1simpH1 <- binbf01(x = x1simpH1, n = n1, p0 = p0, type = "direction", a = a, b = b)
## bf2simpH1 <- binbf01(x = xsimpH1, n = n2, p0 = p0, type = "direction", a = a, b = b)
## mean(bf1simpH1 < kf & bf2simpH1 < k)
## pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0, a = a, b = b,
##             dp = p1, type = "direction")
## ## point prior H1 final only
## mean(bf2simpH1 < k)
## pbinbf01(n = n2, k = k, p0 = p0, a = a, b = b, dp = p1, type = "direction")

## ## truncateted beta prior H0 sequential
## psim <- rbeta(n = nsim*10, shape1 = a, shape2 = b)
## psimH0 <- psim[psim <= p0][1:nsim]
## x1simtbH0 <- rbinom(n = nsim, size = n1, prob = psimH0)
## x2simtbH0 <- rbinom(n = nsim, size = n2 - n1, prob = psimH0)
## xsimtbH0 <- x1simtbH0 + x2simtbH0
## bf1simtbH0 <- binbf01(x = x1simtbH0, n = n1, p0 = p0, type = "direction", a = a, b = b)
## bf2simtbH0 <- binbf01(x = xsimtbH0, n = n2, p0 = p0, type = "direction", a = a, b = b)
## mean(bf1simtbH0 < kf & bf2simtbH0 < k) # does not work yet
## pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0, a = a, b = b,
##             da = da, db = db, dl = 0, du = p0, type = "direction")
## ## truncateted beta prior H0 final only
## x3simtbH0 <- rbinom(n = nsim, size = n2, prob = psimH0)
## bfsimtbH0 <- binbf01(x = x3simtbH0, n = n2, p0 = p0, type = "direction", a = a, b = b)
## mean(bfsimtbH0 < k) # this works
## pbinbf01(n = n2, k = k, p0 = p0, a = a, b = b, da = da, db = db, dl = 0,
##          du = p0, type = "direction")

## ## truncateted beta prior H1 sequential
## psim <- rbeta(n = nsim*10, shape1 = a, shape2 = b)
## psimH1 <- psim[psim > p0][1:nsim]
## x1simtbH1 <- rbinom(n = nsim, size = n1, prob = psimH1)
## x2simtbH1 <- rbinom(n = nsim, size = n2 - n1, prob = psimH1)
## xsimtbH1 <- x1simtbH1 + x2simtbH1
## bf1simtbH1 <- binbf01(x = x1simtbH1, n = n1, p0 = p0, type = "direction", a = a, b = b)
## bf2simtbH1 <- binbf01(x = xsimtbH1, n = n2, p0 = p0, type = "direction", a = a, b = b)
## mean(bf1simtbH1 < kf & bf2simtbH1 < k) # does not work yet
## pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0, a = a, b = b,
##             da = da, db = db, dl = p0, du = 1)
## ## truncateted beta prior H1 final only
## x3simtbH1 <- rbinom(n = nsim, size = n2, prob = psimH1)
## bfsimtbH1 <- binbf01(x = x3simtbH1, n = n2, p0 = p0, type = "direction", a = a, b = b)
## mean(bfsimtbH1 < k) # this works
## pbinbf01(n = n2, k = k, p0 = p0, a = a, b = b, da = da, db = db, dl = p0,
##          du = 1, type = "direction")
