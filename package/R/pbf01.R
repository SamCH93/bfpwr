pbf01. <- function(k, n, sd, null = 0, pm, psd, dpm = pm, dpsd = psd,
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

        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        0 < sd,

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

        length(pm) == 1,
        is.numeric(pm),
        is.finite(pm),

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        0 <= psd,

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        0 <= dpsd,

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail)
    )

    ## variance of the data based on the design prior
    v <- sd^2/n + dpsd^2

    ## point prior in analysis
    if (psd == 0) {
        Z <- (sd^2*log(k)/n/(null - pm) + (null + pm)/2 - dpm)/sqrt(v)
        if (sign(null - pm) >= 0) {
            tail <- TRUE
        } else {
            tail <- FALSE
        }
        pow <- stats::pnorm(q = Z, mean = 0, sd = 1, lower.tail = tail)
    } else {
        ## normal prior in the analysis
        X <- (log(1 + n*psd^2/sd^2) + (null - pm)^2/psd^2 - log(k^2))*
            (1 + sd^2/n/psd^2)*sd^2/n/v
        if (X < 0) {
            pow <- 1
        } else {
            M <- (dpm - null - sd^2/n/psd^2*(null - pm))/sqrt(v)
            pow <- stats::pnorm(q = -sqrt(X) - M) + stats::pnorm(q = -sqrt(X) + M)
        }
    }

    if (lower.tail) return(pow)
    else return(1 - pow)
}


#' @title Cumulative distribution function of the Bayes factor under normality
#'
#' @description This function computes the probability of obtaining a Bayes
#'     factor (\link{bf01}) smaller (or larger) than a threshold \code{k} with a
#'     specified sample size.
#'
#' @param k Bayes factor threshold
#' @param n Sample size
#' @param sd Standard deviation of one unit
#' @param null Parameter value under the point null hypothesis. Defaults to 0
#' @param pm Mean of the normal prior assigned to the parameter under the
#'     alternative in the analysis
#' @param psd Standard deviation of the normal prior assigned to the parameter
#'     under the alternative in the analysis. Set to 0 to obtain a point prior
#'     at the prior mean
#' @param dpm Mean of the normal design prior assigned to the parameter.
#'     Defaults to the same value as the analysis prior \code{pm}
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     parameter. Defaults to the same value as the analysis prior \code{psd}
#' @param lower.tail Logical indicating whether Pr(BF \eqn{\leq} \code{k})
#'     (\code{TRUE}) or Pr(BF \eqn{>} \code{k}) (\code{FALSE}) should be
#'     computed. Defaults to \code{TRUE}
#'
#' @return The probability that the Bayes factor is less or greater (depending
#'     on the specified \code{lower.tail}) than the specified threshold \code{k}
#'
#' @author Samuel Pawel
#'
#' @seealso \link{nbf01}, \link{powerbf01}, \link{bf01}
#'
#' @examples
#' ## point alternative (psd = 0)
#' pbf01(k = 1/10, n = 200, sd = 2, null = 0, pm = 0.5, psd = 0)
#'
#' ## normal alternative (psd > 0)
#' pbf01(k = 1/10, n = 100, sd = 2, null = 0, pm = 0.5, psd = 2)
#'
#' ## design prior is the null hypothesis (dpm = 0, dpsd = 0)
#' pbf01(k = 10, n = 1000, sd = 2, null = 0, pm = 0.3, psd = 2, dpm = 0, dpsd = 0, lower.tail = FALSE)
#'
#' ## draw a power curve
#' nseq <- round(exp(seq(log(10), log(10000), length.out = 100)))
#' plot(nseq, pbf01(k = 1/10, n = nseq, sd = 2, null = 0, pm = 0.3, psd = 0), type = "l",
#'      xlab = "n", ylab = bquote("Pr(BF"["01"] <= 1/10 * ")"), ylim = c(0, 1),
#'      log = "x", las = 1)
#'
#' @export
pbf01 <- Vectorize(FUN = pbf01.)

## ## verify with simulation that CDF correct
## set.seed(44)
## sd <- 2
## n <- 100
## pm <- 1
## psd <- 0
## dpm <- 2
## dpsd <- 3
## null <- 0.5
## estsim <- stats::rnorm(n = 100000, mean = dpm, sd = sqrt(dpsd^2 + sd^2/n))
## bf01sim <- bf01(estimate = estsim, se = sd/sqrt(n), null = null, pm = pm, psd = psd)
## kseq <- c(1/1000, 1/100, 1/10, 1, 3, 10)
## ## simulation
## round(sapply(X = kseq, FUN = function(k) mean(bf01sim <= k)), 4)
## ## theoretical
## round(sapply(X = kseq, FUN = function(k) {
##     pbf01(k = k, n = n, sd = sd, null = null, pm = pm, psd = psd, dpm = dpm,
##           dpsd = dpsd)
## }), 4)
