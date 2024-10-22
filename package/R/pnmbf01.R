pnmbf01. <- function(k, n, usd, null = 0, psd, dpm, dpsd, lower.tail = TRUE) {
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

        length(usd) == 1,
        is.numeric(usd),
        is.finite(usd),
        0 < usd,

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        0 < psd,

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
    v <- usd^2/n + dpsd^2

    ## compute power with closed-form formula
    Y <- (2*lamW::lambertW0(x = (1 + n*psd^2/usd^2)^1.5*sqrt(exp(1))/2/k) - 1)*
        (1 + usd^2/(n*psd^2))/(1 + n*dpsd^2/usd^2)
    A <- (dpm - null)/sqrt(v)
    if (Y < 0) {
        pow <- 1
    } else {
        pow <- stats::pnorm(-sqrt(Y) - A) + stats::pnorm(-sqrt(Y) + A)
    }

    if (lower.tail == TRUE) {
        return(pow)
    } else {
        return(1 - pow)
    }
}


#' @title Cumulative distribution function of the normal moment prior Bayes factor
#'
#' @description This function computes the probability of obtaining a normal
#'     moment prior Bayes factor (\link{nmbf01}) more extreme than a threshold
#'     \code{k} with a specified sample size.
#'
#' @inheritParams pbf01
#' @param psd Spread of the normal moment prior assigned to the parameter under
#'     the alternative in the analysis. The modes of the prior are located at
#'     \eqn{\pm\sqrt{2}\,\code{psd}}{+-sqrt(2)*\code{psd}}
#' @param dpm Mean of the normal design prior assigned to the parameter
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     parameter. Set to 0 to obtain a point prior at the design prior mean
#'
#' @inherit pbf01 details
#'
#' @inherit pbf01 return
#'
#' @author Samuel Pawel
#'
#' @seealso \link{nmbf01}, \link{nnmbf01}, \link{powernmbf01}
#'
#' @examples
#' ## point desing prior (psd = 0)
#' pnmbf01(k = 1/10, n = 200, usd = 2, null = 0, psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0)
#'
#' ## normal design prior to incorporate parameter uncertainty (psd > 0)
#' pnmbf01(k = 1/10, n = 200, usd = 2, null = 0, psd = 0.5/sqrt(2), dpm = 0.5, dpsd = 0.25)
#'
#' ## design prior is the null hypothesis (dpm = 0, dpsd = 0)
#' pnmbf01(k = 10, n = 200, usd = 2, null = 0, psd = 0.5/sqrt(2), dpm = 0, dpsd = 0,
#'         lower.tail = FALSE)
#'
#' @export
pnmbf01 <- Vectorize(FUN = pnmbf01.)


## ## look at example
## nsim <- 1000
## k <- 1/6
## psd <- 1
## dpm <- 1 # assume that design prior is centered around 1
## dpv <- 0.5 # assume that design prior variance is 1/2
## nseq <- seq(1, 100, 1)
## usd <- 2
## power <- t(sapply(X = nseq, FUN = function(n) {
##     se <- usd/sqrt(n) # standard error with unit variance 4
##     v <- se^2 + dpv # marginal variance of y under design prior
##     ysim <- rnorm(n = nsim, mean = dpm, usd = sqrt(v)) # simulate for comparison
##     bf01sim <- nmbf01(estimate = ysim, se = se, psd = psd)
##     powsim <- mean(bf01sim < k)
##     pow <- pnmbf01(k = k, n = n, usd = usd, psd = psd, dpm = dpm, dpsd = sqrt(dpv)) # exact
##     c("exact" = pow, "simulation" = powsim)
## }))

## ## verify power curve
## matplot(nseq, power, xlab = "n", ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * ")"),
##         ylim = c(0, 1), type = "s", las = 1, col = c(1, 4), lty = 1, lwd = 2)
## mcse <- sqrt(power[,2]*(1 - power[,2])/nsim)
## polygon(x = c(nseq, rev(nseq)),
##         y = c(power[,2] + mcse, rev(power[,2] - mcse)),
##         border = FALSE, col = adjustcolor(col = 4, alpha.f = 0.1),
##         lty = 2)
