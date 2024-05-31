ptbf01. <- function(k = 1/10, n, n1 = n, n2 = n, null = 0,
                    plocation = 0, pscale = 1/sqrt(2), pdf = 1, dpm = plocation,
                    dpsd = pscale,
                    type = c("two.sample", "one.sample", "paired"),
                    alternative = c("two.sided", "less", "greater"),
                    lower.tail = TRUE,
                    strict = FALSE) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(n1) == 1,
        is.numeric(n1),
        is.finite(n1),
        0 < n1,

        length(n2) == 1,
        is.numeric(n2),
        is.finite(n2),
        0 < n2,

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

        length(strict) == 1,
        is.logical(strict),
        !is.na(strict)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    if (type != "two.sample") {
        if (n1 != n2) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
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
    ## HACK this can probably be improved...
    se <- 1/sqrt(neff) # standard error of SMD assuming variance is known
    estsd <- sqrt(se^2 + dpsd^2) # standard deviation of SMD under design prior
    rootFun <- function(est) {
        tbf01(t = (est - null)/se, n1 = n1, n2 = n2, plocation = plocation,
              pscale = pscale, pdf = pdf, type = type,
              alternative = alternative) - k
    }
    upper <- try(stats::uniroot(f = rootFun, interval = c(0, 100))$root,
                 silent = TRUE)
    lower <- try(stats::uniroot(f = rootFun, interval = c(-100, 0))$root,
                 silent = TRUE)
    if (inherits(upper, "try-error")) {
        powup <- 0
    } else {
        powup <- stats::pnorm(q = upper, mean = dpm, sd = estsd, lower.tail = FALSE)
    }
    if (inherits(lower, "try-error")) {
        powlow <- 0
    } else {
        powlow <- stats::pnorm(q = lower, mean = dpm, sd = estsd, lower.tail = TRUE)
    }
    if (strict == TRUE) {
        if (alternative == "greater") powlow <- 0
        if (alternative == "less") powup <- 0
    }
    pow <- powup + powlow

    if (lower.tail == TRUE) return(pow)
    else return(1 - pow)
}


#' @title Cumulative distribution function of the t-test Bayes factor
#'
#' @description This function computes the probability of obtaining a JZS Bayes
#'     (\link{bf01}) smaller (or larger) than a threshold \code{k} with a
#'     specified sample size.
#'
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param n Sample size (per group)
#' @param n1 Sample size in group 1 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param n2 Sample size in group 2 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param null Standardized mean difference under the point null hypothesis.
#'     Defaults to \code{0}
#' @param plocation Analysis \eqn{t} prior location. Defaults to \code{0}
#' @param pscale Analysis \eqn{t} prior scale. Defaults to \code{1/sqrt(2)}
#' @param pscale Analysis \eqn{t} prior degrees of freedom. Defaults to \code{1}
#' @param pdf Analysis \eqn{t} prior degrees of freedom. Defaults to \code{1}
#' @param type Type of \eqn{t}-test associated with \eqn{t}-statistic. Can be
#'     \code{"two.sample"} (default), \code{"one.sample"}, or \code{"paired"}
#' @param alternative Direction of the test. Can be either \code{"two.sided"}
#'     (default), \code{"less"}, or \code{"greater"}
#' @param dpm Mean of the normal design prior assigned to the standardized mean
#'     difference. Defaults to the analysis prior location
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     standardized mean difference. Set to \code{0} to obtain a point prior at
#'     the design prior mean. Defaults to the analysis prior scale
#' @param type The type of test. One of \code{"two.sample"},
#'     \code{"one.sample"}, \code{"paired"}. Defaults to \code{"two.sample"}
#' @param lower.tail Logical indicating whether Pr(BF <= k) (\code{TRUE}) or
#'     Pr(BF > k) (\code{FALSE}) should be computed. Defaults to \code{TRUE}
#' @param strict Logical indicating whether in case of one-sided alternatives
#'     the power should be computed also in the opposite direction. Defaults to
#'     \code{FALSE}
#'
#' @return The probability that the Bayes factor is less or greater (depending
#'     on the specified \code{lower.tail}) than the specified threshold \code{k}
#'
#' @author Samuel Pawel
#'
#' @seealso \link{tbf01}
#'
#' @examples
#' ptbf01(k = 1/6, n = 146, dpm = 0.5, dpsd = 0, alternative = "greater")
#' ptbf01(k = 6, n = 146, dpm = 0, dpsd = 0, alternative = "greater",
#'        lower.tail = FALSE)
#'
#' @export
ptbf01 <- Vectorize(FUN = ptbf01.,
                    vectorize.args = c("k", "n", "n1", "n2", "null",
                                       "plocation", "pscale", "pdf", "type",
                                       "alternative", "dpm", "dpsd", "type",
                                       "lower.tail", "strict"))

## ## verify with simulation
## set.seed(1000)
## nsim <- 10000
## dpm <- 0.5
## dpsd <- 0.1
## n1 <- 100
## n2 <- 120
## r <- 1/sqrt(2)
## bfs <- replicate(n = nsim, expr = {
##     if (dpsd == 0) m <- dpm
##     else m <- rnorm(n = 1, mean = dpm, sd = dpsd)
##     x <- rnorm(n = n1, mean = 0, sd = 1)
##     y <- rnorm(n = n2, mean = m, sd = 1)
##     t <- t.test(y, x)$statistic
##     ## se <- sqrt(1/n1 + 1/n2)
##     ## est <- rnorm(n = 1, mean = dpm, sd = sqrt(se^2 + dpsd^2))
##     ## t <- est/se
##     tbf01(t = t, n1 = n1, n2 = n2, pscale = r, type = "two.sample", alternative = "greater")
## })
## mean(bfs < 1/10)
## ptbf01(k = 1/10, n1 = n1, n2 = n2, pscale = r, dpm = dpm, dpsd = dpsd,
##        type = "two.sample", alternative = "greater")
