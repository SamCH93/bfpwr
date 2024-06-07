ntbf01. <- function(k = 1/10, power, null = 0, plocation = 0,
                    pscale = 1/sqrt(2), pdf = 1,
                    type = c("two.sample", "one.sample", "paired"),
                    alternative = c("two.sided", "less", "greater"),
                    dpm = plocation, dpsd = pscale, lower.tail = TRUE,
                    integer = TRUE, nrange = c(2, 10^4), ...) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

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

         length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1],
        nrange[1] > 1,

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail),

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)


    ## define function for numerical root-finding
    rootFun <- function(n) {
        suppressWarnings({
            ptbf01(k = k, n = n, null = null, plocation = plocation,
                   pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd,
                   type = type, alternative = alternative,
                   lower.tail = lower.tail, ... = ...) - power
        })
    }

    ## check boundaries of sample size search range
    if (rootFun(n = nrange[1]) > 0) {
        warning("lower bound of sample size search range ('nrange') leads to higher power than specified")
        n <- NaN
    } else if (rootFun(n = nrange[2]) < 0) {
        warning("upper bound of sample size search range ('nrange') leads to lower power than specified")
        n <- NaN
    } else {
        ## perform root-finding
        res <- try(stats::uniroot(f = rootFun, interval = nrange, ... = ...)$root)
        if (inherits(res, "try-error")) {
            n <- NaN
        } else {
            n <- res
        }
    }

    if (integer) return(ceiling(n))
    else return(n)
}


#' @title Sample size calculations for \eqn{t}-test Bayes factor
#'
#' @description This function computes the required sample size to obtain a
#'     \eqn{t}-test Bayes factor (\link{tbf01}) less or greater than a threshold
#'     \code{k} with a specified target power.
#'
#' @note An error message will be displayed in case that the specified target
#'     power is not achievable under the specified analysis and design priors.
#'
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param power Target power
#' @param null Standardized mean difference under the point null hypothesis.
#'     Defaults to \code{0}
#' @param plocation Analysis \eqn{t} prior location. Defaults to \code{0}
#' @param pscale Analysis \eqn{t} prior scale. Defaults to \code{1/sqrt(2)}
#' @param pdf Analysis \eqn{t} prior degrees of freedom. Defaults to \code{1} (a
#'     Cauchy prior)
#' @param type Type of \eqn{t}-test associated with \eqn{t}-statistic. Can be
#'     \code{"two.sample"} (default), \code{"one.sample"}, or \code{"paired"}
#' @param alternative Direction of the test. Can be either \code{"two.sided"}
#'     (default), \code{"less"} , or \code{"greater"}. The latter two truncate
#'     the analysis prior to negative and positive effects, respectively
#' @param dpm Mean of the normal design prior assigned to the standardized mean
#'     difference. Defaults to the analysis prior location
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     standardized mean difference. Set to \code{0} to obtain a point prior at
#'     the design prior mean. Defaults to the analysis prior scale
#' @param type The type of test. One of \code{"two.sample"},
#'     \code{"one.sample"}, \code{"paired"}. Defaults to \code{"two.sample"}
#' @param lower.tail Logical indicating whether Pr(BF \eqn{\leq} \code{k})
#'     (\code{TRUE}) or Pr(BF \eqn{>} \code{k}) (\code{FALSE}) should be
#'     computed. Defaults to \code{TRUE}
#' @param integer Logical indicating whether only integer valued sample sizes
#'     should be returned. If \code{TRUE} the required sample size is rounded to
#'     the next larger integer. Defaults to \code{TRUE}
#' @param nrange Sample size search range over which numerical search is
#'     performed. Defaults to \code{c(2, 10^4)}
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return Object of class \code{"power.bftest"}, a list of the arguments
#'     (including the computed one) augmented with \code{method} and \code{note}
#'     elements
#'
#' @author Samuel Pawel
#'
#' @seealso \link{plot.power.bftest}, \link{ptbf01}, \link{powertbf01}
#'
#' @examples
#'  ## example from Schönbrodt and Wagenmakers (2018, p.135)
#'  ntbf01(k = 1/6, power = 0.95, dpm = 0.5, dpsd = 0, alternative = "greater")
#'  ntbf01(k = 1/6, power = 0.95, dpm = 0.5, dpsd = 0.1, alternative = "greater")
#'  ntbf01(k = 6, power = 0.95, dpm = 0.5, dpsd = 0, alternative = "greater",
#'         lower.tail = FALSE, nrange = c(2, 10000))
#'
#' @export
ntbf01 <- Vectorize(FUN = ntbf01.,
                    vectorize.args = c("k", "power", "null", "plocation",
                                       "pscale", "pdf", "type", "alternative",
                                       "dpm", "dpsd", "lower.tail", "integer"))
