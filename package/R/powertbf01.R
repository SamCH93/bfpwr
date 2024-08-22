#' @title Power and sample size calculations for \eqn{t}-test Bayes factor
#'
#' @description Compute probability that \eqn{t}-test Bayes factor is smaller
#'     than a specified threshold (the power), or determine sample size to
#'     obtain a target power
#'
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param n Sample size (per group for two-sample tests). Has to be \code{NULL}
#'     if \code{power} is specified. Defaults to \code{NULL}
#' @param power Target power. Has to be \code{NULL} if \code{n} is specified.
#'     Defaults to \code{NULL}
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
#'     the analysis prior to negative and positive effects, respectively. If set
#'     to \code{"less"} or \code{"greater"}, the power is only computed based on
#'     data with effect estimates in the direction of the alternative
#' @param dpm Mean of the normal design prior assigned to the standardized mean
#'     difference. Defaults to the analysis prior location
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     standardized mean difference. Set to \code{0} to obtain a point prior at
#'     the design prior mean. Defaults to the analysis prior scale
#' @param type The type of test. One of \code{"two.sample"},
#'     \code{"one.sample"}, \code{"paired"}. Defaults to \code{"two.sample"}
#' @param nrange Sample size search range over which numerical search is
#'     performed (only taken into account when \code{n} is \code{NULL}).
#'     Defaults to \code{c(2, 10^4)}
#'
#' @return Object of class \code{"power.bftest"}, a list of the arguments
#'     (including the computed one) augmented with \code{method} and \code{note}
#'     elements
#'
#' @author Samuel Pawel
#'
#' @seealso \link{plot.power.bftest}, \link{ptbf01}, \link{ntbf01}, \link{tbf01}
#'
#' @examples
#' ## determine power
#' powertbf01(k = 1/6, n = 146, dpm = 0.5, dps = 0, alternative = "greater")
#'
#' ## determine sample size
#' powertbf01(k = 1/6, power = 0.95, dpm = 0.5, dps = 0, alternative = "greater")
#'
#' @export
powertbf01 <- function(k = 1/10, n = NULL, power = NULL, null = 0,
                       plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                       type = c("two.sample", "one.sample", "paired"),
                       alternative = c("two.sided", "less", "greater"),
                       dpm = plocation, dpsd = pscale, nrange = c(2, 10^4)) {
    ## input checks
    if (is.null(n) && is.null(power)) {
        stop("exactly one of 'n' and 'power' must be NULL")
    }
    if (is.null(n)) {
        stopifnot(
            length(power) == 1,
            is.numeric(power),
            is.finite(power),
            0 < power, power < 1
        )
    } else {
        stopifnot(
            length(n) == 1,
            is.numeric(n),
            is.finite(n),
            0 < n
        )
    }
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

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
        0 <= dpsd
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)

    ## determine sample size
    if (is.null(n)) {
        n <- ntbf01(k = k, power = power, null = null, plocation = plocation,
                    pscale = pscale, pdf = pdf, type = type,
                    alternative = alternative, dpm = dpm, dpsd = dpsd,
                    integer = FALSE)
    } else {
        ## determine power
        power <- ptbf01(k = k, n = n, null = null, plocation = plocation,
                        pscale = pscale, pdf = pdf, type = type,
                        alternative = alternative, dpm = dpm, dpsd = dpsd)
    }

    ## return object
    structure(list(n = n, power = power, sd = 1, null = null,
                   alternative = alternative, plocation = plocation,
                   pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd, k = k,
                   nrange = nrange, type = type, test = "t"),
              class = "power.bftest")

}
