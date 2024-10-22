#' @title Power and sample size calculations for \eqn{t}-test Bayes factor
#'
#' @description Compute probability that \eqn{t}-test Bayes factor is smaller
#'     than a specified threshold (the power), or determine sample size to
#'     obtain a target power.
#'
#' @details This function provides a similar interface as
#'     \code{stats::power.t.test}. For some users, the low-level functions
#'     \link{ntbf01} (to directly compute the sample size for a fixed power) and
#'     \link{ptbf01} (to directly compute the power for a fixed sample size) may
#'     also be useful.
#'
#' @inheritParams ptbf01
#' @inheritParams powerbf01
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param nrange Sample size search range over which numerical search is
#'     performed (only taken into account when \code{n} is \code{NULL}).
#'     Defaults to \code{c(2, 10^4)}
#'
#' @inherit powerbf01 return
#'
#' @author Samuel Pawel
#'
#' @seealso \link{plot.power.bftest}, \link{ptbf01}, \link{ntbf01}, \link{tbf01}
#'
#' @examples
#' ## determine power
#' powertbf01(n = 146, k = 1/6, dpm = 0.5, dps = 0, alternative = "greater")
#'
#' ## determine sample size
#' powertbf01(power = 0.95, k = 1/6, dpm = 0.5, dps = 0, alternative = "greater")
#'
#' @export
powertbf01 <- function(n = NULL, power = NULL, k = 1/10, null = 0,
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
