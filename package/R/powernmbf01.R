#' @title Power and sample size calculations for normal moment prior Bayes
#'     factor
#'
#' @description Compute probability that normal moment prior Bayes factor is
#'     smaller than a specified threshold (the power), or determine sample size
#'     to obtain a target power
#'
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param n Sample size. Has to be \code{NULL} if \code{power} is specified.
#'     Defaults to \code{NULL}
#' @param power Target power. Has to be \code{NULL} if \code{n} is specified.
#'     Defaults to \code{NULL}
#' @param sd Standard deviation of one observation (for \code{type =
#'     "two.sample"} or \code{type = "one.sample"}) or of one difference within
#'     a pair of observations (\code{type = "paired"}). Is assumed to be known.
#'     Defaults to \code{1}
#' @param null Mean difference under the point null hypothesis. Defaults to
#'     \code{0}
#' @param psd Spread of the normal moment prior assigned to the parameter under
#'     the alternative in the analysis. The modes of the prior are located at
#'     \eqn{\pm\sqrt{2}\,\code{psd}}{+-sqrt(2)*\code{psd}}
#' @param type The type of test. One of \code{"two.sample"},
#'     \code{"one.sample"}, \code{"paired"}. Defaults to \code{"two.sample"}
#' @param dpm Mean of the normal design prior assigned to the parameter
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     parameter. Set to 0 to obtain a point prior at the prior mean
#' @param nrange Sample size search range over which numerical search is
#'     performed (only taken into account when \code{n} is \code{NULL}).
#'     Defaults to \code{c(1, 10^5)}
#'
#' @return Object of class \code{"power.bftest"}, a list of the arguments
#'     (including the computed one) augmented with \code{method} and \code{note}
#'     elements
#'
#' @author Samuel Pawel
#'
#' @seealso \link{plot.power.bftest}, \link{nnmbf01}, \link{pnmbf01}, \link{nmbf01}
#'
#' @examples
#' ## determine power
#' powernmbf01(n = 100, psd = 1, dpm = 0.5, dpsd = 0)
#'
#' ## determine sample size
#' powernmbf01(power = 0.99, psd = 1, dpm = 0.5, dpsd = 0)
#'
#' @export
powernmbf01 <- function(n = NULL, power = NULL, k = 1/10, sd = 1, null = 0, psd,
                        type = c("two.sample", "one.sample", "paired"), dpm,
                        dpsd, nrange = c(1, 10^5)) {
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

        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        0 < sd,

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

        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1]
    )
    type <- match.arg(type)

    ## determine unit variance
    if (type == "two.sample") {
        uv <- 2*sd^2
    } else {
        uv <- sd^2
    }

    ## determine sample size
    if (is.null(n)) {
        n <- nnmbf01(k = k, power = power, sd = sqrt(uv), null = null,
                     psd = psd, dpm = dpm, dpsd = dpsd, nrange = nrange,
                     integer = FALSE)
    } else {
        ## determine power
        power <- pnmbf01(k = k, n = n, sd = sqrt(uv), null = null, psd = psd,
                         dpm = dpm, dpsd = dpsd, lower.tail = TRUE)
    }

    ## return object
    structure(list(n = n, power = power, sd = sd, null = null, psd = psd,
                   dpm = dpm, dpsd = dpsd, k = k, nrange = nrange, type = type,
                   test = "nm"),
              class = "power.bftest")

}
