#' @title Power and sample size calculations for normal moment prior Bayes
#'     factor
#'
#' @description Compute probability that normal moment prior Bayes factor is
#'     smaller than a specified threshold (the power), or determine sample size
#'     to obtain a target power.
#'
#' @details This function provides a similar interface as
#'     \code{stats::power.t.test}. It also assumes that the data are continuous
#'     and that the parameter of interest is either a mean or a (standardized)
#'     mean difference. For some users, the low-level functions \link{nnmbf01}
#'     (to directly compute the sample size for a fixed power) and
#'     \link{pnmbf01} (to directly compute the power for a fixed sample size)
#'     may also be useful because they can be used for other data and parameter
#'     types.
#'
#' @inheritParams pnmbf01
#' @inheritParams powerbf01
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param n Sample size (per group for two-sample tests). Has to be \code{NULL}
#'     if \code{power} is specified. Defaults to \code{NULL}
#'
#' @inherit powerbf01 return
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
        n <- nnmbf01(k = k, power = power, usd = sqrt(uv), null = null,
                     psd = psd, dpm = dpm, dpsd = dpsd, nrange = nrange,
                     integer = FALSE)
    } else {
        ## determine power
        power <- pnmbf01(k = k, n = n, usd = sqrt(uv), null = null, psd = psd,
                         dpm = dpm, dpsd = dpsd, lower.tail = TRUE)
    }

    ## return object
    structure(list(n = n, power = power, sd = sd, null = null, psd = psd,
                   dpm = dpm, dpsd = dpsd, k = k, nrange = nrange, type = type,
                   test = "nm"),
              class = "power.bftest")

}
