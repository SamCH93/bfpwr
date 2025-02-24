## TODO think about interface and naming to select directional or point testing,
## and corresponding prior, an option could be to have an "alternative" argument
## as in binom.test, but that may be confusing because we are also changing the
## null
## TODO should we also add possibility test point vs. one-sided direction?
## should we add possibility to add point hypothesis as alternative?
binbf01. <- function(x, n, p0 = 0.5, type = c("point", "direction"), a = 1,
                     b = 1, log = FALSE) {
    ## input checks
    stopifnot(
        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(n) == 1,
        is.numeric(n),
        is.finite(n),
        0 < n,
        x <= n,

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

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )
    type <- match.arg(arg = type)

    if (type == "point") {
        logbf <- x*log(p0) + (n - x)*log(1 - p0) + lbeta(a = a, b = b) -
            lbeta(a = a + x, b = b + n - x)
    } else { # type == "direction"
        lpost0 <- stats::pbeta(q = p0, shape1 = a + x, shape2 = b + n - x,
                               log.p = TRUE)
        lpost1 <- stats::pbeta(q = p0, shape1 = a + x, shape2 = b + n - x,
                               lower.tail = FALSE, log.p = TRUE)
        lprior0 <- stats::pbeta(q = p0, shape1 = a, shape2 = b, log.p = TRUE)
        lprior1 <- stats::pbeta(q = p0, shape1 = a, shape2 = b,
                                lower.tail = FALSE, log.p = TRUE)
        logbf <- lpost0 - lpost1 + lprior1 - lprior0
    }
    if (log == TRUE) {
        return(logbf)
    } else {
        return(exp(logbf))
    }
}


#' @title Binomial Bayes factor
#'
#' @description This function computes the Bayes factor for testing a binomial
#'     proportion \eqn{p} based on \eqn{x} observed successes out of \eqn{n}
#'     trials. Two types of tests are available:
#'
#' * Test of a point null hypothesis: The Bayes factor quantifies the evidence
#' for \eqn{H_0 \colon p = p_0}{H0: p = p0} against \eqn{H_1 \colon p \neq
#' p_0}{H1: p != p0}. A beta prior is assigned to the proportion \eqn{p} under
#' the alternative hypothesis \eqn{H_1}{H1}.
#'
#' * Test of a directional null hypothesis: The Bayes factor quantifies the
#' evidence for \eqn{H_0 \colon p \leq p_0}{H0: p <= p0} against \eqn{H_1 \colon
#' p > p_0}{H1: p > p0}. A beta prior that is truncated to the range \eqn{[0,
#' p_0]}{[0, p0]} under the null \eqn{H_0}{H0} and to \eqn{(p_0, 1]}{(p0, 1]}
#' under the alternative \eqn{H_1}{H1} is assigned to the proportion \eqn{p}
#' under the corresponding hypothesis.
#'
#' @md
#'
#' @param x Number of successes
#' @param n Number of trials
#' @param p0 Tested binomial proportion. Defaults to \code{0.5}
#' @param type Type of test. Can be \code{"point"} or \code{"directional"}.
#'     Defaults to \code{"point"}
#' @param a Number of successes parameter of the beta prior distribution.
#'     Defaults to \code{1}
#' @param b Number of failures parameter of the beta prior distribution.
#'     Defaults to \code{1}
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#'
#' @inherit bf01 return
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbinbf01}, \link{nbinbf01}
#'
#' @examples
#' ## example on Mendelian inheritance from ?stats::binom.test
#' binbf01(x = 682, n = 925, p0 = 3/4, a = 1, b = 1, type = "point")
#' ## 18.6 => strong evidence for the hypothesized p = 3/4 compared to other p
#'
#' ## with directional hypothesis
#' binbf01(x = 682, n = 925, p0 = 3/4, a = 1, b = 1, type = "direction")
#' ## 1.5 => only anecdotal evidence for p <= 3/4 over p > 3/4
#'
#' ## Particle-counting experiment from Stone (1997) with point null
#' binbf01(x = 106298, n = 527135, p0 = 0.2, a = 1, b = 1, type = "point")
#' ## 8.1 => moderate evidence for the alternative over the null
#'
#' ## Coin flip experiment from Bartos et al. (2023) with point null
#' binbf01(x = 178079, n = 350757 , p0 = 0.5, a = 5100, b = 4900, type = "point")
#' ## => 1/1.72e+17 extreme evidence in favor of the alternative over the null
#'
#' @export
binbf01 <- Vectorize(FUN = binbf01.)
