bf01. <- function(estimate, se, null = 0, pm, psd, log = FALSE) {
    ## input checks
    stopifnot(
        length(estimate) == 1,
        is.numeric(estimate),
        is.finite(estimate),

        length(se) == 1,
        is.numeric(se),
        is.finite(se),
        0 < se,

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

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )

    logbf <- stats::dnorm(x = estimate, mean = null, sd = se, log = TRUE) -
        stats::dnorm(x = estimate, mean = pm, sd = sqrt(se^2 + psd^2), log = TRUE)
    if (log) return(logbf)
    else return(exp(logbf))
}


#' @title z-test Bayes factor
#'
#' @description This function computes the Bayes factor that quantifies the
#'     evidence that the data (in the form of an asymptotically normally
#'     distributed parameter estimate with standard error) provide for a point
#'     null hypothesis with a normal prior assigned to the parameter under the
#'     alternative. The standard error is assumed to be known.
#'
#' @param estimate Parameter estimate
#' @param se Standard error of the parameter estimate
#' @param null Parameter value under the point null hypothesis. Defaults to
#'     \code{0}
#' @param pm Mean of the normal prior assigned to the parameter under the
#'     alternative
#' @param psd Standard deviation of the normal prior assigned to the parameter
#'     under the alternative. Set to \code{0} to obtain a point prior at the
#'     prior mean
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#'
#' @return Bayes factor in favor of the null hypothesis over the alternative (BF
#'     > 1 indicates evidence for the null hypothesis, whereas BF < 1 indicates
#'     evidence for the alternative)
#'
#' @author Samuel Pawel
#'
#' @examples
#' bf01(estimate = 0.2, se = 0.05, null = 0, pm = 0, psd = 2)
#'
#' @export
bf01 <- Vectorize(FUN = bf01.)
