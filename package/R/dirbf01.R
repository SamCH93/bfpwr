dirbf01. <- function(estimate, se, null = 0, pm, psd, log = FALSE) {
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
        0 < psd,

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )

    postsd <- 1/sqrt(1/se^2 + 1/psd^2)
    postm <- (estimate/se^2 + pm/psd^2)*postsd^2

    priorz <- (pm - null)/psd
    postz <- (postm - null)/postsd

    ## Compute directional prior/posterior odds on the log scale; otherwise
    ## very small tail masses collapse to zero or Inf.
    logpriorodds <- stats::pnorm(q = priorz, lower.tail = FALSE,
                                 log.p = TRUE) -
        stats::pnorm(q = priorz, lower.tail = TRUE, log.p = TRUE)
    logpostodds <- stats::pnorm(q = postz, lower.tail = FALSE,
                                log.p = TRUE) -
        stats::pnorm(q = postz, lower.tail = TRUE, log.p = TRUE)

    logbf <- logpostodds - logpriorodds

    if (log) return(logbf)
    else return(exp(logbf))
}


#' @title Directional z-test Bayes factor
#'
#' @description This function computes the Bayes factor that quantifies the
#'     evidence that the data (in the form of an asymptotically normally
#'     distributed parameter estimate with standard error) provide for a
#'     directional null hypothesis that the the parameter value is less than the
#'     null value against the alternative that it is greater than the null
#'     value. A marginal normal prior is assigned to the parameter. The standard
#'     error is assumed to be known.
#'
#' @param estimate Parameter estimate
#' @param se Standard error of the parameter estimate
#' @param null Null value that separates the null from the alternative
#'     hypothesis. Defaults to \code{0}
#' @param pm Mean of the normal prior assigned to the parameter
#' @param psd Standard deviation of the normal prior assigned to the parameter
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#'
#' @return Bayes factor in favor of the null hypothesis over the alternative
#'     (\eqn{\text{BF}_{01}}{BF01} > 1 indicates evidence for the null
#'     hypothesis, whereas \eqn{\text{BF}_{01}}{BF01} < 1 indicates evidence for
#'     the alternative)
#'
#' @author Samuel Pawel
#'
#' @examples
#' dirbf01(estimate = 0.2, se = 0.2, null = 0, pm = 0, psd = 2)
#'
#' @export
dirbf01 <- Vectorize(FUN = dirbf01.)
