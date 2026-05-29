#' @title Sequential Correlations Test Bayes Factor Design
#'
#' @description Computes cumulative probabilities of observing Fisher-\eqn{z}
#'     correlation test Bayes factors that provide evidence for the null
#'     hypothesis \eqn{H_0}{H0}, the alternative hypothesis \eqn{H_1}{H1}, or
#'     remain inconclusive in a sequential design. Optionally, also computes the
#'     expected sample size.
#'
#' @inheritParams pbf01seq
#' @param n Numeric vector of sample sizes. The minimal sample size has to be at
#'     least 4 for estimating the standard error of a z-transformed correlation
#' @param pm Mean of the analysis prior assigned to the \eqn{z}-transformed
#'     correlation. Not taken into account for \code{type = "moment"}
#' @param psd Standard deviation (\code{type = "moment"} and \code{type =
#'     "directional"}) or scale (\code{type = "moment"}) of the analysis prior
#'     assigned the \eqn{z}-transformed correlation
#' @param dpm Mean of the normal design prior assigned the \eqn{z}-transformed
#'     correlation
#' @param dpsd Standard deviation of the normal design prior assigned the
#'     \eqn{z}-transformed correlation. Set \code{dpsd = 0} to obtain a point
#'     prior at \code{dpm}
#'
#' @return An object of class \code{"bfseqdesign"}, which is a list containing
#'     the input arguments, the critical z-values, the expected sample size, the
#'     cumulative probabilities of stopping for \eqn{H_1}{H1} and \eqn{H_0}{H0}
#'     by each stage, and the cumulative probabilities of remaining inconclusive
#'     by each stage.
#'
#' @details The function constructs per-stage integration regions for cumulative
#'     z-statistics based on the Bayes factor thresholds \code{k1} and
#'     \code{k0}, then computes the probability of these regions under a
#'     predictive distribution defined by \code{se} and the normal design prior
#'     with \code{dpm} and \code{dpsd}. Integration is performed via
#'     \code{mvtnorm::lpmvnorm}.
#'
#' @examples
#' n <- seq(30, 90, 30) # sample size per stage
#' res <- pcorbf01seq(k1 = 1/10, k0 = 3, n = n, pm = 0, psd = 1,
#'                    dpm = 0.3, dpsd = 0.05, type = "moment")
#' res # print summary
#' plot(res) # plot summary
#' res$cumpH1 # cumulative probability to stop for H1 by each stage
#' res$cumpH0 # cumulative probability to stop for H0 by each stage
#' res$EN # expected sample size
#'
#' @author Samuel Pawel
#'
#' @export
pcorbf01seq <- function(k1, k0 = 1/k1, n, pm, psd, dpm = pm, dpsd = psd,
                        type = c("normal", "directional", "moment"),
                        strict = TRUE, ...) {

    ## require n > 3 for standard error to be defined
    stopifnot(
        is.numeric(n),
        all(is.finite(n)),
        all(n > 3)
    )

    ## compute standard error
    se <- 1/sqrt(n - 3)

    ## call the generic z-test function with the computed standard error
    res <- pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm, psd, dpm = dpm,
                    dpsd = dpsd, type = type, strict = strict, ...)
    res$test <- "zcor"

    return(res)
}
