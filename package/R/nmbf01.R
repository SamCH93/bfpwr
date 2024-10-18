nmbf01. <- function(estimate, se, null = 0, psd, log = FALSE) {
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

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        0 < psd,

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )

    logbf <- 1.5*log(1 + psd^2/se^2) - 0.5*(estimate - null)^2/se^2/(1 + se^2/psd^2) -
        log(1 + (estimate - null)^2/se^2/(1 + se^2/psd^2))

    if (log) return(logbf)
    else return(exp(logbf))
}


#' @title Normal moment prior Bayes factor
#'
#' @description This function computes the Bayes factor that quantifies the
#'     evidence that the data (in the form of an asymptotically normally
#'     distributed parameter estimate with standard error) provide for a point
#'     null hypothesis with a normal moment prior assigned to the parameter
#'     under the alternative.
#'
#' @details A normal moment prior has density \eqn{f(x \mid \code{null},
#'     \code{psd}) = N(x \mid \code{null}, \code{psd}^2) \times (x -
#'     \code{null})/ \code{psd}^2}{f(x|\code{null},\code{psd}) =
#'     N(x|\code{null},\code{psd})*(x - \code{null})^2/\code{psd}^2} with
#'     \eqn{N(x \mid m, v)}{N(x|m,v)} the normal density with mean \eqn{m} and
#'     variance \eqn{v} evaluated at \eqn{x}.
#'
#' @inheritParams bf01
#' @param psd Spread of the normal moment prior assigned to the parameter under
#'     the alternative. The modes of the prior are located at
#'     \eqn{\pm\sqrt{2}\,\code{psd}}{+-sqrt(2)*\code{psd}}
#'
#' @inherit bf01 return
#'
#' @author Samuel Pawel
#'
#' @references Johnson, V. E. and Rossell, D. (2010). On the use of non-local
#'     prior densities in Bayesian hypothesis tests. Journal of the Royal
#'     Statistical Society: Series B (Statistical Methodology), 72(2):143–170.
#'     \doi{10.1111/j.1467-9868.2009.00730.x}
#'
#' Pramanik, S. and Johnson, V. E. (2024). Efficient alternatives for Bayesian
#'     hypothesis tests in psychology. Psychological Methods, 29(2):243–261.
#'     \doi{10.1037/met0000482}
#'
#' @examples
#' nmbf01(estimate = 0.25, se = 0.05, null = 0, psd = 0.5/sqrt(2)) # mode at 0.5
#'
#' @seealso \link{nmbf01}, \link{pnmbf01}, \link{nnmbf01}, \link{powernmbf01}
#'
#' @export
nmbf01 <- Vectorize(FUN = nmbf01.)
