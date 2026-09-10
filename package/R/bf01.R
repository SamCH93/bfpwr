bf01. <- function(estimate, se, null = 0, pm, psd, log = FALSE,
                  alternative = c("two.sided", "less", "greater")) {
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

    alternative <- match.arg(alternative)
    .bf01_check_alternative(alternative, pm = pm, psd = psd, null = null)
    if (alternative != "two.sided" && psd > 0) {
        direction <- if (alternative == "greater") 1 else -1
        logbf <- .bf01_log_one_sided(z = direction*(estimate - null)/se,
                                    m = direction*(pm - null)/se, r = psd/se)
    } else {
        logbf <- stats::dnorm(x = estimate, mean = null, sd = se, log = TRUE) -
            stats::dnorm(x = estimate, mean = pm, sd = sqrt(se^2 + psd^2), log = TRUE)
    }
    if (log) return(logbf)
    else return(exp(logbf))
}


## A one-sided point prior must have positive mass on the retained side.
.bf01_check_alternative <- function(alternative, pm, psd, null,
                                    type = "normal") {
    if (alternative == "two.sided") return(invisible(NULL))
    if (type != "normal") {
        stop("one-sided 'alternative' requires a normal analysis prior (type = \"normal\")")
    }
    if (psd == 0 && ((alternative == "greater" && pm <= null) ||
                    (alternative == "less" && pm >= null))) {
        stop("a one-sided point prior must lie strictly on the specified side of 'null'")
    }
    invisible(NULL)
}


## log(Phi(x)) + x^2/2, using the normal-tail expansion to avoid subtracting
## nearly equal large numbers when x is far into the negative tail.
.bf01_log_scaled_pnorm <- function(x) {
    if (x >= -20) return(stats::pnorm(x, log.p = TRUE) + x^2/2)
    q <- 1/x^2
    -log(2*pi)/2 - log(-x) +
        log1p(q*(-1 + q*(3 + q*(-15 + q*(105 + q*(-945 + q*10395))))))
}


## Point-null vs. positive truncated normal BF (Table 1, Pawel and Held,
## Bayes Factor Group Sequential Designs). z, m and r are in standard-error
## units. Reflection gives the negative alternative. The BF is the ordinary
## normal BF times prior/posterior probability of the retained half-line.
.bf01_log_one_sided <- function(z, m, r) {
    r2 <- r^2
    priorZ <- m/r
    posteriorZ <- (priorZ + r*z)/sqrt(1 + r2)
    if (priorZ < -20 && posteriorZ < -20) {
        ## Both normal tails are tiny. Combine their leading terms before
        ## taking logs, and factor the difference of the tail polynomials.
        ## This also preserves BF01 - 1 for priors concentrated near the null.
        shift <- r*z/priorZ
        priorQ <- 1/priorZ^2
        posteriorQ <- 1/posteriorZ^2
        deltaQ <- priorQ*(2*shift + shift^2 - r2)/(1 + shift)^2
        coefficients <- c(-1, 3, -15, 105, -945, 10395)
        posteriorSeries <- difference <- 0
        powerSum <- 1
        for (i in seq_along(coefficients)) {
            posteriorSeries <- posteriorSeries + coefficients[i]*posteriorQ^i
            ## (a^i - b^i)/(a - b) = a^(i-1) + ... + b^(i-1).
            difference <- difference + coefficients[i]*powerSum
            powerSum <- posteriorQ*powerSum + priorQ^i
        }
        return(log1p(shift) + log1p(deltaQ*difference/(1 + posteriorSeries)))
    }
    if (priorZ >= 0 && posteriorZ >= 0) {
        ## Expanded normal BF avoids cancellation for narrow shifted priors.
        logbf <- log1p(r2)/2 + (m*(m - 2*z) - r2*z^2)/(2*(1 + r2))
        return(logbf + stats::pnorm(priorZ, log.p = TRUE) -
                   stats::pnorm(posteriorZ, log.p = TRUE))
    }
    log1p(r2)/2 + .bf01_log_scaled_pnorm(priorZ) -
        .bf01_log_scaled_pnorm(posteriorZ)
}


#' @title Point null z-test Bayes factor
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
#' @param alternative Direction of the alternative hypothesis, one of
#'     \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#'     For a one-sided alternative, the normal analysis prior is restricted to
#'     values below \code{null} (\code{"less"}) or above \code{null}
#'     (\code{"greater"}) and renormalized. \code{pm} and \code{psd} describe
#'     the normal distribution before truncation. A point prior
#'     (\code{psd = 0}) must lie strictly on the specified side of \code{null}.
#'     In power and sample-size calculations, the normal design prior is
#'     not truncated.
#'
#' @return Bayes factor in favor of the null hypothesis over the alternative
#'     (\eqn{\text{BF}_{01}}{BF01} > 1 indicates evidence for the null
#'     hypothesis, whereas \eqn{\text{BF}_{01}}{BF01} < 1 indicates evidence for
#'     the alternative)
#'
#' @author Samuel Pawel
#'
#' @examples
#' bf01(estimate = 0.2, se = 0.05, null = 0, pm = 0, psd = 2)
#' bf01(estimate = 0.2, se = 0.05, pm = 0, psd = 2, alternative = "greater")
#'
#' @export
bf01 <- Vectorize(FUN = bf01.)
