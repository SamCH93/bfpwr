jzsbf01. <- function(t, n, n1 = n, n2 = n, r = 1/sqrt(2),
                     type = c("two.sample", "one.sample",  "paired"),
                     log = FALSE, ...) {
    ## input checks
    stopifnot(
        length(t) == 1,
        is.numeric(t),
        is.finite(t),

        length(n1) == 1,
        is.numeric(n1),
        is.finite(n1),
        0 < n1,

        length(n2) == 1,
        is.numeric(n2),
        is.finite(n2),
        0 < n2,

        length(r) == 1,
        is.numeric(r),
        is.finite(r),
        0 < r,

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )
    type <- match.arg(type)
    if (type != "two.sample") {
        if (n1 != n2) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
        }
    }

    ## compute df and effective sample size depending on test type
    if (type == "two.sample") {
        df <- n1 + n2 - 2
        neff <- 1/(1/n1 + 1/n2)
    } else {
        df <- n1 - 1
        neff <- n1
    }

    ## compute JSZ Bayes factor
    f0 <- (1 + t^2/df)^(-0.5*(df + 1))
    intFun <- function(g) {
        (1 + neff*g*r^2)^(-0.5)*(1 + t^2/df/(1 + neff*g*r^2))^(-0.5*(df + 1))*
            (2*pi)^(-0.5)*g^(-1.5)*exp(-0.5/g)
    }
    f1 <- try(stats::integrate(f = intFun,
                               lower = .Machine$double.eps, # doesn't work with 0
                               upper = Inf, ... = ...)$value,
              silent = TRUE)

    if (inherits(f1, "try-error")) {
        bf <- NaN
    } else {
        bf <- f0/f1
    }

    if (log) return(log(bf))
    else return(bf)
}


#' @title Jeffreys-Zellner-Siow (JZS) Bayes factor
#'
#' @description This function computes the Jeffreys-Zellner-Siow Bayes factor
#'     that quantifies the evidence that the data provide for the null
#'     hypothesis that the standardized mean difference is zero vs. that the
#'     alternative that it is non-zero.
#'
#'  The data are summarized by \eqn{t}-statistics and sample sizes. The
#'     following types of \eqn{t}-tests are accepted:
#'
#' - Two-sample \eqn{t}-test where the SMD represents the standardized
#' mean difference between two group means (assuming equal variances in
#' both groups)
#' - One-sample \eqn{t}-test where the SMD represents the standardized
#' mean difference to the null value
#' - Paired \eqn{t}-test where the SMD represents the standardized mean
#' difference score
#'
#' The JZS Bayes factor is implemented as equation (1) in Rouder et al. (2014).
#' Integration is performed numerically with \code{stats::integrate}
#'
#' @param t \eqn{t}-statistic
#' @param n Sample size (per group)
#' @param n1 Sample size in group 1 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param n2 Sample size in group 2 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param r Scale parameter of the Cauchy prior. Defaults to \code{1/sqrt(2)}
#' @param type Type of \eqn{t}-test associated with \eqn{t}-statistic. Can be
#'     `"two.sample"`, `"one.sample"`, `"paired"`. Defaults to `"two.sample"
#' @param log Logical indicating whether natural logarithm of the Bayes factor
#'     should be returned. Defaults to \code{FALSE}
#' @param ... Additional arguments passed to \code{stats::integrate}
#'
#' @return Bayes factor in favor of the null hypothesis over the alternative (BF
#'     > 1 indicates evidence for the null hypothesis, whereas BF < 1 indicates
#'     evidence for the alternative)
#'
#' @author Samuel Pawel
#'
#' @references Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., Iverson,
#'     G. (2014). Bayesian t tests for accepting and rejecting the null
#'     hypothesis. Psychonomic Bulletin & Review, 16(2):225-237.
#'     \doi{10.3758/PBR.16.2.225}
#'
#' @examples
#' ## values from Table 1 in Rouder et al. (2019)
#' jzsbf01(t = c(0.69, 3.20), n = 100, r = 1, type = "one.sample")
#'
#' ## examples from p. 232 in Rouder et al. (2019)
#' jzsbf01(t = c(2.24, 2.03), n = 80, r = 1, type = "one.sample")
#'
#' @export
jzsbf01 <- Vectorize(FUN = jzsbf01.,
                     vectorize.args = c("t", "n", "n1", "n2", "r", "type",
                                        "log"))


## library(BayesFactor)
## n <- 100
## set.seed(44)
## y <- rnorm(n = n)
## t <- unname(t.test(y)$statistic)
## jzsbf01(t = t, n = n, r = 1/sqrt(2), type = "one.sample")
## 1/exp(ttest.tstat(t = t, n1 = n, r = 1/sqrt(2))$bf)

## set.seed(4456)
## n1 <- 100
## n2 <- 50
## x <- rnorm(n = n1)
## y <- rnorm(n = n2, mean = 0.5)
## t <- unname(t.test(x, y)$statistic)
## jzsbf01(t = t, n1 = n1, n2 = n2, r = 1/sqrt(2), type = "two.sample")
## 1/exp(ttest.tstat(t = t, n1 = n1, n2 = n2, r = 1/sqrt(2))$bf)

## set.seed(100)
## n <- 100
## y1 <- rnorm(n = n, mean = 0.5)
## y2 <- rnorm(n = n, mean = 0)
## t <- unname(t.test(y1, y2, paired = TRUE)$statistic)
## jzsbf01(t = t, n = n, r = 1/sqrt(2), type = "paired")
## 1/exp(ttest.tstat(t = t, n1 = n, r = 1/sqrt(2))$bf)
