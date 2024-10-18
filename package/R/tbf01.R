tbf01. <- function(t, n, n1 = n, n2 = n, plocation = 0, pscale = 1/sqrt(2),
                   pdf = 1, type = c("two.sample", "one.sample",  "paired"),
                   alternative = c("two.sided", "less", "greater"), log = FALSE,
                   ...) {
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

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
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

    ## marginal likelihood under the null hypothesis
    f0 <- stats::dt(x = t, df = df, ncp = 0)

    ## marginal likelihood under the alternative hypothesis
    if (alternative == "two.sided") {
        normConst <- 1
        lower <- -Inf
        upper <- Inf
    } else if (alternative == "greater") {
        normConst <- 1 - stats::pt(q = (0 - plocation)/pscale, df = pdf)
        lower <- 0
        upper <- Inf
    } else {
        normConst <- stats::pt(q = (0 - plocation)/pscale, df = pdf)
        lower <- -Inf
        upper <- 0
    }
    dpriorH1 <- function(d) {
        stats::dt(x = (d - plocation)/pscale, df = pdf, ncp = 0)/(pscale*normConst)
    }
    intFun <- function(d) {
        suppressWarnings({
            stats::dt(x = t, df = df, ncp = sqrt(neff)*d)*dpriorH1(d)
        })
    }
    f1 <- try(stats::integrate(f = intFun, lower = lower, upper = upper,
                               ... = ...)$value, silent = TRUE)

    ## compute Bayes factor
    if (inherits(f1, "try-error")) {
        bf <- NaN
    } else {
        bf <- f0/f1
    }
    if (log) return(log(bf))
    else return(bf)
}


#' @title t-test Bayes factor
#'
#' @description This function computes the Bayes factor that forms the basis of
#'     the informed Bayesian \eqn{t}-test from Gronau et al. (2020). The Bayes
#'     factor quantifies the evidence that the data provide for the null
#'     hypothesis that the standardized mean difference (SMD) is zero against
#'     the alternative that the SMD is non-zero. A location-scale
#'     \eqn{t}-distribution is assumed for the SMD under the alternative
#'     hypothesis. The Jeffreys-Zellner-Siow (JZS) Bayes factor (Rouder et al.,
#'     2009) is obtained as a special case by setting the location of the prior
#'     to zero and the prior degrees of freedom to one, which is the default.
#'
#'  The data are summarized by \eqn{t}-statistics and sample sizes. The
#'     following types of \eqn{t}-statistics are accepted:
#'
#' - Two-sample \eqn{t}-test where the SMD represents the standardized
#'   mean difference between two group means (assuming equal variances in
#'   both groups)
#' - One-sample \eqn{t}-test where the SMD represents the standardized
#'    mean difference to the null value
#' - Paired \eqn{t}-test where the SMD represents the standardized mean
#'   change score
#'
#' @md
#'
#' @details The Bayes factor is implemented as in equation (5) in Gronau et al.
#'     (2020), and using suitable truncation in case of one-sided alternatives.
#'     Integration is performed numerically with \code{stats::integrate}.
#'
#' @param t \eqn{t}-statistic
#' @param n Sample size (per group)
#' @param n1 Sample size in group 1 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param n2 Sample size in group 2 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param plocation \eqn{t} prior location. Defaults to \code{0}
#' @param pscale \eqn{t} prior scale. Defaults to \code{1/sqrt(2)}
#' @param pdf \eqn{t} prior degrees of freedom. Defaults to \code{1} (a Cauchy
#'     prior)
#' @param type Type of \eqn{t}-test. Can be \code{"two.sample"} (default),
#'     \code{"one.sample"}, or \code{"paired"}
#' @param alternative Direction of the test. Can be either \code{"two.sided"}
#'     (default), \code{"less"}, or \code{"greater"}. The latter two truncate
#'     the analysis prior to negative and positive effects, respectively.
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#' @param ... Additional arguments passed to \code{stats::integrate}
#'
#' @inherit bf01 return
#'
#' @author Samuel Pawel
#'
#' @references Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., Iverson,
#'     G. (2009). Bayesian \eqn{t} tests for accepting and rejecting the null
#'     hypothesis. Psychonomic Bulletin & Review, 16(2):225-237.
#'     \doi{10.3758/PBR.16.2.225}
#'
#' Gronau, Q. F., Ly., A., Wagenmakers, E.J. (2020). Informed Bayesian
#'     \eqn{t}-Tests. The American Statistician, 74(2):137-143.
#'     \doi{10.1080/00031305.2018.1562983}
#'
#' @seealso \link{powertbf01}, \link{ptbf01}, \link{ntbf01}
#'
#' @examples
#' ## analyses from Rouder et al. (2009):
#' ## values from Table 1
#' tbf01(t = c(0.69, 3.20), n = 100, pscale = 1, type = "one.sample")
#' ## examples from p. 232
#' tbf01(t = c(2.24, 2.03), n = 80, pscale = 1, type = "one.sample")
#'
#' ## analyses from Gronau et al. (2020) section 3.2:
#' ## informed prior
#' tbf01(t = -0.90, n1 = 53, n2 = 57, plocation = 0.350, pscale = 0.102, pdf = 3,
#'       alternative = "greater", type = "two.sample")
#' ## default (one-sided) prior
#' tbf01(t = -0.90, n1 = 53, n2 = 57, plocation = 0, pscale = 1/sqrt(2), pdf = 1,
#'       alternative = "greater", type = "two.sample")
#'
#' @export
tbf01 <- Vectorize(FUN = tbf01.,
                   vectorize.args = c("t", "n", "n1", "n2", "plocation",
                                      "pscale", "pdf", "type", "alternative",
                                      "log"))

## ## verify results with BayesFactor package
## library(BayesFactor)
## n <- 100
## set.seed(44)
## y <- rnorm(n = n)
## t <- unname(t.test(y)$statistic)
## tbf01(t = t, n = n, r = 1/sqrt(2), type = "one.sample")
## 1/exp(ttest.tstat(t = t, n1 = n, r = 1/sqrt(2))$bf)

## set.seed(4456)
## n1 <- 100
## n2 <- 50
## x <- rnorm(n = n1)
## y <- rnorm(n = n2, mean = 0.5)
## t <- unname(t.test(x, y)$statistic)
## tbf01(t = t, n1 = n1, n2 = n2, r = 1/sqrt(2), type = "two.sample")
## 1/exp(ttest.tstat(t = t, n1 = n1, n2 = n2, r = 1/sqrt(2))$bf)

## set.seed(100)
## n <- 100
## y1 <- rnorm(n = n, mean = 0.5)
## y2 <- rnorm(n = n, mean = 0)
## t <- unname(t.test(y1, y2, paired = TRUE)$statistic)
## tbf01(t = t, n = n, r = 1/sqrt(2), type = "paired")
## 1/exp(ttest.tstat(t = t, n1 = n, r = 1/sqrt(2))$bf)
