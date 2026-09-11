## Route clearly wrong-tail one-sided statistics through stable quadrature.
.tbf01_tail_quadrature_cutoff <- 4

.tbf01_valid_tail_nquad <- function(tail.nquad) {
    length(tail.nquad) == 1 && is.numeric(tail.nquad) &&
        is.finite(tail.nquad) && tail.nquad >= 2 &&
        tail.nquad == floor(tail.nquad)
}

.tbf01_valid_drange <- function(drange) {
    (is.numeric(drange) && length(drange) == 2 && all(is.finite(drange)) &&
     drange[2] > drange[1]) || (is.character(drange) &&
                                length(drange) == 1 && !is.na(drange) &&
                                drange == "adaptive")
}

.tbf01_pars <- function(n1, n2, type) {
    ## Effective sample size for the noncentrality parameter sqrt(neff)*d.
    if (type == "two.sample") {
        list(df = n1 + n2 - 2, neff = 1/(1/n1 + 1/n2))
    } else {
        list(df = n1 - 1, neff = n1)
    }
}

.tbf01_prior_region <- function(plocation, pscale, pdf, alternative) {
    ## One-sided alternatives truncate and renormalize the analysis prior.
    q0 <- (0 - plocation)/pscale
    if (alternative == "two.sided") {
        list(lower = -Inf, upper = Inf, log_norm_const = 0)
    } else if (alternative == "greater") {
        list(lower = 0, upper = Inf,
             log_norm_const = stats::pt(q = q0, df = pdf, lower.tail = FALSE,
                                        log.p = TRUE))
    } else {
        list(lower = -Inf, upper = 0,
             log_norm_const = stats::pt(q = q0, df = pdf, lower.tail = TRUE,
                                        log.p = TRUE))
    }
}

.tbf01_log_fast <- function(t, df, neff, plocation, pscale, pdf, region,
                            rel.tol = .bfpwr_defaults$rel.tol, abs.tol = rel.tol,
                            subdivisions = .bfpwr_defaults$subdivisions, ...) {
    ## Original one-dimensional integral, evaluated on a centered log scale.
    if (!is.finite(region$log_norm_const)) {
        return(NaN)
    }

    eta <- sqrt(neff)
    log_f0 <- stats::dt(x = t, df = df, log = TRUE)

    log_prior <- function(d) {
        stats::dt(x = (d - plocation)/pscale, df = pdf, log = TRUE) -
            log(pscale) - region$log_norm_const
    }
    log_integrand <- function(d) {
        suppressWarnings({
            stats::dt(x = t, df = df, ncp = eta*d, log = TRUE) + log_prior(d)
        })
    }

    centers <- c(t/eta, plocation, 0)
    centers <- centers[is.finite(centers) & centers >= region$lower &
                       centers <= region$upper]
    if (length(centers) == 0) {
        centers <- if (is.finite(region$lower)) region$lower else region$upper
    }
    log_centers <- vapply(centers, log_integrand, numeric(1))
    log_center <- max(log_centers, na.rm = TRUE)
    if (!is.finite(log_center)) {
        return(NaN)
    }

    intfun <- function(d) {
        z <- log_integrand(d) - log_center
        z[!is.finite(z)] <- -Inf
        exp(z)
    }
    f1 <- try(stats::integrate(f = intfun, lower = region$lower,
                               upper = region$upper, rel.tol = rel.tol,
                               abs.tol = abs.tol, subdivisions = subdivisions,
                               ...)$value,
              silent = TRUE)
    if (inherits(f1, "try-error") || !is.finite(f1) || f1 <= 0) {
        return(NaN)
    }

    log_f0 - (log(f1) + log_center)
}

.tbf01_log_tail_quadrature <- function(t, df, neff, plocation, pscale, pdf,
                                       region, tail.nquad) {
    ## Wrong-tail one-sided calls can underflow in the noncentral-t density used
    ## by the direct integral. Conditional on the observed statistic under H0,
    ## its chi-square mixing variable has the gamma distribution below. The
    ## alternative-to-null density ratio is then averaged over this posterior
    ## and directly over the (possibly truncated) t prior. Integrating the
    ## truncated prior on probability scale avoids dividing a poorly resolved
    ## tail integral by a very small prior mass.
    if (!is.finite(region$log_norm_const)) {
        return(NaN)
    }

    nquad <- as.integer(tail.nquad)
    eta <- sqrt(neff)
    shapeV <- (df + 1)/2
    rateV <- (1 + t^2/df)/2

    quadrature <- .bfpwr_gauss_legendre(nquad)
    u <- quadrature$x
    w <- quadrature$w
    v <- stats::qgamma(p = u, shape = shapeV, rate = rateV)
    ## In the wrong tail the likelihood is concentrated close to the truncation
    ## point d = 0. A power-transformed probability scale places
    ## substantially more nodes there while leaving the integral exact after
    ## its Jacobian is included. This matters for narrow shifted priors, where
    ## the relevant prior probability can be far below the smallest ordinary
    ## Gauss-Legendre node.
    boundaryPower <- 16
    priorU <- u
    priorLogJacobian <- rep(0, nquad)
    if (is.finite(region$lower)) {
        priorU <- u^boundaryPower
        priorLogJacobian <- log(boundaryPower) +
            (boundaryPower - 1)*log(u)
    } else if (is.finite(region$upper)) {
        priorU <- 1 - (1 - u)^boundaryPower
        priorLogJacobian <- log(boundaryPower) +
            (boundaryPower - 1)*log1p(-u)
    }
    priorQuantiles <- if (is.infinite(region$lower)) {
        if (is.infinite(region$upper)) {
            stats::qt(p = priorU, df = pdf)
        } else {
            ## Compute log(priorU) without forming 1 - a tiny number.
            logPriorU <- log1p(-(1 - u)^boundaryPower)
            stats::qt(p = logPriorU + region$log_norm_const, df = pdf,
                      lower.tail = TRUE, log.p = TRUE)
        }
    } else {
        stats::qt(p = log1p(-priorU) + region$log_norm_const, df = pdf,
                  lower.tail = FALSE, log.p = TRUE)
    }
    d <- plocation + pscale*priorQuantiles

    vv <- rep(v, each = nquad)
    dd <- rep(d, times = nquad)
    ok <- is.finite(vv) & vv > 0 & is.finite(dd)
    if (!any(ok)) {
        return(NaN)
    }
    vv <- vv[ok]
    dd <- dd[ok]
    logWeights <- rep(log(w), each = nquad) +
        rep(log(w) + priorLogJacobian, times = nquad)
    logWeights <- logWeights[ok]

    y <- t*sqrt(vv/df)
    noncentrality <- eta*dd
    logRatio <- noncentrality*y - noncentrality^2/2
    logTerms <- logWeights + logRatio
    logTerms <- logTerms[is.finite(logTerms)]
    if (length(logTerms) == 0) {
        return(NaN)
    }

    -.bfpwr_logspace_sum(logTerms)
}

.tbf01_needs_exact_path <- function(t, alternative, log_bf) {
    ## Use the slower path only where the direct integral is unreliable. The
    ## cutoff is empirical: it catches the observed wrong-tail underflow cases
    ## without slowing down ordinary one-sided calculations.
    if (!is.finite(log_bf)) {
        return(TRUE)
    }
    if (alternative == "greater" && t <= -.tbf01_tail_quadrature_cutoff) {
        return(TRUE)
    }
    if (alternative == "less" && t >= .tbf01_tail_quadrature_cutoff) {
        return(TRUE)
    }
    FALSE
}

tbf01. <- function(t, n, n1 = n, n2 = n, plocation = 0, pscale = 1/sqrt(2),
                   pdf = 1, type = c("two.sample", "one.sample",  "paired"),
                   alternative = c("two.sided", "less", "greater"), log = FALSE,
                   tail.nquad = .bfpwr_defaults$tail.nquad,
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
        !is.na(log),

        .tbf01_valid_tail_nquad(tail.nquad)
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    if (type != "two.sample") {
        if (n1 != n2) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
        }
    }

    pars <- .tbf01_pars(n1 = n1, n2 = n2, type = type)
    region <- .tbf01_prior_region(plocation = plocation, pscale = pscale,
                                  pdf = pdf, alternative = alternative)

    ## Skip the direct integral where it is known to be unreliable, rather
    ## than spending its integration budget before using stable quadrature.
    log_bf <- NaN
    if (!.tbf01_needs_exact_path(t, alternative, log_bf = 0)) {
        log_bf <- .tbf01_log_fast(t = t, df = pars$df, neff = pars$neff,
                                  plocation = plocation, pscale = pscale,
                                  pdf = pdf, region = region, ...)
    }
    if (.tbf01_needs_exact_path(t = t, alternative = alternative,
                                log_bf = log_bf)) {
        log_bf <- .tbf01_log_tail_quadrature(
            t = t, df = pars$df, neff = pars$neff, plocation = plocation,
            pscale = pscale, pdf = pdf, region = region,
            tail.nquad = tail.nquad
        )
    }

    if (log) return(log_bf)
    else return(exp(log_bf))
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
#'     Most calculations use \code{stats::integrate}. Extreme wrong-tail
#'     one-sided calculations use fixed Gauss-Legendre quadrature to avoid
#'     noncentral-\eqn{t} underflow.
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
#' @param tail.nquad Number of Gauss-Legendre quadrature nodes used for stable
#'     wrong-tail one-sided calculations. Larger values are more accurate but
#'     slower. Defaults to \code{512} nodes per dimension (the fallback uses
#'     a two-dimensional product rule).
#' @param ... Additional arguments passed to \code{stats::integrate} for the
#'     direct one-dimensional Bayes factor integral. Defaults are
#'     \code{rel.tol = 1e-8}, \code{abs.tol = rel.tol}, and
#'     \code{subdivisions = 1000}.
#'
#' @inherit bf01 return
#'
#' @author Samuel Pawel, František Bartoš
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
