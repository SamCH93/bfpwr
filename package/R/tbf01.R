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

## Gauss-Laguerre nodes for the opposite-sign noncentral-t density; 64 nodes
## agree with adaptive integration to about 1e-12 on the log scale for
## df from 1 to 2e4, |t| up to 40 and |ncp| up to 100.
.tbf01_dnt_nquad <- 64

.tbf01_log_dnt_opposite <- function(t, df, ncp) {
    ## Exact noncentral-t log density at t for a vector of ncp values whose
    ## sign is opposite to the sign of t. Write
    ##   dnt(t; df, ncp) = dt(t; df) exp(-ncp^2/2) E[exp(s Y)],
    ##   Y ~ chi with df + 1 degrees of freedom, s = ncp t/sqrt(df + t^2) < 0.
    ## With k = df + 1, E[exp(s Y)] = I(s)/I(0) for
    ## I(s) = integral_0^Inf y^(k - 1) exp(-y^2/2 + s y) dy, where
    ## I(0) = 2^(k/2 - 1) Gamma(k/2). Writing
    ## -y^2/2 = a^2/2 - a y - (y - a)^2/2 with a (a - s) = k gives
    ##   I(s) = exp(a^2/2) Gamma(k) (a - s)^(-k) E[exp(-(G - a)^2/2)],
    ## G ~ Gamma(k, rate = a - s). G has mean a and, for s < 0, standard
    ## deviation a/sqrt(k) < 1, so the remaining bounded factor is smooth and
    ## generalized Gauss-Laguerre quadrature evaluates it accurately.
    ## t/sqrt(df + t^2) is written as sign(t)/sqrt(1 + df/t^2) so that t^2
    ## cannot overflow.
    s <- ncp*sign(t)/sqrt(1 + df/t^2)
    k <- df + 1
    a <- (s + sqrt(s^2 + 4*k))/2
    rule <- .bfpwr_gauss_laguerre(shape = k, n = .tbf01_dnt_nquad)
    G <- outer(1/(a - s), rule$x)
    log_correction <- log(drop(exp(-(G - a)^2/2) %*% rule$w))
    log_I0 <- (k/2 - 1)*log(2) + lgamma(k/2)
    stats::dt(x = t, df = df, log = TRUE) - ncp^2/2 + a^2/2 + lgamma(k) -
        k*log(a - s) - log_I0 + log_correction
}

## Gauss-Hermite nodes and series terms for the same-sign noncentral-t density
## (see .tbf01_log_dnt_same()).
.tbf01_dnt_nhermite <- 64
.tbf01_dnt_nseries <- 400

.tbf01_log_dnt_same_hermite <- function(t, df, ncp, s) {
    ## Same-sign log density (s = ncp sign(t)/sqrt(1 + df/t^2) > 0) from the
    ## chi mixture of .tbf01_log_dnt_opposite(). Its Gauss-Laguerre
    ## evaluation relies on s < 0 (G concentrated), so write
    ##   I(s) = integral_0^Inf y^(k - 1) exp(-y^2/2 + s y) dy
    ##        = exp(s^2/2) J(s),  J(s) = integral_0^Inf y^(k - 1) exp(-(y - s)^2/2) dy,
    ## so that, with k = df + 1 and I(0) = 2^(k/2 - 1) Gamma(k/2),
    ##   log dnt = log dt(t; df) - ncp^2 df/(2 (df + t^2)) + log J(s) - log I(0)
    ## (-ncp^2/2 + s^2/2 = -ncp^2 df/(2 (df + t^2))). In z = log y the
    ## integrand exp(k z - (e^z - s)^2/2) is unimodal with its Laplace mode at
    ## ys = (s + sqrt(s^2 + 4k))/2 (where y (y - s) = k) and curvature
    ## ys sqrt(s^2 + 4k) there, so a Gauss-Hermite rule centred at log(ys)
    ## with scale sig = 1/sqrt(ys sqrt(s^2 + 4k)) integrates it accurately,
    ## for large s and for large k (relative error below 1e-9 on the log
    ## scale, checked against adaptive integration, also at log densities of
    ## order -1e11). Every term is written to stay finite for t up to the
    ## largest double: s is only squared in 4k/s^2 when s is large, and
    ## ncp^2 df/(2 (df + t^2)) is evaluated as (ncp/t)^2 df/(2 (df/t^2 + 1))
    ## for |t| > 1.
    k <- df + 1
    big <- s > sqrt(k)
    root <- ys <- numeric(length(s))
    q <- sqrt(1 + 4*k/s[big]^2)
    root[big] <- s[big]*q
    ys[big] <- s[big]*((1 + q)/2)
    root[!big] <- sqrt(s[!big]^2 + 4*k)
    ys[!big] <- (s[!big] + root[!big])/2
    ## ys - s without cancellation
    shift <- 2*k/(root + s)
    sig <- 1/(sqrt(ys)*sqrt(root))

    rule <- .bfpwr_gauss_hermite(.tbf01_dnt_nhermite)
    delta <- outer(sqrt(2)*sig, rule$x)
    h <- k*(log(ys) + delta) - (ys*expm1(delta) + shift)^2/2
    terms <- sweep(h, 2, log(rule$w) + rule$x^2, "+")
    top <- apply(terms, 1, max)
    log_J <- top + log(rowSums(exp(terms - top))) + log(sqrt(2)*sig)

    log_quad <- if (abs(t) <= 1) {
        ncp^2*df/(2*(df + t^2))
    } else {
        (ncp/t)^2*df/(2*(df/t^2 + 1))
    }
    stats::dt(x = t, df = df, log = TRUE) - log_quad + log_J -
        ((k/2 - 1)*log(2) + lgamma(k/2))
}

.tbf01_log_dnt_same_series <- function(t, df, ncp, s) {
    ## Same-sign log density for small k = df + 1 and 0 < s < 10, where the
    ## Gauss-Hermite rule of .tbf01_log_dnt_same_hermite() is not accurate
    ## (the integrand is not yet close to normal in log y). The even part of
    ## the density is a confluent hypergeometric function,
    ##   dnt(t) + dnt(-t) = 2 dt(t; df) exp(-ncp^2/2) 1F1(k/2; 1/2; s^2/2),
    ## with 1F1 = sum_j (k/2)_j/(1/2)_j (s^2/2)^j/j!, a series of positive
    ## terms. For k/2 < 5 and s^2/2 < 50 the terms peak below j = 55 and are
    ## negligible by j = 400. The mirrored opposite-sign density dnt(-t) is
    ## subtracted from the even part, and since dnt(-t) <= dnt(t) at most a
    ## factor two cancels.
    k <- df + 1
    j <- seq_len(.tbf01_dnt_nseries + 1) - 1
    base <- lgamma(k/2 + j) - lgamma(k/2) + lgamma(1/2) - lgamma(1/2 + j) -
        lgamma(j + 1)
    terms <- outer(log(s^2/2), j) + rep(base, each = length(s))
    terms[, 1] <- base[1]
    top <- apply(terms, 1, max)
    log_1F1 <- top + log(rowSums(exp(terms - top)))
    log_even <- log(2) + stats::dt(x = t, df = df, log = TRUE) - ncp^2/2 +
        log_1F1
    log_opp <- .tbf01_log_dnt_opposite(t = -t, df = df, ncp = ncp)
    log_even + log1p(-exp(log_opp - log_even))
}

.tbf01_log_dnt <- function(t, df, ncp) {
    ## Noncentral-t log density of the statistic t for a vector of ncp values.
    ## stats::dt is unreliable in the tails, where it subtracts nearly equal
    ## distribution functions and returns noise (t = -4.9, df = 24, ncp = 6
    ## gives log density -29.2 instead of -44.8; the same-sign density is off
    ## by up to 650 for extreme t and very large samples, and for df above
    ## about 4e5). Noisy integrands make stats::integrate fail or miss the
    ## likelihood and send informed priors and large samples to a fallback
    ## that is only accurate in the one-sided tail. The density therefore
    ## comes from exact representations, according to the sign of
    ## s = ncp sign(t)/sqrt(1 + df/t^2) (that of ncp t):
    ## - s < 0 (t and ncp of opposite sign): the chi mixture of
    ##   .tbf01_log_dnt_opposite(), evaluated with Gauss-Laguerre quadrature.
    ## - s > 0 with s >= 10 or df + 1 >= 10: the same mixture, evaluated with
    ##   Gauss-Hermite quadrature in log y, .tbf01_log_dnt_same_hermite().
    ## - s > 0 with s < 10 and df + 1 < 10: the hypergeometric series of the
    ##   even part, .tbf01_log_dnt_same_series().
    ## - s = 0 (t = 0 or ncp = 0): the central case, stats::dt is accurate.
    ## (The noncentral F density of t^2 would give the even part as well, but
    ## stats::df returns NaN once ncp^2 t^2/2 exceeds 2^31 and is slow for
    ## large noncentrality.)
    s <- ncp*sign(t)/sqrt(1 + df/t^2)
    out <- numeric(length(ncp))
    opposite <- is.finite(s) & s < 0
    same <- is.finite(s) & s > 0
    other <- !(opposite | same)
    if (any(other)) {
        out[other] <- suppressWarnings(
            stats::dt(x = t, df = df, ncp = ncp[other], log = TRUE)
        )
    }
    if (any(opposite)) {
        out[opposite] <- .tbf01_log_dnt_opposite(t = t, df = df,
                                                 ncp = ncp[opposite])
    }
    if (any(same)) {
        hermite <- same & (s >= 10 | df + 1 >= 10)
        series <- same & !hermite
        if (any(hermite)) {
            out[hermite] <- .tbf01_log_dnt_same_hermite(
                t = t, df = df, ncp = ncp[hermite], s = s[hermite]
            )
        }
        if (any(series)) {
            out[series] <- .tbf01_log_dnt_same_series(
                t = t, df = df, ncp = ncp[series], s = s[series]
            )
        }
    }
    out
}

.tbf01_log_fast <- function(t, df, neff, plocation, pscale, pdf, region,
                            rel.tol = getOption("bfpwr.rel.tol", 1e-8), abs.tol = rel.tol,
                            subdivisions = getOption("bfpwr.subdivisions", 1000), ...) {
    ## Original one-dimensional integral, evaluated on a centered log scale.
    if (!is.finite(region$log_norm_const)) {
        return(NaN)
    }

    ## The noncentral-t integral is unreliable in this one-sided tail. Reject
    ## it as a scout so callers use the stable quadrature path directly.
    if ((region$lower == 0 && t <= -.tbf01_tail_quadrature_cutoff) ||
        (region$upper == 0 && t >= .tbf01_tail_quadrature_cutoff)) {
        return(NaN)
    }

    eta <- sqrt(neff)
    log_f0 <- stats::dt(x = t, df = df, log = TRUE)

    ## Where the density cannot be evaluated even at the likelihood peak
    ## ncp = t (it is not finite there), leave the Bayes factor to the
    ## fallback.
    if (!is.finite(.tbf01_log_dnt(t = t, df = df, ncp = t))) {
        return(NaN)
    }

    log_prior <- function(d) {
        stats::dt(x = (d - plocation)/pscale, df = pdf, log = TRUE) -
            log(pscale) - region$log_norm_const
    }
    log_integrand <- function(d) {
        suppressWarnings({
            .tbf01_log_dnt(t = t, df = df, ncp = eta*d) + log_prior(d)
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

    ## A single integrate() over an infinite range can miss a likelihood spike
    ## that is narrow at large n (n = 1e7, t = 20 gives log BF01 -106.5
    ## instead of -192.06; splitting only at the peak still loses one half of
    ## the spike). The likelihood in d is roughly normal around t/eta with
    ## standard deviation width, so integrate piecewise between the region
    ## limits, the null, and +-8 width around the likelihood peak and +-8
    ## prior scales around the prior location.
    width <- sqrt(1 + t^2/(2*df))/eta
    breaks <- c(region$lower, region$upper, 0,
                t/eta + c(-8, 0, 8)*width, plocation + c(-8, 0, 8)*pscale)
    breaks <- breaks[!is.na(breaks)]
    breaks <- sort(unique(breaks[breaks >= region$lower &
                                 breaks <= region$upper]))
    pieces <- vapply(seq_len(length(breaks) - 1), function(i) {
        value <- try(stats::integrate(f = intfun, lower = breaks[i],
                                      upper = breaks[i + 1], rel.tol = rel.tol,
                                      abs.tol = abs.tol,
                                      subdivisions = subdivisions,
                                      ...)$value,
                     silent = TRUE)
        if (inherits(value, "try-error")) NaN else value
    }, numeric(1))
    f1 <- sum(pieces)
    if (!all(is.finite(pieces)) || f1 <= 0) {
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

.tbf01_wrong_tail <- function(t, alternative) {
    ## Use the slower path only where the direct integral is unreliable. The
    ## cutoff is empirical: it catches the observed wrong-tail underflow cases
    ## without slowing down ordinary one-sided calculations.
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
                   tail.nquad = getOption("bfpwr.tail.nquad", 512),
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
    .bfpwr_validate_controls(list(...), root = FALSE)
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
    if (!.tbf01_wrong_tail(t, alternative)) {
        log_bf <- .tbf01_log_fast(t = t, df = pars$df, neff = pars$neff,
                                  plocation = plocation, pscale = pscale,
                                  pdf = pdf, region = region, ...)
    }
    if (!is.finite(log_bf)) {
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
#'     Most calculations use \code{stats::integrate}, with the integral over
#'     the prior split around the likelihood peak and the prior location.
#'     Noncentral-\eqn{t} densities come from exact representations instead of
#'     \code{stats::dt}, which is unreliable in the tails: a chi mixture
#'     evaluated with Gauss-Laguerre quadrature when the \eqn{t}-statistic and
#'     the noncentrality parameter have opposite signs, and Gauss-Hermite
#'     quadrature or a hypergeometric series otherwise. Extreme wrong-tail
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
#'     Session defaults can be changed with \link{bfpwrOptions}.
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
