## Helper functions for sequential BF design calculations
## -----------------------------------------------------------------------------

#' @title Predictive Distribution Parameters
#'
#' @description Compute mean vector and covariance matrix and covariance of
#'     predictive distribution of z-statistics
#'
#' @param se Vector of standard errors
#' @param dpm Design prior mean
#' @param dpsd Design prior standard deviation
#'
#' @return A list with the predictive mean vector and covariance matrix
#'
#' @author Samuel Pawel
#'
#' @noRd
#'
#' @keywords internal
#'
#' @examples
#' ## regions to stop with two-sided p < 0.05 in first or second stage
#' predpars(se = sqrt(2/seq(10, 50, 10)), dpm = 0.5, dpsd = 0.1)
predpars <- function(se, dpm, dpsd) {
    m <- length(se)
    inf <- 1/se^2 # information levels
    mean <- dpm/se # mean vector
    sigma <- matrix(nrow = m, ncol = m)
    for (i in seq_len(m)) {
        for (j in seq_len(m)) {
            sigma[i,j] <- sqrt(pmin(inf[i], inf[j])/pmax(inf[i], inf[j]))
        }
    }
    sigma <- sigma + dpsd^2 * sqrt(inf) %*% t(sqrt(inf))
    list("mean" = mean, "sigma" = sigma)
}

#' @title Integrate Success Regions of Cumulative Z-statistics
#'
#' @description Compute per-stage probabilities that cumulative z-statistics
#'     fall into given success regions
#'
#' @param intregions A list of length `m = length(mean)`. Each element `i` is a
#'     list of regions for stage `i`. Each region is a 2 x i matrix, where the
#'     first row gives lower bounds and the second row gives upper bounds. If
#'     any element of a region is `NaN`, that region's contribution is treated
#'     as zero
#' @param mean Numeric vector of means for the cumulative z-statistic
#' @param sigma Covariance matrix of the cumulative z-statistic
#' @param method Method to compute the integral. Either \code{lpmvnorm}
#'     (default) or \code{"pmvnorm"}
#' @param ... Other arguments passed to \code{mvtnorm::lpmvnorm} or
#'     \code{mvtnorm::pmvnorm}
#'
#' @return Numeric vector of per-stage probabilities
#'
#' @author Samuel Pawel
#'
#' @noRd
#'
#' @keywords internal
#'
#' @examples
#' ## regions to stop with two-sided p < 0.05 in first or second stage
#' intregions <- list(list(matrix(c(1.96,
#'                                  Inf), byrow = TRUE, ncol = 1),
#'                         matrix(c(-Inf,
#'                                  -1.96), byrow = TRUE, ncol = 1)),
#'                    list(matrix(c(-1.96, 1.96,
#'                                  1.96, Inf), byrow = TRUE, ncol = 2),
#'                         matrix(c(-1.96, -Inf,
#'                                  1.96, -1.96), byrow = TRUE, ncol = 2)))
#' n1 <- 50
#' n2 <- 100
#' mu <- 0.2
#' mean <- mu*sqrt(c(n1, n2))
#' sigma <- matrix(c(1, sqrt(n1/n2),
#'                   sqrt(n1/n2), 1), nrow = 2, byrow = TRUE)
#' intstages(intregions = intregions, mean = mean, sigma = sigma)

intstages <- function(intregions, mean, sigma, method = "lpmvnorm", ...) {
    stopifnot(
        is.list(intregions),
        is.numeric(mean),
        is.matrix(sigma),
        length(mean) == nrow(sigma),
        nrow(sigma) == ncol(sigma)
    )

    m <- length(mean)
    probs <- numeric(m)

    if (m > 1 & method == "lpmvnorm") {
        ## for more stable multivariate normal integration
        C <- t(chol(sigma))
        Ct <- mvtnorm::ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE)
        ## use a fixed grid instead of Monte Carlo approach
        ngrid <- 1000
        w <- withr::with_seed(seed = 42, code = {
            t(qrng::ghalton(n = ngrid, d = m - 1))
        })
    }

    for (i in seq_len(m)) {
        stageregions <- intregions[[i]]
        stopifnot(is.list(stageregions))

        regionprobs <- vapply(stageregions,
                              FUN.VALUE = numeric(1),
                              FUN = function(region) {
            ## NaN encodes that critical value doesn't exist => probability = 0
            if (any(is.nan(region))) {
                p <- 0
            } else {
                if (i == 1) {
                    p <- exp(.bfpwr_lpnorm_interval(lower = region[1,],
                                                    upper = region[2,],
                                                    mean = mean[1],
                                                    sd = sqrt(sigma[1:1])))
                } else if (method == "lpmvnorm") {
                    p <- exp(mvtnorm::lpmvnorm(lower = region[1, ],
                                               upper = region[2, ],
                                               mean = mean[1:i],
                                               chol = Ct[,1:i],
                                               M = ngrid,
                                               w = w[1:(i - 1),,drop = FALSE],
                                               ...))
                } else {
                    p <- mvtnorm::pmvnorm(lower = region[1, ],
                                          upper = region[2, ],
                                          mean  = mean[1:i],
                                          sigma = sigma[1:i, 1:i],
                                          seed = 42,
                                          keepAttr = FALSE,
                                          ...)
                }
            }
            return(p)
        })
        probs[i] <- sum(regionprobs, na.rm = TRUE)
    }

    return(probs)
}


#' @title Generate Integration Regions for Cumulative Z-Statistics Based on One
#'     Critical Value
#'
#' @description Constructs the per-stage integration regions for cumulative
#'     z-statistics given one lower and one upper critical value.
#'
#' @param zcrit0 Numeric vector of lower critical values (H0 boundaries).
#' @param zcrit1 Numeric vector of upper critical values (H1 boundaries).
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{H1}{List of length `m`; per-stage regions where evidence supports H1.}
#'   \item{H0}{List of length `m`; per-stage regions where evidence supports H0.}
#' }
#'
#' Each element is a list containing one 2 x i matrix representing lower and upper
#' bounds for the first `i` stages.
#'
#' @noRd
#'
#' @keywords internal
#'
#' @examples
#' zcrit0 <- c(-1.2, -1.1, -1)
#' zcrit1 <- c(2,  1.96, 1.9)
#' genregions1(zcrit0, zcrit1)
#'

genregions1 <- function(zcrit0, zcrit1) {
    stopifnot(all(is.numeric(zcrit0)),
              all(is.numeric(zcrit1)),
              length(zcrit0) == length(zcrit1))

    H0nan <- is.nan(zcrit0)
    finite <- !H0nan & !is.nan(zcrit1)
    finiteH1 <- !is.nan(zcrit1)
    if (any(finite) && all(zcrit1[finite] >= zcrit0[finite])) {
        direction <- "positive"
    } else if (any(finite) && all(zcrit1[finite] < zcrit0[finite])) {
        direction <- "negative"
    } else if (!any(finite) && any(finiteH1) && all(zcrit1[finiteH1] >= 0)) {
        direction <- "positive"
    } else if (!any(finite) && any(finiteH1) && all(zcrit1[finiteH1] <= 0)) {
        direction <- "negative"
    } else {
        stop("Inconsistent critical values: direction cannot be inferred.")
    }

    m <- length(zcrit1)
    intregionsH1 <- vector("list", m)
    intregionsH0 <- vector("list", m)

    for (i in seq_len(m)) {
        ## region where evidence for H1 in stage i
        matH1 <- matrix(nrow = 2, ncol = i)
        for (j in seq_len(i)) {
            if (i == j) {
                ## evidence for H1
                if (direction == "positive") {
                    lower <- zcrit1[j]
                    upper <- Inf
                } else {
                    lower <- -Inf
                    upper <- zcrit1[j]
                }
            } else {
                ## continue (no stop yet)
                if (direction == "positive") {
                    if (H0nan[j]) lower <- -Inf
                    else lower <- zcrit0[j]
                    upper <- zcrit1[j]
                } else {
                    lower <- zcrit1[j]
                    if (H0nan[j]) upper <- Inf
                    else upper <- zcrit0[j]
                }
            }
            matH1[, j] <- c(lower, upper)
        }
        intregionsH1[[i]] <- list(matH1)

        ## region where evidence for H0 in stage i
        matH0 <- matH1
        if (direction == "positive") {
            matH0[, i] <- c(-Inf, zcrit0[i])
        } else {
            matH0[, i] <- c(zcrit0[i], Inf)
        }
        intregionsH0[[i]] <- list(matH0)
    }

    return(list(H1 = intregionsH1, H0 = intregionsH0))
}


#' @title Generate Integration Regions for Cumulative Z-Statistics Based on Two
#'     Critical Values
#'
#' @description Constructs per-stage integration regions for cumulative
#'     z-statistics given lower and upper critical boundaries. This function is
#'     intended for situations with two success regions for H1 (both tails) and
#'     corresponding stopping regions for H0.
#'
#' @param zcrit0 2 x m numeric matrix of H0 boundaries. Each column corresponds
#'     to one stage. The first row gives the lower bound and the second row the
#'     upper bound for H0. Specify NaN if no H0 boundary exists at a stage
#' @param zcrit1 2 x m numeric matrix of H1 boundaries. Each column corresponds
#'     to one stage. The first row gives the upper bound of the lower region
#'     (extending from -Inf to this upper boudn) and the second row the lower
#'     bound of the upper region for H1 (extending from this lower bound to Inf)
#' @param strict Logical. If \code{TRUE}, return all possible region
#'     combinations (slow but exact). If \code{FALSE}, only returns the main
#'     regions where the sign of the z-statistics does not change across stages
#'     (faster, recommended when many interim analyses are performed). Defaults
#'     to \code{FALSE}
#'
#' @return
#' A list with two components:
#' \describe{
#'   \item{H1}{List of length `m`; per-stage regions where evidence supports H1.}
#'   \item{H0}{List of length `m`; per-stage regions where evidence supports H0.}
#' }
#'
#' Each element is a list of 2 x i matrices, representing lower and upper bounds
#' for the first `i` stages.
#'
#' @noRd
#'
#' @keywords internal
#'
#' @examples
#' zcrit0 <- matrix(c(-1, -1, -1,
#'                     1,  1, 1),
#'                  nrow = 2, byrow = TRUE)
#' zcrit1 <- matrix(c(-1.96, -1.96, -1.96,
#'                     1.96,  1.96, 1.96),
#'                  nrow = 2, byrow = TRUE)
#' genregions2(zcrit0, zcrit1)
#'

genregions2 <- function(zcrit0, zcrit1, strict = FALSE) {
    stopifnot(
        is.matrix(zcrit0),
        is.matrix(zcrit1),
        all(dim(zcrit0) == dim(zcrit1)),
        nrow(zcrit0) == 2
    )

    ## if (strict == TRUE & ncol(zcrit0) > 10) {
    ##     warning("strict = TRUE with many stages may cause numerical problems")
    ## }

    m <- ncol(zcrit0)
    intregionsH1 <- vector("list", m)
    intregionsH0 <- vector("list", m)

    ## identify stages where evidence for H0 impossible
    H0nan <- apply(zcrit0, 2, function(x) any(is.nan(x)))

    for (i in seq_len(m)) {
        ## build stage-wise regions
        H1i <- lapply(seq_len(i), function(j) {
            if (i == j) {
                ## stopping regions (evidence for H1)
                list(
                    c(-Inf, zcrit1[1, j]), # lower
                    c(zcrit1[2, j], Inf)   # upper
                )
            } else {
                ## continuation regions (no evidence for H1 or H0)
                if (H0nan[j] == TRUE) {
                    list(c(zcrit1[1, j], zcrit1[2, j]))
                } else {
                    list(
                        ## between H1 lower and H0 lower
                        c(zcrit1[1, j], zcrit0[1, j]),
                        ## between H0 upper and H1 upper
                        c(zcrit0[2, j], zcrit1[2, j])
                    )
                }
            }
        })

        H0i <- H1i
        H0i[[i]] <- list(c(zcrit0[1, i], zcrit0[2, i])) # H0 stop region

       ## build all region combinations
        if (strict == TRUE) {
            combosH1i <- expand.grid(lapply(H1i, seq_along))
            combosH0i <- expand.grid(lapply(H0i, seq_along))
        } else {
            ## integrating all regions is usually not worth it because the
            ## probability of regions where the sign of z_i flips is almost
            ## zero), hence, only take the two regions where the sign of z_i
            ## doesn't flip (the first and last)
            combosH1i <- rbind(rep(1, length(H1i)),
                               sapply(H1i, length))
            if (i <= 1 + sum(H0nan)) {
                ## only one H0 region in the first non-NaN stage
                combosH0i <- matrix(rep(1, length(H0i)), nrow = 1)
            } else {
                combosH0i <- rbind(rep(1, length(H0i)),
                                   sapply(H0i, length))
            }
        }

        makeregions <- function(combos, Hi) {
            apply(combos, 1, function(row) {
                sapply(seq_along(row), function(j) {
                    region_index <- as.numeric(row[[j]])
                    Hi[[j]][[region_index]]
                })
            }, simplify = FALSE)
        }

        intregionsH1[[i]] <- makeregions(combosH1i, H1i)
        intregionsH0[[i]] <- makeregions(combosH0i, H0i)
    }

    list(H1 = intregionsH1, H0 = intregionsH0)
}

#' @title Compute Critical Z-Values for Bayes Factors
#'
#' @description Computes critical z-values for Bayes factors using normal,
#'     directional normal, or normal moment priors.
#'
#' @param k Positive numeric. Bayes factor threshold (BF01 oriented in favor of
#'     H0)
#' @param se Positive numeric. Standard error
#' @param mu Numeric. Prior location. Not taken into account for \code{type =
#'     "moment"}
#' @param tau Non-negative numeric. Prior scale
#' @param type Character. One of "normal" (point null vs. normal alternative),
#'     "directional" (directional null vs. directional alternative with marginal
#'     normal prior), or "moment" (point null vs. normal moment alternative)
#'
#' @return Numeric vector of critical z-value(s)
#'
#' @examples
#' se <- 0.05
#' pm <- 0.2
#' psd <- 0.1
#' zcrit1 <- zcrit(k = 3, se = se, mu = pm, tau = psd, type = "normal")
#' bf01(estimate = zcrit1*se, se = se, null = 0, pm = pm, psd = psd)
#'
#' ## tau = 0 leads to point alternative
#' zcrit2 <- zcrit(k = 5, se = se, mu = pm, tau = 0, type = "normal")
#' bf01(estimate = zcrit2*se, se = se, null = 0, pm = pm, psd = 0)
#'
#' zcrit3 <- zcrit(k = 5, se = se, mu = pm, tau = psd, type = "directional")
#' dirbf01(estimate = zcrit3*se, se = se, null = 0, pm = pm, psd = psd)
#'
#' zcrit4 <- zcrit(k = 1/10, se = se, mu = 0, tau = psd, type = "moment")
#' nmbf01(estimate = zcrit4*se, se = se, null = 0, psd = psd)
#'
#' @noRd
#'
#' @keywords internal
zcrit <- function(k, se, mu = NULL, tau, type = c("normal", "directional", "moment")) {

    type <- match.arg(type)

    if (type == "normal") {
        if (tau == 0) {
            ## point prior under the alternative
            zcrit <- (mu^2/se^2 - 2*log(k))/(2*mu/se)
        } else {
            ## normal prior under the alternative
            X <- (mu^2/tau^2 + log(1 + tau^2/se^2) - 2*log(k))*
                (1 + se^2/tau^2)
            if (X < 0) {
                zcrit <- c(NaN, NaN)
            } else {
                M <- -mu*se/tau^2
                zcrit <- M + c(-1, 1)*sqrt(X)
            }
        }
    }

    if (type == "directional") {
        logpriorodds <- stats::pnorm(q = mu/tau, lower.tail = FALSE,
                                     log.p = TRUE) -
            stats::pnorm(q = mu/tau, lower.tail = TRUE, log.p = TRUE)
        postq <- .bfpwr_qnorm_logistic_inverse(log(k) + logpriorodds)
        zcrit <- (postq*sqrt(1/se^2 + 1/tau^2) -
                  mu/tau^2)*se
    }

    if (type == "moment") {
        Y <- (2*lamW::lambertW0(x = exp(0.5)*(1 + tau^2/se^2)^1.5/(2*k)) -
              1)*(1 + se^2/tau^2)
        if (Y < 0) {
            zcrit <- c(NaN, NaN)
        } else {
            zcrit <- c(-1, 1)*sqrt(Y)
        }
    }

    return(zcrit)
}



#' @title Compute Critical T-Values for T-Test Bayes Factors
#'
#' @description Computes critical t-values for T-Test Bayes factors
#'
#' @param k Positive numeric. Bayes factor threshold (BF01 oriented in favor of
#'     H0)
#' @param n1 Sample size in group 1
#' @param n2 Sample size in group 2 (is ignored for one-sample \eqn{t}-tests)
#' @param plocation \eqn{t} prior location
#' @param pscale \eqn{t} prior scale
#' @param pdf \eqn{t} prior degrees of freedom
#' @param type Type of \eqn{t}-test. Can be \code{"two.sample"},
#'     \code{"one.sample"}, or \code{"paired"}
#' @param alternative Direction of the test. Can be either \code{"two.sided"},
#'     \code{"less"}, or \code{"greater"}. The latter two truncate the analysis
#'     prior to negative and positive effects, respectively
#' @param drange Numerical search strategy. Can be either \code{"adaptive"}
#'     (default) or an interval
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return Numeric vector of critical t-value(s)
#'
#' @examples
#' tseq <- seq(-10, 10, length.out = 100)
#' n1 <- 50
#' n2 <- 60
#' type <- "two.sample"
#' alternative <- "two.sided"
#' plocation <- 0
#' pscale <- 1/sqrt(2)
#' pdf <- 1
#' k <- 3
#' tcrit1 <- tcrit(k = k, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
#'                 pdf = pdf, alternative = alternative, type = type)
#' plot(tseq, tbf01(t = tseq, n1 = n1, n2 = n2, plocation = plocation,
#'                  pscale = pscale, pdf = pdf, alternative = alternative,
#'                  type = type),
#'      type = "l", xlab = "t-statistic", ylab = bquote("BF"["01"]), log = "y")
#' abline(h = k, lty = 2)
#' abline(v = tcrit1, lty = 2)
#'
#' n1 <- n2 <- 100
#' alternative <- "greater"
#' k <- 6
#' tcrit2 <- tcrit(k = k, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
#'                 pdf = pdf, alternative = alternative, type = type)
#' plot(tseq, tbf01(t = tseq, n1 = n1, n2 = n2, plocation = plocation,
#'                  pscale = pscale, pdf = pdf, alternative = alternative,
#'                  type = type),
#'      type = "l", xlab = "t-statistic", ylab = bquote("BF"["01"]), log = "y")
#' abline(h = k, lty = 2)
#' abline(v = tcrit2, lty = 2)
#'
#' @noRd
#'
#' @keywords internal
tcrit <- function(k, n1, n2, plocation, pscale, pdf, type, alternative,
                  drange = "adaptive", ...) {

    ## determine t-statistic for which BF = k
    rootFun <- function(t) {
        tbf01(t = t, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
              pdf = pdf, type = type, alternative = alternative,
              log = TRUE) - log(k)
    }

    if (k > 1) {
        ## find maximum BF to see whether BF = k is possible
        opt <- stats::optim(par = 0, fn = rootFun, control = list(fnscale = -1),
                            method = "BFGS")
        if (opt$convergence != 0) {
            warning("numerical problems finding maximum BF")
            return(NaN)
        } else {
            if (opt$value < 0) {
                warning("maximum BF is less than k; BF01 = k impossible")
                if (alternative == "two.sided") {
                    return(c(NaN, NaN))
                } else {
                    return(NaN)
                }
            }
        }
    }

    if (alternative == "two.sided") {
        ## guess search range based on search range from z-test BF
        if (!is.numeric(drange) && drange == "adaptive") {
            if (type == "two.sample") {
                neff <- 1/(1/n1 + 1/n2)
            } else {
                neff <- n1
            }
            se <- 1/sqrt(neff)
            X <- (plocation^2/pscale^2 + log(1 + pscale^2/se^2) -
                  2*log(k))*(1 + se^2/pscale^2)
            if (X <= 0) {
                X <- 5/k
            }
            zcrit <- -plocation*se/pscale^2 + c(-1, 1)*sqrt(X)
            meant <- mean(zcrit)
            if (zcrit[1] < zcrit[2]) {
                searchIntLow <- c(zcrit[1] - 2, meant)
                searchIntUp <- c(meant, zcrit[2] + 2)
            } else {
                searchIntLow <- c(zcrit[2] - 2, meant)
                searchIntUp <- c(meant, zcrit[1] + 2)
            }
        } else {
            meant <- mean(drange)
            searchIntLow <- c(drange[1], meant)
            searchIntUp <- c(meant, drange[2])
        }
        ## search for critical values
        tcrit <- c(NaN, NaN)
        lower <- try(stats::uniroot(f = rootFun, interval = searchIntLow,
                                    extendInt = "upX", ...)$root,
                     silent = TRUE)
        upper <- try(stats::uniroot(f = rootFun, interval = searchIntUp,
                                    extendInt = "downX", ...)$root,
                     silent = TRUE)
        if (inherits(lower, "try-error") || inherits(upper, "try-error")) {
            warning("Numerical problems: Could not find 2 t-roots")
        } else {
            tcrit <- c(lower, upper)
        }
    } else { # one-sided cases
        if (!is.numeric(drange) && drange == "adaptive") {
            ## extend the search range if critical value not contained
            if (alternative == "greater") {
                ## want to first find the critical value on the positive side
                searchint <- c(0, 0.1)
                extend <- "downX"
            } else {
                ## want to first find the critical value on the negative side
                searchint <- c(-0.1, 0)
                extend <- "upX"
            }
        } else {
            searchint <- drange
            extend <- "no"

        }
        suppressWarnings({
            res <- try(stats::uniroot(f = rootFun, interval = searchint,
                                      extendInt = extend, ...)$root,
                       silent = TRUE)
        })
        if (inherits(res, "try-error")) {
            warning("Numerical problems finding critical value")
            tcrit <- NaN
        } else {
            tcrit <- res
        }
    }
    return(tcrit)
}
