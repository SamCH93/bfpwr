## Helper functions for sequential BF design calculations
## -----------------------------------------------------------------------------

#' @title Predictive Distribution Parameters
#'
#' @description Compute mean vector and covariance matrix and covariance of
#'     predictive distribution of z-statistics
#'
#' @param se Vector of standard errors
#' @param null Parameter value under the point null hypothesis. Defaults to
#'     \code{0}
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
predpars <- function(se, null = 0, dpm, dpsd) {
    m <- length(se)
    inf <- 1/se^2 # information levels
    mean <- (dpm - null)/se # mean vector
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
#' @author Samuel Pawel, František Bartoš
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
        if (i > 1 && method == "lpmvnorm") {
            probs[i] <- .bfseq_intstage_sum(
                stageregions = intregions[[i]], mean = mean[1:i],
                sigma = sigma[1:i, 1:i], method = method,
                cholFactor = Ct[,1:i], w = w[1:(i - 1),,drop = FALSE],
                ngrid = ngrid, ...
            )
        } else {
            probs[i] <- .bfseq_intstage_sum(
                stageregions = intregions[[i]], mean = mean[1:i],
                sigma = sigma[1:i, 1:i], method = method, ...
            )
        }
    }

    return(probs)
}

## Integrate the stopping regions for one terminal stage, preparing the
## quasi-Monte Carlo grid when lpmvnorm is used.
.bfseq_intstage <- function(stageregions, mean, sigma, method = "lpmvnorm",
                            ...) {
    stopifnot(
        is.list(stageregions),
        is.numeric(mean),
        is.matrix(sigma),
        length(mean) == nrow(sigma),
        nrow(sigma) == ncol(sigma)
    )

    i <- length(mean)
    if (i > 1 && method == "lpmvnorm") {
        C <- t(chol(sigma))
        Ct <- mvtnorm::ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE)
        ngrid <- 1000
        w <- withr::with_seed(seed = 42, code = {
            t(qrng::ghalton(n = ngrid, d = i - 1))
        })
        return(.bfseq_intstage_sum(stageregions = stageregions,
                                   mean = mean, sigma = sigma,
                                   method = method, cholFactor = Ct,
                                   w = w, ngrid = ngrid, ...))
    }

    .bfseq_intstage_sum(stageregions = stageregions, mean = mean,
                        sigma = sigma, method = method, ...)
}

## Sum the probability mass over all disjoint stopping regions for one stage.
.bfseq_intstage_sum <- function(stageregions, mean, sigma,
                                method = "lpmvnorm", cholFactor = NULL,
                                w = NULL, ngrid = 1000, ...) {
    stopifnot(is.list(stageregions))

    i <- length(mean)
    regionprobs <- vapply(stageregions,
                          FUN.VALUE = numeric(1),
                          FUN = function(region) {
        ## NaN in generated regions encodes an empty stopping event.
        if (any(is.nan(region))) {
            p <- 0
        } else if (any(is.na(region))) {
            stop("Sequential integration received NA bounds for a non-empty region",
                 call. = FALSE)
        } else if (i == 1) {
            p <- exp(.bfpwr_lpnorm_interval(lower = region[1,],
                                            upper = region[2,],
                                            mean = mean[1],
                                            sd = sqrt(sigma[1:1])))
        } else if (method == "lpmvnorm") {
            p <- exp(mvtnorm::lpmvnorm(lower = region[1, ],
                                       upper = region[2, ],
                                       mean = mean,
                                       chol = cholFactor,
                                       M = ngrid,
                                       w = w,
                                       ...))
        } else {
            p <- mvtnorm::pmvnorm(lower = region[1, ],
                                  upper = region[2, ],
                                  mean  = mean,
                                  sigma = sigma,
                                  seed = 42,
                                  keepAttr = FALSE,
                                  ...)
        }
        return(p)
    })

    if (any(is.na(regionprobs))) {
        stop("Sequential integration returned NA/NaN for a non-empty region",
             call. = FALSE)
    }

    sum(regionprobs)
}


#' @title Generate Integration Regions for Cumulative Z-Statistics Based on One
#'     Critical Value
#'
#' @description Constructs the per-stage integration regions for cumulative
#'     z-statistics given one lower and one upper critical value.
#'
#' @param zcrit0 Numeric vector of lower critical values (H0 boundaries).
#' @param zcrit1 Numeric vector of upper critical values (H1 boundaries).
#' @param direction Optional one-sided boundary direction, either
#'     \code{"positive"} or \code{"negative"}. If omitted, inferred from the
#'     supplied finite boundaries.
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

## NaN marks an empty boundary; plain NA means an invalid critical value.
.bfseq_has_plain_na <- function(x) {
    any(is.na(x) & !is.nan(x))
}

genregions1 <- function(zcrit0, zcrit1, direction = NULL) {
    stopifnot(all(is.numeric(zcrit0)),
              all(is.numeric(zcrit1)),
              length(zcrit0) == length(zcrit1))
    if (.bfseq_has_plain_na(zcrit0) || .bfseq_has_plain_na(zcrit1)) {
        stop("Critical values cannot contain NA; use NaN for empty boundaries.",
             call. = FALSE)
    }

    m <- length(zcrit1)
    if (is.null(direction)) {
        direction <- .bfseq_one_critical_direction(zcrit0 = zcrit0,
                                                   zcrit1 = zcrit1)
    }
    stopifnot(length(direction) == 1,
              direction %in% c("positive", "negative"))
    stages <- lapply(seq_len(m), function(i) {
        .bfseq_genregions1_stage(zcrit0 = zcrit0[seq_len(i)],
                                 zcrit1 = zcrit1[seq_len(i)],
                                 direction = direction)
    })

    return(list(H1 = lapply(stages, `[[`, "H1"),
                H0 = lapply(stages, `[[`, "H0")))
}

## Infer whether a one-critical-value design stops in the positive or negative
## direction from the finite H0/H1 boundaries.
.bfseq_one_critical_direction <- function(zcrit0, zcrit1) {
    H0nan <- is.nan(zcrit0)
    finiteH0 <- !H0nan
    finite <- !H0nan & !is.nan(zcrit1)
    finiteH1 <- !is.nan(zcrit1)
    if (any(finite) && all(zcrit1[finite] >= zcrit0[finite])) {
        return("positive")
    }
    if (any(finite) && all(zcrit1[finite] < zcrit0[finite])) {
        return("negative")
    }
    if (!any(finite) && any(finiteH1) && all(zcrit1[finiteH1] >= 0)) {
        return("positive")
    }
    if (!any(finite) && any(finiteH1) && all(zcrit1[finiteH1] <= 0)) {
        return("negative")
    }
    if (!any(finite) && !any(finiteH1) &&
        any(finiteH0) && all(zcrit0[finiteH0] >= 0)) {
        return("positive")
    }
    if (!any(finite) && !any(finiteH1) &&
        any(finiteH0) && all(zcrit0[finiteH0] <= 0)) {
        return("negative")
    }
    stop("Inconsistent critical values: direction cannot be inferred.")
}

## Build the H1/H0 integration regions for the final stage of a one-sided or
## one-critical-value schedule.
.bfseq_genregions1_stage <- function(zcrit0, zcrit1, direction = NULL) {
    stopifnot(all(is.numeric(zcrit0)),
              all(is.numeric(zcrit1)),
              length(zcrit0) == length(zcrit1))

    if (is.null(direction)) {
        direction <- .bfseq_one_critical_direction(zcrit0 = zcrit0,
                                                   zcrit1 = zcrit1)
    }
    stopifnot(length(direction) == 1,
              direction %in% c("positive", "negative"))
    H0nan <- is.nan(zcrit0)
    H1nan <- is.nan(zcrit1)
    i <- length(zcrit1)
    matH1 <- matrix(nrow = 2, ncol = i)
    for (j in seq_len(i)) {
        if (i == j) {
            if (direction == "positive") {
                lower <- zcrit1[j]
                upper <- Inf
            } else {
                lower <- -Inf
                upper <- zcrit1[j]
            }
        } else {
            if (direction == "positive") {
                if (H0nan[j]) lower <- -Inf
                else lower <- zcrit0[j]
                if (H1nan[j]) upper <- Inf
                else upper <- zcrit1[j]
            } else {
                if (H1nan[j]) lower <- -Inf
                else lower <- zcrit1[j]
                if (H0nan[j]) upper <- Inf
                else upper <- zcrit0[j]
            }
        }
        matH1[, j] <- c(lower, upper)
    }

    matH0 <- matH1
    if (direction == "positive") {
        matH0[, i] <- c(-Inf, zcrit0[i])
    } else {
        matH0[, i] <- c(zcrit0[i], Inf)
    }

    list(H1 = list(matH1), H0 = list(matH0))
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
#'     (extending from -Inf to this upper bound) and the second row the lower
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
    if (.bfseq_has_plain_na(zcrit0) || .bfseq_has_plain_na(zcrit1)) {
        stop("Critical values cannot contain NA; use NaN for empty boundaries.",
             call. = FALSE)
    }

    ## if (strict == TRUE & ncol(zcrit0) > 10) {
    ##     warning("strict = TRUE with many stages may cause numerical problems")
    ## }

    m <- ncol(zcrit0)
    stages <- lapply(seq_len(m), function(i) {
        .bfseq_genregions2_stage(zcrit0 = zcrit0[, seq_len(i), drop = FALSE],
                                 zcrit1 = zcrit1[, seq_len(i), drop = FALSE],
                                 strict = strict)
    })

    list(H1 = lapply(stages, `[[`, "H1"),
         H0 = lapply(stages, `[[`, "H0"))
}

## Build the final-stage integration regions when both tails can stop for H1.
## strict = TRUE keeps all continuation paths; FALSE keeps the dominant paths.
.bfseq_genregions2_stage <- function(zcrit0, zcrit1, strict = FALSE) {
    stopifnot(
        is.matrix(zcrit0),
        is.matrix(zcrit1),
        all(dim(zcrit0) == dim(zcrit1)),
        nrow(zcrit0) == 2
    )

    i <- ncol(zcrit0)
    H0nan <- apply(zcrit0, 2, function(x) any(is.nan(x)))

    H1i <- lapply(seq_len(i), function(j) {
        if (i == j) {
            list(
                c(-Inf, zcrit1[1, j]),
                c(zcrit1[2, j], Inf)
            )
        } else {
            if (H0nan[j] == TRUE) {
                list(c(zcrit1[1, j], zcrit1[2, j]))
            } else {
                list(
                    c(zcrit1[1, j], zcrit0[1, j]),
                    c(zcrit0[2, j], zcrit1[2, j])
                )
            }
        }
    })

    H0i <- H1i
    H0i[[i]] <- list(c(zcrit0[1, i], zcrit0[2, i]))

    if (strict == TRUE) {
        combosH1i <- expand.grid(lapply(H1i, seq_along))
        combosH0i <- expand.grid(lapply(H0i, seq_along))
    } else {
        combosH1i <- rbind(rep(1, length(H1i)), sapply(H1i, length))
        if (i <= 1 + sum(H0nan)) {
            combosH0i <- matrix(rep(1, length(H0i)), nrow = 1)
        } else {
            combosH0i <- rbind(rep(1, length(H0i)), sapply(H0i, length))
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

    list(H1 = makeregions(combosH1i, H1i),
         H0 = makeregions(combosH0i, H0i))
}

## Count the integration regions induced by strict two-sided continuation paths.
.count_strict_two_sided_regions <- function(zcrit0) {
    stopifnot(
        is.matrix(zcrit0),
        nrow(zcrit0) == 2
    )

    m <- ncol(zcrit0)
    H0nan <- apply(zcrit0, 2, function(x) any(is.nan(x)))
    H1 <- H0 <- numeric(m)
    finiteH0 <- 0L

    for (i in seq_len(m)) {
        ## Previous finite H0 boundaries split the continuation region into
        ## lower and upper paths; strict = TRUE integrates all combinations.
        npaths <- 2^finiteH0
        H1[i] <- 2*npaths
        H0[i] <- if (H0nan[i]) 0 else npaths

        if (!H0nan[i]) {
            finiteH0 <- finiteH0 + 1L
        }
    }

    list(
        total = sum(H1 + H0),
        perStage = H1 + H0,
        H0nan = H0nan,
        firstH0 = match(FALSE, H0nan)
    )
}

#' @title Compute Critical Z-Values for Bayes Factors
#'
#' @description Computes critical z-values for Bayes factors using normal,
#'     directional normal, or normal moment priors.
#'
#' @param k Positive numeric. Bayes factor threshold (BF01 oriented in favor of
#'     H0)
#' @param se Positive numeric. Standard error
#' @param null Parameter value under the point null hypothesis. Defaults to
#'     \code{0}
#' @param mu Numeric. Prior location. Not taken into account for \code{type =
#'     "moment"}
#' @param tau Non-negative numeric. Prior scale
#' @param type Character. One of "normal" (point null vs. normal alternative),
#'     "directional" (directional null \eqn{\theta \leq 0}{theta <= 0} vs.
#'     directional alternative \eqn{\theta > 0}{theta > 0} with marginal normal
#'     prior), or "moment" (point null vs. normal moment alternative)
#'
#' @return Numeric vector of critical z-value(s)
#'
#' @examples
#' se <- 0.05
#' null <- 0
#' pm <- 0.2
#' psd <- 0.1
#' zcrit1 <- zcrit(k = 3, se = se, null = null, mu = pm, tau = psd,
#'                 type = "normal")
#' bf01(estimate = null + zcrit1*se, se = se, null = null, pm = pm, psd = psd)
#'
#' ## tau = 0 leads to point alternative
#' zcrit2 <- zcrit(k = 5, se = se, null = null, mu = pm, tau = 0,
#'                 type = "normal")
#' bf01(estimate = null + zcrit2*se, se = se, null = null, pm = pm, psd = 0)
#'
#' zcrit3 <- zcrit(k = 5, se = se, null = null, mu = pm, tau = psd,
#'                 type = "directional")
#' dirbf01(estimate = null + zcrit3*se, se = se, null = null, pm = pm,
#'         psd = psd)
#'
#' zcrit4 <- zcrit(k = 1/10, se = se, null = null, tau = psd,
#'                 type = "moment")
#' nmbf01(estimate = null + zcrit4*se, se = se, null = null, psd = psd)
#'
#' @noRd
#'
#' @keywords internal
zcrit <- function(k, se, null = 0, mu = NULL, tau,
                  type = c("normal", "directional", "moment")) {

    type <- match.arg(type)
    if (type != "moment") {
        mu <- mu - null
    }

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
