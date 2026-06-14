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



.bfpwr_root_value <- function(f, x) {
    ans <- try(suppressWarnings(f(x)), silent = TRUE)
    if (inherits(ans, "try-error") || length(ans) != 1 ||
        !is.numeric(ans) || !is.finite(ans)) {
        return(NaN)
    }
    ans
}

.bfpwr_certified_root <- function(f, x0, x1, f0 = NaN, f1 = NaN,
                                  final_fun = f, tolerance = 1e-5, ...) {
    if (x0 == x1) {
        return(structure("degenerate root interval", class = "try-error"))
    }

    if (!is.finite(f0)) {
        f0 <- .bfpwr_root_value(f = f, x = x0)
    }
    if (!is.finite(f1)) {
        f1 <- .bfpwr_root_value(f = f, x = x1)
    }
    if (!is.finite(f0) || !is.finite(f1) || f0*f1 > 0) {
        return(structure("root not bracketed", class = "try-error"))
    }
    if (f0 == 0) {
        if (identical(final_fun, f)) {
            return(x0)
        }
        final0 <- .bfpwr_root_value(f = final_fun, x = x0)
        if (is.finite(final0) && abs(final0) <= tolerance) {
            return(x0)
        }
    }
    if (f1 == 0) {
        if (identical(final_fun, f)) {
            return(x1)
        }
        final1 <- .bfpwr_root_value(f = final_fun, x = x1)
        if (is.finite(final1) && abs(final1) <= tolerance) {
            return(x1)
        }
    }

    root <- try(stats::uniroot(f = f, interval = sort(c(x0, x1)),
                               extendInt = "no", ...)$root,
                silent = TRUE)
    if (inherits(root, "try-error")) {
        return(root)
    }
    if (identical(final_fun, f)) {
        return(root)
    }

    residual <- .bfpwr_root_value(f = final_fun, x = root)
    if (is.finite(residual) && abs(residual) <= tolerance) {
        return(root)
    }

    final_f0 <- .bfpwr_root_value(f = final_fun, x = x0)
    final_f1 <- .bfpwr_root_value(f = final_fun, x = x1)
    if (!is.finite(final_f0) || !is.finite(final_f1) ||
        final_f0*final_f1 > 0) {
        return(structure("root not certified by final function",
                         class = "try-error"))
    }
    if (final_f0 == 0) return(x0)
    if (final_f1 == 0) return(x1)

    try(stats::uniroot(f = final_fun, interval = sort(c(x0, x1)),
                       extendInt = "no", ...)$root,
        silent = TRUE)
}

.bfpwr_integrate_dots <- function(dots, rel.tol.default = NULL) {
    if (length(dots) == 0) {
        out <- list()
    } else {
        dot_names <- names(dots)
        if (is.null(dot_names)) {
            out <- list()
        } else {
            integrate_names <- c("subdivisions", "rel.tol", "abs.tol",
                                 "stop.on.error", "keep.xy")
            keep <- nzchar(dot_names) & dot_names %in% integrate_names
            out <- dots[keep]
        }
    }
    if (!is.null(rel.tol.default) && !("rel.tol" %in% names(out))) {
        out$rel.tol <- rel.tol.default
    }
    out
}

.bfpwr_uniroot_dots <- function(dots) {
    if (length(dots) == 0) {
        return(list())
    }
    dot_names <- names(dots)
    if (is.null(dot_names)) {
        return(list())
    }
    keep_names <- c("tol", "maxiter", "trace", "check.conv")
    keep <- nzchar(dot_names) & dot_names %in% keep_names
    dots[keep]
}

.bfpwr_tcrit_call <- function(args) {
    tcritWarnings <- character()
    value <- withCallingHandlers(
        do.call(tcrit, args),
        warning = function(w) {
            tcritWarnings <<- c(tcritWarnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    list(value = value, warnings = tcritWarnings)
}

.bfpwr_tcrit_has_max_bf_warning <- function(warnings) {
    any(grepl("maximum BF is less than k", warnings, fixed = TRUE))
}

.bfpwr_tcrit_result <- function(...) {
    args <- list(...)
    trange <- args$trange
    if (is.null(trange)) {
        trange <- "adaptive"
    }

    result <- .bfpwr_tcrit_call(args = args)
    status <- .bfpwr_tcrit_status(value = result$value,
                                  warnings = result$warnings,
                                  trange = trange)
    if (identical(status, "search_failed") &&
        identical(args$alternative, "two.sided") &&
        is.numeric(trange) &&
        .bfpwr_tcrit_has_max_bf_warning(result$warnings)) {
        adaptiveArgs <- args
        adaptiveArgs$trange <- "adaptive"
        adaptiveResult <- .bfpwr_tcrit_call(args = adaptiveArgs)
        adaptiveStatus <- .bfpwr_tcrit_status(
            value = adaptiveResult$value,
            warnings = adaptiveResult$warnings,
            trange = "adaptive"
        )
        if (identical(adaptiveStatus, "impossible")) {
            status <- "impossible"
        } else {
            result$warnings <- unique(c(
                result$warnings,
                paste0(
                    "numeric 'trange' does not contain the two-sided ",
                    "BF01 = k boundary"
                )
            ))
        }
    }
    list(
        value = result$value,
        status = status,
        warnings = result$warnings
    )
}

.bfpwr_tcrit_status <- function(value, warnings, trange = "adaptive") {
    if (.bfpwr_tcrit_has_max_bf_warning(warnings)) {
        if (is.numeric(trange)) {
            return("search_failed")
        }
        return("impossible")
    }
    if (any(grepl("Adaptive t critical-value search reached",
                  warnings, fixed = TRUE))) {
        return("tail_cutoff")
    }
    if (any(grepl("requires 'search_limit'", warnings, fixed = TRUE))) {
        return("error")
    }
    if (any(grepl("Numerical problems", warnings, fixed = TRUE))) {
        return("search_failed")
    }
    if (length(value) > 0 && all(is.finite(value))) {
        return("ok")
    }
    if (any(!is.finite(value))) {
        return("search_failed")
    }
    "search_failed"
}

.bfpwr_tcrit_unhandled_warnings <- function(results) {
    warnings <- unlist(lapply(results, `[[`, "warnings"), use.names = FALSE)
    if (length(warnings) == 0) {
        return(character())
    }
    known <- .bfpwr_tcrit_has_max_bf_warning(warnings) |
        grepl("Adaptive t critical-value search reached", warnings,
              fixed = TRUE) |
        grepl("requires 'search_limit'", warnings, fixed = TRUE) |
        grepl("Numerical problems", warnings, fixed = TRUE)
    unique(warnings[!known])
}

.bfseq_t_boundary_statuses <- function(results) {
    vapply(results, `[[`, character(1), "status")
}

.bfseq_t_boundary_status_message <- function(results, boundary,
                                             looks = seq_along(results)) {
    stopifnot(boundary %in% c("H0", "H1"))
    if (length(looks) != length(results)) {
        stop("internal error: boundary status look labels do not match results",
             call. = FALSE)
    }

    statuses <- .bfseq_t_boundary_statuses(results)
    allowed <- if (boundary == "H0") {
        c("ok", "impossible", "tail_cutoff")
    } else {
        c("ok", "tail_cutoff")
    }
    invalid <- !(statuses %in% allowed)
    if (!any(invalid)) {
        return(NULL)
    }

    look <- which(invalid)[1]
    warningText <- results[[look]]$warnings
    detail <- if (length(warningText) > 0) {
        paste(unique(warningText), collapse = "; ")
    } else {
        paste0("status: ", statuses[[look]])
    }
    paste0(
        "Failed to compute ", boundary,
        " sequential t stopping boundary at look ", looks[[look]], ": ",
        detail,
        ". Widen numeric 'trange' or use adaptive 'trange' with a smaller ",
        "'tail.eps'."
    )
}

.bfseq_validate_t_boundary_statuses <- function(results, boundary,
                                               looks = seq_along(results)) {
    msg <- .bfseq_t_boundary_status_message(results = results,
                                            boundary = boundary,
                                            looks = looks)
    if (!is.null(msg)) {
        stop(msg, call. = FALSE)
    }
    invisible(.bfseq_t_boundary_statuses(results))
}

.bfseq_warn_t_boundary_statuses <- function(results0, results1, tail.eps) {
    statuses0 <- .bfseq_t_boundary_statuses(results0)
    statuses1 <- .bfseq_t_boundary_statuses(results1)

    unhandled <- .bfpwr_tcrit_unhandled_warnings(c(results0, results1))
    for (msg in unhandled) {
        warning(msg, call. = FALSE)
    }

    impossible <- sum(statuses0 == "impossible")
    if (impossible > 0) {
        warning(paste0(
            "No H0 sequential t stopping boundary exists in ",
            impossible,
            " boundary search(es); the corresponding H0 stopping regions ",
            "are treated as empty."
        ), call. = FALSE)
    }

    tailCutoff <- sum(c(statuses0, statuses1) == "tail_cutoff")
    if (tailCutoff > 0) {
        warning(paste0(
            "Adaptive t critical-value search reached the predictive tail ",
            "cutoff in ",
            tailCutoff,
            " sequential boundary search(es); each unresolved boundary has ",
            "marginal tail probability <= ", format(tail.eps),
            ". Pass a wider numeric 'trange' interval to search exact bounds."
        ), call. = FALSE)
    }
}

.bfpwr_one_sided_direction <- function(alternative, f_origin) {
    stopifnot(
        alternative %in% c("greater", "less"),
        length(f_origin) == 1,
        is.numeric(f_origin),
        is.finite(f_origin),
        f_origin != 0
    )

    if (alternative == "greater") {
        if (f_origin > 0) 1 else -1
    } else {
        if (f_origin > 0) -1 else 1
    }
}

.bfpwr_one_sided_tail_limit <- function(direction, origin, step_scale, mean,
                                        sd, tail.eps) {
    stopifnot(
        length(direction) == 1,
        direction %in% c(-1, 1),
        length(origin) == 1,
        is.numeric(origin),
        is.finite(origin),
        length(step_scale) == 1,
        is.numeric(step_scale),
        is.finite(step_scale),
        step_scale > 0,
        length(mean) == 1,
        is.numeric(mean),
        is.finite(mean),
        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        sd > 0,
        length(tail.eps) == 1,
        is.numeric(tail.eps),
        is.finite(tail.eps),
        tail.eps > 0,
        tail.eps < 0.5
    )

    if (direction > 0) {
        limit <- stats::qnorm(p = tail.eps, mean = mean, sd = sd,
                              lower.tail = FALSE)
        if (!is.finite(limit) || limit <= origin) {
            limit <- origin
        }
        tail_probability <- stats::pnorm(q = limit, mean = mean, sd = sd,
                                         lower.tail = FALSE)
    } else {
        limit <- stats::qnorm(p = tail.eps, mean = mean, sd = sd,
                              lower.tail = TRUE)
        if (!is.finite(limit) || limit >= origin) {
            limit <- origin
        }
        tail_probability <- stats::pnorm(q = limit, mean = mean, sd = sd,
                                         lower.tail = TRUE)
    }

    list(
        direction = direction,
        limit = limit,
        search_limit = abs(limit - origin)/step_scale,
        tail_probability = tail_probability,
        tail.eps = tail.eps
    )
}

.bfpwr_one_sided_tail_limits <- function(origin, step_scale, mean, sd,
                                         tail.eps) {
    list(
        positive = .bfpwr_one_sided_tail_limit(
            direction = 1, origin = origin, step_scale = step_scale,
            mean = mean, sd = sd, tail.eps = tail.eps
        ),
        negative = .bfpwr_one_sided_tail_limit(
            direction = -1, origin = origin, step_scale = step_scale,
            mean = mean, sd = sd, tail.eps = tail.eps
        ),
        tail.eps = tail.eps
    )
}

.bfpwr_select_one_sided_search_limit <- function(search_limit, direction,
                                                 origin, step_scale) {
    stopifnot(
        length(direction) == 1,
        direction %in% c(-1, 1),
        length(origin) == 1,
        is.numeric(origin),
        is.finite(origin),
        length(step_scale) == 1,
        is.numeric(step_scale),
        is.finite(step_scale),
        step_scale > 0
    )

    if (is.null(search_limit)) {
        return(NULL)
    }

    if (is.numeric(search_limit) && length(search_limit) == 1) {
        stopifnot(is.finite(search_limit), search_limit >= 0)
        return(list(
            direction = direction,
            limit = origin + direction*step_scale*search_limit,
            search_limit = search_limit,
            tail_probability = NA_real_,
            tail.eps = NA_real_
        ))
    }

    stopifnot(is.list(search_limit))
    selected <- if (direction > 0) search_limit$positive else search_limit$negative
    stopifnot(
        is.list(selected),
        length(selected$search_limit) == 1,
        is.numeric(selected$search_limit),
        is.finite(selected$search_limit),
        selected$search_limit >= 0
    )
    selected
}

.bfpwr_one_sided_adaptive_result <- function(root, search_limit_reached,
                                             selected_limit = NULL) {
    if (is.null(selected_limit)) {
        selected_limit <- list(
            direction = NA_real_,
            limit = NA_real_,
            search_limit = NA_real_,
            tail_probability = NA_real_,
            tail.eps = NA_real_
        )
    }
    list(
        root = root,
        search_limit_reached = search_limit_reached,
        direction = selected_limit$direction,
        limit = selected_limit$limit,
        search_limit = selected_limit$search_limit,
        tail_probability = selected_limit$tail_probability,
        tail.eps = selected_limit$tail.eps
    )
}

.bfpwr_residual_certified_root <- function(scout_fun, certify_fun, x0, x1,
                                           final_fun = certify_fun,
                                           tolerance = 1e-5, ...) {
    root <- try(stats::uniroot(f = scout_fun, interval = sort(c(x0, x1)),
                               extendInt = "no", ...)$root,
                silent = TRUE)
    if (inherits(root, "try-error") || !is.numeric(root) ||
        length(root) != 1 || !is.finite(root)) {
        return(structure("scout root failed", class = "try-error"))
    }

    residual <- .bfpwr_root_value(f = certify_fun, x = root)
    if (!is.finite(residual) || abs(residual) > tolerance) {
        return(structure("scout root not certified", class = "try-error"))
    }

    final_residual <- .bfpwr_root_value(f = final_fun, x = root)
    if (is.finite(final_residual) && abs(final_residual) <= tolerance) {
        return(root)
    }

    structure("scout root not certified by final function",
              class = "try-error")
}

## One-sided adaptive boundary search. The fast scout function is used only to
## locate candidate brackets. search_fun performs stable bracket checks, while
## returned roots and tail cutoffs must be certified by certify_fun.
## The return value is a list with root and search_limit_reached.
.bfpwr_one_sided_adaptive_root <- function(certify_fun, scout_fun, alternative,
                                           search_fun = certify_fun,
                                           origin = 0, step_scale = 1,
                                           try_opposite = FALSE,
                                           search_limit = NULL,
                                           steps = c(0.1, 0.25, 0.5, 1, 1.5,
                                                     2, 2.5, 3, 3.5),
                                           tail_steps = c(4, 8, 16, 32, 64,
                                                          128, 256),
                                           scout_tail_steps = c(4, 5, 6, 7, 8,
                                                                16, 32, 64,
                                                                128, 256),
                                           scout_tolerance = 1e-5,
                                           ...) {
    f_origin <- .bfpwr_root_value(f = certify_fun, x = origin)
    if (!is.finite(f_origin)) {
        return(.bfpwr_one_sided_adaptive_result(
            root = structure("non-finite root start", class = "try-error"),
            search_limit_reached = FALSE
        ))
    }
    if (f_origin == 0) {
        return(.bfpwr_one_sided_adaptive_result(
            root = origin, search_limit_reached = FALSE
        ))
    }
    f_origin_search <- .bfpwr_root_value(f = search_fun, x = origin)
    if (!is.finite(f_origin_search)) {
        return(.bfpwr_one_sided_adaptive_result(
            root = structure("non-finite search start", class = "try-error"),
            search_limit_reached = FALSE
        ))
    }

    expected_direction <- .bfpwr_one_sided_direction(alternative = alternative,
                                                     f_origin = f_origin)
    if (is.null(search_limit)) {
        return(.bfpwr_one_sided_adaptive_result(
            root = structure("missing finite search limit", class = "try-error"),
            search_limit_reached = FALSE
        ))
    }

    directions <- if (try_opposite) {
        c(expected_direction, -expected_direction)
    } else {
        expected_direction
    }
    search_limit_reached <- FALSE
    selected_limit <- NULL
    for (direction in directions) {
        selected_limit <- .bfpwr_select_one_sided_search_limit(
            search_limit = search_limit, direction = direction,
            origin = origin, step_scale = step_scale
        )
        if (is.null(selected_limit)) {
            return(.bfpwr_one_sided_adaptive_result(
                root = structure("missing finite search limit",
                                 class = "try-error"),
                search_limit_reached = FALSE
            ))
        }
        current_search_limit <- selected_limit$search_limit
        if (current_search_limit <= 0) {
            search_limit_reached <- TRUE
            next
        }

        scan_steps <- sort(unique(c(steps, scout_tail_steps)))
        scan_steps <- scan_steps[is.finite(scan_steps) & scan_steps > 0 &
                                 scan_steps <= current_search_limit]
        if (!current_search_limit %in% scan_steps) {
            scan_steps <- sort(c(scan_steps, current_search_limit))
        }
        exact_tail_steps <- sort(unique(c(tail_steps, current_search_limit)))
        exact_tail_steps <- exact_tail_steps[
            is.finite(exact_tail_steps) & exact_tail_steps > 0 &
                exact_tail_steps <= current_search_limit
        ]

        last_finite_step <- 0
        last_finite_x <- origin
        finite_x <- numeric(0)
        scout_prev_x <- origin
        scout_prev_f <- f_origin_search

        ## First use the fast direct integral only as a scout. A finite scout
        ## sign change is never returned unless the stable BF path certifies it.
        for (step in scan_steps) {
            x1 <- origin + direction*step_scale*step
            f1 <- .bfpwr_root_value(f = scout_fun, x = x1)
            if (!is.finite(f1)) {
                break
            }
            if (is.finite(scout_prev_f) && scout_prev_f*f1 <= 0) {
                root <- .bfpwr_residual_certified_root(
                    scout_fun = scout_fun, certify_fun = search_fun,
                    final_fun = certify_fun,
                    x0 = scout_prev_x, x1 = x1,
                    tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
                root <- .bfpwr_certified_root(
                    f = search_fun, x0 = scout_prev_x, x1 = x1,
                    final_fun = certify_fun, tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
            }
            scout_prev_x <- x1
            scout_prev_f <- f1
            last_finite_step <- step
            last_finite_x <- x1
            finite_x <- c(finite_x, x1)
        }

        fprev <- f_origin_search
        xprev <- origin
        certified_limit_finite <- FALSE
        if (last_finite_step > 0) {
            f_last <- .bfpwr_root_value(f = search_fun, x = last_finite_x)
            if (is.finite(f_last) && f_origin_search*f_last <= 0) {
                x_bracket0 <- origin
                f_bracket0 <- f_origin_search
                x_bracket1 <- last_finite_x
                f_bracket1 <- f_last
                if (length(finite_x) > 1) {
                    for (x_candidate in rev(finite_x[-length(finite_x)])) {
                        f_candidate <- .bfpwr_root_value(f = search_fun,
                                                         x = x_candidate)
                        if (!is.finite(f_candidate)) {
                            next
                        }
                        if (f_candidate*f_bracket1 <= 0) {
                            x_bracket0 <- x_candidate
                            f_bracket0 <- f_candidate
                            break
                        }
                        x_bracket1 <- x_candidate
                        f_bracket1 <- f_candidate
                    }
                }
                root <- .bfpwr_certified_root(
                    f = search_fun, x0 = x_bracket0, x1 = x_bracket1,
                    f0 = f_bracket0, f1 = f_bracket1,
                    final_fun = certify_fun, tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
            }
            if (is.finite(f_last)) {
                fprev <- f_last
                xprev <- last_finite_x
                if (last_finite_step >= current_search_limit) {
                    f_limit_final <- .bfpwr_root_value(f = certify_fun,
                                                       x = last_finite_x)
                    if (is.finite(f_limit_final) &&
                        f_origin*f_limit_final <= 0) {
                        root <- .bfpwr_certified_root(
                            f = certify_fun, x0 = origin, x1 = last_finite_x,
                            f0 = f_origin, f1 = f_limit_final, ...
                        )
                        if (!inherits(root, "try-error")) {
                            return(.bfpwr_one_sided_adaptive_result(
                                root = root, search_limit_reached = FALSE,
                                selected_limit = selected_limit
                            ))
                        }
                    }
                    certified_limit_finite <- is.finite(f_limit_final) &&
                        f_origin*f_limit_final > 0
                }
            }
        }
        if (last_finite_step >= current_search_limit) {
            search_limit_reached <- search_limit_reached ||
                certified_limit_finite
            next
        }

        direction_limit_reached <- FALSE
        exact_steps <- exact_tail_steps[exact_tail_steps > last_finite_step]
        f_limit <- NaN
        limit_checked <- FALSE
        limit_ruled_out <- FALSE
        for (step in exact_steps) {
            x1 <- origin + direction*step_scale*step
            f1 <- if (step >= current_search_limit && is.finite(f_limit)) {
                f_limit
            } else {
                .bfpwr_root_value(f = search_fun, x = x1)
            }
            if (is.finite(f1) && fprev*f1 <= 0) {
                root <- .bfpwr_certified_root(
                    f = search_fun, x0 = xprev, x1 = x1, f0 = fprev,
                    f1 = f1, final_fun = certify_fun,
                    tolerance = scout_tolerance, ...
                )
                if (!inherits(root, "try-error")) {
                    return(.bfpwr_one_sided_adaptive_result(
                        root = root, search_limit_reached = FALSE,
                        selected_limit = selected_limit
                    ))
                }
            }
            if (is.finite(f1)) {
                xprev <- x1
                fprev <- f1
                if (step >= current_search_limit) {
                    f_limit_final <- .bfpwr_root_value(f = certify_fun,
                                                       x = x1)
                    if (is.finite(f_limit_final) &&
                        f_origin*f_limit_final <= 0) {
                        root <- .bfpwr_certified_root(
                            f = certify_fun, x0 = origin, x1 = x1,
                            f0 = f_origin, f1 = f_limit_final, ...
                        )
                        if (!inherits(root, "try-error")) {
                            return(.bfpwr_one_sided_adaptive_result(
                                root = root, search_limit_reached = FALSE,
                                selected_limit = selected_limit
                            ))
                        }
                    }
                    direction_limit_reached <- is.finite(f_limit_final) &&
                        f_origin*f_limit_final > 0
                }
            }
            ## If the exact path is still far from the threshold, check the
            ## finite search limit before walking every tail step. Near-root
            ## cases continue locally instead of jumping to the search limit.
            if (!limit_checked && step < current_search_limit &&
                is.finite(fprev) && abs(fprev) > 0.1) {
                x_limit <- origin + direction*step_scale*current_search_limit
                f_limit <- .bfpwr_root_value(f = search_fun, x = x_limit)
                limit_checked <- TRUE
                if (is.finite(f_limit)) {
                    f_limit_final <- .bfpwr_root_value(f = certify_fun,
                                                       x = x_limit)
                    if (is.finite(f_limit_final) &&
                        f_origin*f_limit_final <= 0) {
                        root <- .bfpwr_certified_root(
                            f = certify_fun, x0 = origin, x1 = x_limit,
                            f0 = f_origin, f1 = f_limit_final, ...
                        )
                        if (!inherits(root, "try-error")) {
                            return(.bfpwr_one_sided_adaptive_result(
                                root = root, search_limit_reached = FALSE,
                                selected_limit = selected_limit
                            ))
                        }
                    }
                    direction_limit_reached <- is.finite(f_limit_final) &&
                        f_origin*f_limit_final > 0
                    if (direction_limit_reached && fprev*f_limit > 0) {
                        limit_ruled_out <- TRUE
                        break
                    }
                }
            }
        }
        if (limit_ruled_out) {
            search_limit_reached <- TRUE
            next
        }

        search_limit_reached <- search_limit_reached || direction_limit_reached
    }

    .bfpwr_one_sided_adaptive_result(
        root = structure("root not bracketed", class = "try-error"),
        search_limit_reached = search_limit_reached,
        selected_limit = selected_limit
    )
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
#' @param trange Numerical search strategy. Can be either \code{"adaptive"}
#'     (default) or an interval. One-sided adaptive searches require a finite
#'     \code{search_limit} supplied by the caller.
#' @param search_limit Finite one-sided adaptive search limit, either as a
#'     scalar in \eqn{t}-statistic units or as the directional object returned
#'     by \code{.bfpwr_one_sided_tail_limits()}.
#' @param ... Optional numerical controls. For numeric ranges and two-sided
#'     adaptive searches, arguments are passed to \code{stats::uniroot}. In
#'     adaptive one-sided searches, \code{subdivisions}, \code{rel.tol},
#'     \code{abs.tol}, \code{stop.on.error}, and \code{keep.xy} are used for BF
#'     integration, while \code{tol}, \code{maxiter}, \code{trace}, and
#'     \code{check.conv} are passed to \code{stats::uniroot}.
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
#' search_limit <- .bfpwr_one_sided_tail_limits(origin = 0, step_scale = 1,
#'                                              mean = 0, sd = 2,
#'                                              tail.eps = 1e-3)
#' tcrit2 <- tcrit(k = k, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
#'                 pdf = pdf, alternative = alternative, type = type,
#'                 search_limit = search_limit)
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
                  trange = "adaptive", search_limit = NULL, ...) {

    ## determine t-statistic for which BF = k
    dots <- list(...)
    searchDots <- .bfpwr_integrate_dots(dots = dots,
                                        rel.tol.default = 1e-2)
    rootDots <- .bfpwr_uniroot_dots(dots = dots)
    rootFun <- function(t) {
        tbf01(t = t, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
              pdf = pdf, type = type, alternative = alternative,
              log = TRUE) - log(k)
    }
    rootFunSearch <- function(t) {
        do.call(tbf01, c(list(
            t = t, n1 = n1, n2 = n2, plocation = plocation, pscale = pscale,
            pdf = pdf, type = type, alternative = alternative, log = TRUE
        ), searchDots)) - log(k)
    }
    if (type == "two.sample") {
        pars <- .tbf01_pars(n1 = n1, n2 = n2, type = type)
    } else {
        pars <- .tbf01_pars(n1 = n1, n2 = n1, type = type)
    }
    region <- .tbf01_prior_region(plocation = plocation, pscale = pscale,
                                  pdf = pdf, alternative = alternative)
    rootFunFast <- function(t) {
        do.call(.tbf01_log_fast, c(list(
            t = t, df = pars$df, neff = pars$neff,
            plocation = plocation, pscale = pscale, pdf = pdf,
            region = region
        ), searchDots)) - log(k)
    }

    if (alternative == "two.sided") {
        ## guess search range based on search range from z-test BF
        if (!is.numeric(trange) && trange == "adaptive") {
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
            meant <- mean(trange)
            searchIntLow <- c(trange[1], meant)
            searchIntUp <- c(meant, trange[2])
        }
        if (k > 1) {
            ## Check impossible H0 boundaries only in the interval searched below.
            ## This avoids unconstrained wrong-tail evaluations in tbf01().
            maxInt <- c(searchIntLow[1], searchIntUp[2])
            opt <- try(stats::optimize(f = function(t) {
                                           ans <- suppressWarnings(rootFun(t))
                                           if (is.finite(ans)) ans else -Inf
                                       },
                                       interval = maxInt,
                                       maximum = TRUE),
                       silent = TRUE)
            if (!inherits(opt, "try-error") &&
                is.finite(opt$objective) && opt$objective < 0) {
                warning("maximum BF is less than k; BF01 = k impossible")
                return(c(NaN, NaN))
            }
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
        search_limit_reached <- FALSE
        search_limit_missing <- FALSE
        if (!is.numeric(trange) && trange == "adaptive") {
            if (is.null(search_limit)) {
                search_limit_missing <- TRUE
                warning(paste0(
                    "Adaptive one-sided t critical-value search requires ",
                    "'search_limit'"
                ))
                res <- structure("missing finite search limit",
                                 class = "try-error")
            } else {
                search <- do.call(.bfpwr_one_sided_adaptive_root, c(list(
                    certify_fun = rootFun, search_fun = rootFunSearch,
                    scout_fun = rootFunFast, alternative = alternative,
                    origin = 0, step_scale = 1, try_opposite = FALSE,
                    search_limit = search_limit
                ), rootDots))
                res <- search$root
                search_limit_reached <- search$search_limit_reached
            }
        } else {
            searchint <- trange
            extend <- "no"

            suppressWarnings({
                res <- try(stats::uniroot(f = rootFun, interval = searchint,
                                          extendInt = extend, ...)$root,
                           silent = TRUE)
            })
        }
        if (inherits(res, "try-error")) {
            if (search_limit_missing) {
                tcrit <- NaN
            } else if (search_limit_reached) {
                warning(paste0(
                    "Adaptive t critical-value search reached the predictive ",
                    "tail cutoff without bracketing BF01 = k; pass a wider ",
                    "numeric 'trange' interval to search exact bounds."
                ))
            } else {
                warning("Numerical problems finding critical value")
            }
            tcrit <- NaN
        } else {
            tcrit <- res
        }
    }
    return(tcrit)
}
