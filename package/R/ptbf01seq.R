#' @title Sequential T-Test Bayes Factor Design
#'
#' @description Computes cumulative probabilities of observing \eqn{t}-test
#'     Bayes factors that provide evidence for the null hypothesis
#'     \eqn{H_0}{H0}, the alternative hypothesis \eqn{H_1}{H1}, or remain
#'     inconclusive in a sequential design. Also computes the expected sample
#'     size.
#'
#' @inheritParams ptbf01
#' @inheritParams pbf01seq
#' @param trange Critical \eqn{t}-statistic search strategy for the sequential
#'     stopping boundaries. Can be either \code{"adaptive"} (default) or a
#'     numeric interval. For one-sided adaptive searches, \code{tail.eps}
#'     determines the predictive tail cutoff used to stop finite scans. Use a
#'     numeric interval to force exact boundary searches over that interval.
#' @param tail.eps One-sided adaptive sequential boundary searches stop once the
#'     marginal predictive tail probability beyond the searched \eqn{t}
#'     statistic is at most \code{tail.eps}. The same cutoff is used for each
#'     adaptive boundary search; it is not split across looks or boundaries.
#'     If no boundary is found before the cutoff, that boundary is treated as
#'     unresolved with omitted marginal mass bounded by \code{tail.eps};
#'     smaller values search farther at additional computational cost. Defaults
#'     to \code{1e-3}.
#' @param ... Additional arguments passed to \code{mvtnorm::lpmvnorm}
#'
#' @inherit pbf01seq return
#'
#' @details The function constructs per-stage integration regions for cumulative
#'     z-statistics based on the Bayes factor thresholds \code{k1} and
#'     \code{k0}, then computes the probability of these regions under a
#'     predictive distribution defined by the asymptotic variance of the
#'     \eqn{t}-statistic and the normal design prior with \code{dpm} and
#'     \code{dpsd}. Integration is performed via \code{mvtnorm::lpmvnorm}.
#'
#' @examples
#' ## similar to example from Schönbrodt and Wagenmakers (2018, p. 138)
#' k0 <- 6
#' k1 <- 1/30
#' dpm <- 0.5
#' dpsd <- 0.1
#' plocation <- 0
#' pscale <- 1/sqrt(2)
#' pdf <- 1
#' type <- "two.sample"
#' alternative <- "greater"
#' n <- seq(40, 100, 10) # sample size (per group) per stage
#' res <- ptbf01seq(k1 = k1, k0 = k0, n = n, plocation = plocation,
#'                  pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd,
#'                  alternative = alternative, type = type)
#' res
#' plot(res) # show stopping probabilities
#' plot(res, zplot = TRUE) # show critical z-values
#'
#' @author Samuel Pawel
#'
#' @export
ptbf01seq <- function(k1, k0 = 1/k1, n, n1 = n, n2 = n, plocation = 0,
                      pscale = 1/sqrt(2), pdf = 1, dpm = plocation,
                      dpsd = pscale,
                      type = c("two.sample", "one.sample", "paired"),
                      alternative = c("two.sided", "less", "greater"),
                      strict = TRUE, trange = "adaptive",
                      tail.eps = 1e-3,
                      tail.nquad = 128, ...) {

    ## input checks
    stopifnot(
        length(k1) == 1,
        is.numeric(k1),
        is.finite(k1),
        k1 <= 1,

        length(k0) == 1,
        is.numeric(k0),
        is.finite(k0),
        k0 >= 1,

        length(n1) >= 1,
        is.numeric(n1),
        all(is.finite(n1)),
        all(1 < n1),

        length(n2) >= 1,
        length(n1) == length(n2),
        is.numeric(n2),
        all(is.finite(n2)),
        all(1 < n2),

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

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        0 <= dpsd,

        length(tail.eps) == 1,
        is.numeric(tail.eps),
        is.finite(tail.eps),
        tail.eps > 0,
        tail.eps < 0.5,

        .tbf01_valid_tail_nquad(tail.nquad),

        (is.numeric(trange) && length(trange) == 2 && all(is.finite(trange)) &&
         trange[2] > trange[1]) || (is.character(trange) && length(trange) == 1 &&
                                    !is.na(trange) && trange == "adaptive")
    )
    type <- match.arg(type)
    alternative <- match.arg(alternative)
    if (type != "two.sample") {
        if (all(n1 != n2)) {
            warning(paste0('different n1 and n2 supplied but type set to "', type,
                           '", using n = n1'))
            n2 <- n1
        }
    }

    ## effective sample size
    if (type == "two.sample") {
        neff <- 1/(1/n1 + 1/n2)
    } else {
        neff <- n1
    }

    ## get marginal mean and covariance matrix
    se <- 1/sqrt(neff) # standard error of SMD assuming variance is known
    pars <- predpars(se = se, dpm = dpm, dpsd = dpsd)
    mean <- pars$mean
    sigma <- pars$sigma

    ## get integration regions
    adaptiveOneSided <- alternative != "two.sided" &&
        !is.numeric(trange) && trange == "adaptive"
    tSearchLimits <- if (adaptiveOneSided) {
        lapply(seq_along(n1), function(i) {
            .bfpwr_one_sided_tail_limits(
                origin = 0, step_scale = 1, mean = mean[i],
                sd = sqrt(sigma[i, i]), tail.eps = tail.eps
            )
        })
    } else {
        vector("list", length(n1))
    }
    zk0Results <- lapply(seq_along(n1), function(i) {
        .bfpwr_tcrit_result(
            k = k0, n1 = n1[i], n2 = n2[i], plocation = plocation,
            pscale = pscale, pdf = pdf, alternative = alternative,
            type = type, trange = trange, search_limit = tSearchLimits[[i]],
            tail.nquad = tail.nquad
        )
    })
    zk1Results <- lapply(seq_along(n1), function(i) {
        .bfpwr_tcrit_result(
            k = k1, n1 = n1[i], n2 = n2[i], plocation = plocation,
            pscale = pscale, pdf = pdf, alternative = alternative,
            type = type, trange = trange, search_limit = tSearchLimits[[i]],
            tail.nquad = tail.nquad
        )
    })
    .bfseq_validate_t_boundary_statuses(zk0Results, boundary = "H0")
    .bfseq_validate_t_boundary_statuses(zk1Results, boundary = "H1")
    .bfseq_warn_t_boundary_statuses(zk0Results, zk1Results,
                                    tail.eps = tail.eps)
    zk0 <- sapply(zk0Results, `[[`, "value")
    zk1 <- sapply(zk1Results, `[[`, "value")
    bounds <- lapply(seq_along(n1), function(i) {
        list(n1 = n1[[i]], n2 = n2[[i]], se = se[[i]],
             zk0 = zk0Results[[i]]$value,
             zk1 = zk1Results[[i]]$value)
    })
    if (alternative == "two.sided" && strict) {
        regionCount <- .count_strict_two_sided_regions(zk0)
        if (is.infinite(regionCount$total) || regionCount$total > 1000) {
            nregions <- if (is.finite(regionCount$total)) {
                format(regionCount$total, big.mark = ",", scientific = FALSE,
                       trim = TRUE)
            } else {
                "more than 1e308"
            }
            firstH0 <- if (is.na(regionCount$firstH0)) {
                "no finite H0 boundary"
            } else {
                paste0("first finite H0 boundary at look ", regionCount$firstH0)
            }
            warning(paste0(
                "strict = TRUE with two-sided sequential t testing will ",
                "integrate ", nregions, " regions across ", length(n1),
                " looks (", firstH0, "); this can be slow. Consider ",
                "strict = FALSE for the sign-preserving approximation."
            ), immediate. = TRUE, call. = FALSE)
        }
    }
    .bfseq_build_t_design(k1 = k1, k0 = k0, bounds = bounds,
                          dpm = dpm, dpsd = dpsd,
                          plocation = plocation, pscale = pscale, pdf = pdf,
                          alternative = alternative, type = type,
                          trange = trange, strict = strict,
                          tail.eps = tail.eps, tail.nquad = tail.nquad,
                          dots = list(...))
}

## ## compare to simulation-based probabilities
## ## TODO implement as real tests
## set.seed(142)
## n <- seq(40, 100, 10)
## dpm <- 0.5
## dpsd <- 0.1
## k0 <- 6
## k1 <- 1/30
## alternative <- "two.sided"
## type <- "two.sample"
## plocation <- 0
## pscale <- 1/sqrt(2)
## pdf <- 1
## nsim <- 10000
## results <- replicate(n = nsim, expr = {
##     smd <- rnorm(n = 1, mean = dpm, sd = dpsd)
##     y1 <- rnorm(n = max(n), mean = 0, sd = 1)
##     y2 <- rnorm(n = max(n), mean = smd, sd = 1)
##     t <- sapply(seq_along(n), FUN = function(i) {
##         ttest <- t.test(y2[1:n[i]], y1[1:n[i]], var.equal = TRUE,
##                         alternative = "two.sided")$statistic
##         ## (mean(y2[1:n[i]]) - mean(y1[1:n[i]]))*sqrt(n[i]/2)
##         })
##     bf <- sapply(seq_along(n), FUN = function(i) {
##         tbf01(t = t[i], n = n[i], plocation = plocation,
##               pscale = pscale, pdf = pdf, type = type,
##               alternative = alternative)
##     })
##     result <- "inconclusive"
##     for (i in seq_along(bf)) {
##         if (bf[i] >= k0) {
##             result <- "H0"
##             break
##         }
##         if (bf[i] <= k1) {
##         ## if ((bf[i] < k1) & t[i] > 0) {
##             result <- "H1"
##             break
##         }
##     }
##     result
## })

## ptbf01seq(k1 = k1, k0 = k0, n1 = n, n2 = n, plocation = plocation,
##           pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd,
##           alternative = alternative, type = type, strict = FALSE)
## mean(results == "H1")
## mean(results == "H0")
## mean(results == "inconclusive")
