## Small log-scale utilities used where ordinary tail probabilities can
## underflow or where 1 - p would lose all meaningful digits.

.bfpwr_log1pexp <- function(x) {
    ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
}

.bfpwr_logspace_sub <- function(logx, logy) {
    if (logy > logx) {
        return(NaN)
    }
    if (logy == -Inf) {
        return(logx)
    }
    if (logx == logy) {
        return(-Inf)
    }
    logx + log1p(-exp(logy - logx))
}

.bfpwr_logspace_sum <- function(logx) {
    logx <- logx[is.finite(logx)]
    if (length(logx) == 0) {
        return(-Inf)
    }
    m <- max(logx)
    m + log(sum(exp(logx - m)))
}

.bfpwr_lpnorm_interval <- function(lower, upper, mean = 0, sd = 1) {
    ## log P(lower < X < upper), computed from the more stable tail.
    if (is.infinite(lower) && lower < 0 && is.infinite(upper) && upper > 0) {
        return(0)
    }
    if (is.infinite(lower) && lower < 0) {
        return(stats::pnorm(q = upper, mean = mean, sd = sd,
                            lower.tail = TRUE, log.p = TRUE))
    }
    if (is.infinite(upper) && upper > 0) {
        return(stats::pnorm(q = lower, mean = mean, sd = sd,
                            lower.tail = FALSE, log.p = TRUE))
    }

    log_lower <- stats::pnorm(q = lower, mean = mean, sd = sd,
                              lower.tail = TRUE, log.p = TRUE)
    log_upper <- stats::pnorm(q = upper, mean = mean, sd = sd,
                              lower.tail = TRUE, log.p = TRUE)
    lower_scale <- .bfpwr_logspace_sub(log_upper, log_lower)

    log_lower_tail <- stats::pnorm(q = lower, mean = mean, sd = sd,
                                   lower.tail = FALSE, log.p = TRUE)
    log_upper_tail <- stats::pnorm(q = upper, mean = mean, sd = sd,
                                   lower.tail = FALSE, log.p = TRUE)
    upper_scale <- .bfpwr_logspace_sub(log_lower_tail, log_upper_tail)

    max(lower_scale, upper_scale, na.rm = TRUE)
}

.bfpwr_lpbeta_interval <- function(lower, upper, shape1, shape2) {
    ## log P(lower < X < upper), computed from the more stable tail.
    if (lower <= 0 && upper >= 1) {
        return(0)
    }
    if (lower <= 0) {
        return(stats::pbeta(q = upper, shape1 = shape1, shape2 = shape2,
                            lower.tail = TRUE, log.p = TRUE))
    }
    if (upper >= 1) {
        return(stats::pbeta(q = lower, shape1 = shape1, shape2 = shape2,
                            lower.tail = FALSE, log.p = TRUE))
    }

    log_lower <- stats::pbeta(q = lower, shape1 = shape1, shape2 = shape2,
                              lower.tail = TRUE, log.p = TRUE)
    log_upper <- stats::pbeta(q = upper, shape1 = shape1, shape2 = shape2,
                              lower.tail = TRUE, log.p = TRUE)
    lower_scale <- .bfpwr_logspace_sub(log_upper, log_lower)

    log_lower_tail <- stats::pbeta(q = lower, shape1 = shape1, shape2 = shape2,
                                   lower.tail = FALSE, log.p = TRUE)
    log_upper_tail <- stats::pbeta(q = upper, shape1 = shape1, shape2 = shape2,
                                   lower.tail = FALSE, log.p = TRUE)
    upper_scale <- .bfpwr_logspace_sub(log_lower_tail, log_upper_tail)

    max(lower_scale, upper_scale, na.rm = TRUE)
}

.bfpwr_qnorm_logistic_inverse <- function(log_odds) {
    ## qnorm(p), where p = 1/(1 + exp(log_odds)).
    if (log_odds >= 0) {
        logp <- -.bfpwr_log1pexp(log_odds)
        stats::qnorm(p = logp, lower.tail = TRUE, log.p = TRUE)
    } else {
        log_upper <- log_odds - .bfpwr_log1pexp(log_odds)
        stats::qnorm(p = log_upper, lower.tail = FALSE, log.p = TRUE)
    }
}
