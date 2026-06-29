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
    if (any(is.na(logx) | is.nan(logx) | logx == Inf)) {
        stop("log-space sum received invalid non-finite values",
             call. = FALSE)
    }
    logx <- logx[logx != -Inf]
    if (length(logx) == 0) {
        return(-Inf)
    }
    m <- max(logx)
    m + log(sum(exp(logx - m)))
}

.bfpwr_lpnorm_interval <- function(lower, upper, mean = 0, sd = 1) {
    ## log P(lower < X < upper), computed from the more stable tail.
    n <- max(length(lower), length(upper), length(mean), length(sd))
    lower <- rep(lower, length.out = n)
    upper <- rep(upper, length.out = n)
    mean <- rep(mean, length.out = n)
    sd <- rep(sd, length.out = n)

    ans <- rep(NA_real_, n)
    lower_inf <- is.infinite(lower) & lower < 0
    upper_inf <- is.infinite(upper) & upper > 0

    all_real <- lower_inf & upper_inf
    ans[all_real] <- 0
    left_tail <- lower_inf & !upper_inf
    ans[left_tail] <- stats::pnorm(q = upper[left_tail],
                                   mean = mean[left_tail], sd = sd[left_tail],
                                   lower.tail = TRUE, log.p = TRUE)
    right_tail <- !lower_inf & upper_inf
    ans[right_tail] <- stats::pnorm(q = lower[right_tail],
                                    mean = mean[right_tail],
                                    sd = sd[right_tail],
                                    lower.tail = FALSE, log.p = TRUE)

    finite_interval <- is.na(ans)
    if (any(finite_interval)) {
        log_lower <- stats::pnorm(q = lower[finite_interval],
                                  mean = mean[finite_interval],
                                  sd = sd[finite_interval],
                                  lower.tail = TRUE, log.p = TRUE)
        log_upper <- stats::pnorm(q = upper[finite_interval],
                                  mean = mean[finite_interval],
                                  sd = sd[finite_interval],
                                  lower.tail = TRUE, log.p = TRUE)
        lower_scale <- log_upper + log1p(-exp(log_lower - log_upper))
        lower_scale[log_lower == log_upper] <- -Inf

        log_lower_tail <- stats::pnorm(q = lower[finite_interval],
                                       mean = mean[finite_interval],
                                       sd = sd[finite_interval],
                                       lower.tail = FALSE, log.p = TRUE)
        log_upper_tail <- stats::pnorm(q = upper[finite_interval],
                                       mean = mean[finite_interval],
                                       sd = sd[finite_interval],
                                       lower.tail = FALSE, log.p = TRUE)
        upper_scale <- log_lower_tail +
            log1p(-exp(log_upper_tail - log_lower_tail))
        upper_scale[log_lower_tail == log_upper_tail] <- -Inf

        ans[finite_interval] <- pmax(lower_scale, upper_scale, na.rm = TRUE)
    }

    if (n == 1) unname(ans) else ans
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

.bfpwr_gauss_legendre <- local({
    ## Cache nodes by order because t critical-value searches repeatedly reuse
    ## the same quadrature rule. Keeping the cache in this closure avoids a
    ## separate package-level mutable object.
    cache <- new.env(parent = emptyenv())

    function(n) {
        key <- as.character(n)
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }

        i <- seq_len(n - 1)
        beta <- i/sqrt(4*i^2 - 1)
        J <- matrix(0, nrow = n, ncol = n)
        J[cbind(i, i + 1)] <- beta
        J[cbind(i + 1, i)] <- beta
        eig <- eigen(J, symmetric = TRUE)
        o <- order(eig$values)
        ans <- list(
            x = (eig$values[o] + 1)/2,
            w = eig$vectors[1, o]^2
        )
        assign(key, ans, envir = cache)
        ans
    }
})
