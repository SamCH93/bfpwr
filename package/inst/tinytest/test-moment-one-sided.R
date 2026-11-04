library(tinytest)
library(bfpwr)

## Independent integration of the normalized half-line moment density.
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (psd in c(0.2, 1)) {
        for (estimate in c(-1, 0.2, 1.5)) {
            null <- 0.2
            se <- 0.7
            marginal <- integrate(function(distance) {
                theta <- null + direction*distance
                dnorm(estimate, theta, se)*2*dnorm(distance, 0, psd)*
                    (distance/psd)^2
            }, 0, Inf, rel.tol = 1e-11)$value
            actual <- nmbf01(estimate, se, null, psd,
                              alternative = alternative)
            expect_equal(actual, dnorm(estimate, null, se)/marginal,
                         tolerance = 1e-9)
            expect_equal(log(actual), nmbf01(estimate, se, null, psd,
                         log = TRUE, alternative = alternative))
            expect_equal(actual, nmbf01(direction*(estimate - null), se,
                         psd = psd, alternative = "greater"))
        }
    }
}

## Rescaled likelihood-ratio integration stays accurate in the wrong tail,
## where the ordinary normal-CDF formula would subtract nearly equal terms.
for (q in c(-1e6, -100, -20, -5.001, -5, -4.999, -1, 0, 1)) {
    for (r in c(1e-6, 0.5, 10)) {
        scale <- 1/(1 + abs(q))
        logIntegral <- 3*log(scale) + log(integrate(function(t) {
            t^2*exp(q*scale*t - (scale*t)^2/2)
        }, 0, Inf, rel.tol = 1e-12)$value)
        reference <- 1.5*log1p(r^2) + log(2*pi)/2 - log(2) - logIntegral
        actual <- nmbf01(q*sqrt(1 + 1/r^2), se = 1, psd = r,
                          alternative = "greater", log = TRUE)
        expect_true(is.finite(actual))
        expect_true(abs(actual - reference) < 1e-10)
    }
}

## The two half-line marginals average to the original two-sided marginal.
for (estimate in c(-3, 0, 2)) {
    halves <- nmbf01(estimate, se = 0.5, psd = 0.7,
                     alternative = c("less", "greater"))
    expect_equal(mean(1/halves), 1/nmbf01(estimate, se = 0.5, psd = 0.7))
}
expect_equal(nmbf01(0, 1, psd = 1, alternative = "greater"),
             nmbf01(0, 1, psd = 1))
expect_true(all(diff(nmbf01(seq(-20, 20, length.out = 201), 1, psd = 1,
                           alternative = "greater", log = TRUE)) < 0))

## Boundaries recover the requested BF, and probabilities use the correct
## predictive tail for both point and untruncated normal design priors.
for (alternative in c("less", "greater")) {
    direction <- if (alternative == "greater") 1 else -1
    for (psd in c(1e-8, 0.2, 2)) {
        for (k in c(0.1, 1, 10)) {
            critical <- bfpwr:::zcrit(k, se = 0.2, null = 2, tau = psd,
                                      type = "moment", alternative = alternative)
            expect_equal(nmbf01(critical*0.2, se = 0.2, psd = psd,
                                alternative = alternative, log = TRUE),
                         log(k), tolerance = 1e-8)
            for (dpsd in c(0, 0.3)) {
                actual <- pnmbf01(k, n = 25, usd = 1, null = 2, psd = psd,
                    dpm = 2 + direction*0.4, dpsd = dpsd,
                    alternative = alternative, lower.tail = c(TRUE, FALSE))
                reference <- pnorm((0.2*critical - direction*0.4)/
                    sqrt(0.2^2 + dpsd^2),
                    lower.tail = alternative == "less")
                reference <- c(reference, pnorm((0.2*critical - direction*0.4)/
                    sqrt(0.2^2 + dpsd^2), lower.tail = alternative != "less"))
                expect_equal(actual, reference, tolerance = 1e-10)
            }
        }
    }
    for (lower.tail in c(TRUE, FALSE)) {
        k <- if (lower.tail) 0.1 else 3
        dpm <- if (lower.tail) direction*0.5 else 0
        n <- nnmbf01(k, power = 0.8, usd = 1, psd = 0.7, dpm = dpm,
                     dpsd = 0, lower.tail = lower.tail, alternative = alternative)
        achieved <- pnmbf01(k, n = c(n - 1, n), usd = 1, psd = 0.7,
                            dpm = dpm, dpsd = 0, lower.tail = lower.tail,
                            alternative = alternative)
        expect_true(achieved[1] < 0.8 && achieved[2] >= 0.8)
    }
}
tails <- pnmbf01(0.1, n = 100, usd = 1, psd = 1, dpm = 2, dpsd = 0,
                  alternative = "greater", lower.tail = c(TRUE, FALSE))
expect_true(tails[2] > 0 && tails[2] < 1e-40)
expect_equal(sum(tails), 1)
expect_equal(nnmbf01(0.1, 0.8, usd = 1, psd = 0.7, dpm = c(-0.5, 0.5),
    dpsd = 0, alternative = c("less", "greater")),
    rep(nnmbf01(0.1, 0.8, usd = 1, psd = 0.7, dpm = 0.5, dpsd = 0,
                  alternative = "greater"), 2))

## High-level power, sample size, printing, and both plotted power curves
## retain the requested alternative. Objects saved before this addition work.
for (alternative in c("less", "greater")) {
    direction <- if (alternative == "greater") 1 else -1
    design <- powernmbf01(n = 50, psd = 0.7, dpm = direction*0.5, dpsd = 0,
                           alternative = alternative)
    expect_equal(design$alternative, alternative)
    printed <- capture.output(print(design))
    expect_true(any(grepl(alternative, printed)))
    expect_true(any(grepl("prior parameters are before truncation", printed, fixed = TRUE)))
    curves <- plot(design, ngrid = 10, plot = FALSE)
    for (target in c("H1", "H0")) {
        curve <- curves[[paste0("powDF", target)]]
        expect_equal(curve$power, pnmbf01(
            k = if (target == "H1") design$k else 1/design$k,
            n = curve$n, usd = sqrt(2), psd = 0.7,
            dpm = ifelse(curve$prior == "Design prior", direction*0.5, 0),
            dpsd = 0, lower.tail = target == "H1", alternative = alternative))
    }
    found <- powernmbf01(power = 0.8, psd = 0.7, dpm = direction*0.5,
                          dpsd = 0, alternative = alternative)
    expect_equal(pnmbf01(0.1, found$n, usd = sqrt(2), psd = 0.7,
        dpm = direction*0.5, dpsd = 0, alternative = alternative),
        0.8, tolerance = 1e-8)
}
design <- powernmbf01(n = 50, psd = 1, dpm = 0.5, dpsd = 0)
legacy <- design
legacy$alternative <- NULL
expect_equal(plot(legacy, ngrid = 10, plot = FALSE),
             plot(design, ngrid = 10, plot = FALSE))
expect_equal(pnmbf01(0.1, 50, 1, psd = 1, dpm = 0.5, dpsd = 0),
             pnmbf01(0.1, 50, 1, psd = 1, dpm = 0.5, dpsd = 0,
                     alternative = "two.sided"))
expect_error(nmbf01(0, 1, psd = 1, alternative = "invalid"), "should be one of")
expect_error(pnmbf01(0.1, 50, 1, psd = 1, dpm = 0.5, dpsd = 0,
                     alternative = "invalid"), "should be one of")
expect_error(nmbf01(0, 1, psd = 0, alternative = "greater"), "0 < psd")
