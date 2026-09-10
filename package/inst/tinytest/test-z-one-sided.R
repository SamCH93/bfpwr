library(tinytest)
library(bfpwr)

## Independent marginal-likelihood integration, including shifted priors whose
## mean is on the excluded side and a nonzero point null.
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (pm in c(-0.5, 0, 0.8)) {
        for (estimate in c(-1, 0.2, 1.5)) {
            null <- 0.2
            se <- 0.7
            psd <- 0.6
            marginal <- integrate(function(distance) {
                theta <- null + direction*distance
                dnorm(estimate, theta, se)*dnorm(theta, pm, psd)
            }, 0, Inf, rel.tol = 1e-11)$value/
                pnorm(direction*(pm - null)/psd)
            actual <- bf01(estimate, se, null, pm, psd,
                            alternative = alternative)
            expect_equal(actual, dnorm(estimate, null, se)/marginal,
                         tolerance = 1e-9)
            expect_equal(actual, bf01(direction*(estimate - null), se,
                         pm = direction*(pm - null), psd = psd,
                         alternative = "greater"), tolerance = 1e-12)
        }
    }
}

## Tail reference via rescaled integrals of exp(a*x - b*x^2/2). Rescaling
## keeps quadrature accurate when the retained prior probability underflows.
logIntegral <- function(a, b) {
    scale <- 1/(sqrt(b) + abs(a))
    log(scale) + log(integrate(function(x) {
        exp(a*scale*x - b*scale^2*x^2/2)
    }, 0, Inf, rel.tol = 1e-12)$value)
}
for (case in list(c(-1000, 0, 1), c(-100, -30, 1), c(5, -40, 1),
                 c(1, -100, 0.1), c(10, -1, 1e-6))) {
    z <- case[1]
    m <- case[2]
    r <- case[3]
    reference <- logIntegral(m/r^2, 1/r^2) -
        logIntegral(m/r^2 + z, 1 + 1/r^2)
    actual <- bf01(z, se = 1, pm = m, psd = r,
                    alternative = "greater", log = TRUE)
    expect_true(is.finite(actual))
    expect_true(abs(actual - reference) < 1e-10)
}

## Roots must recover the BF threshold even for very narrow priors. Power
## uses the corresponding single tail and preserves its small complement.
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (psd in c(1e-12, 1e-4, 0.5)) {
        for (k in c(0.1, 1, 10)) {
            pm <- direction*0.5
            critical <- bfpwr:::zcrit(k, se = 0.1, mu = pm, tau = psd,
                                      alternative = alternative)
            expect_equal(bf01(0.1*critical, se = 0.1, pm = pm, psd = psd,
                              log = TRUE, alternative = alternative),
                         log(k), tolerance = 1e-8)
            actual <- pbf01(k, n = 100, usd = 1, pm = pm, psd = psd,
                             dpm = direction*0.3, dpsd = 0.2,
                             alternative = alternative)
            expect_equal(actual, pnorm(0.1*critical, mean = direction*0.3,
                         sd = sqrt(0.1^2 + 0.2^2),
                         lower.tail = alternative == "less"), tolerance = 1e-10)
        }
    }
    n <- nbf01(k = 0.1, power = 0.8, usd = 1, pm = 0, psd = 1,
                dpm = direction*0.5, dpsd = 0, alternative = alternative)
    achieved <- pbf01(k = 0.1, n = c(n - 1, n), usd = 1, pm = 0, psd = 1,
                      dpm = direction*0.5, dpsd = 0, alternative = alternative)
    expect_true(achieved[1] < 0.8 && achieved[2] >= 0.8)
}
tails <- pbf01(k = 0.1, n = 100, usd = 1, pm = 0, psd = 1,
                dpm = 2, dpsd = 0, lower.tail = c(TRUE, FALSE),
                alternative = "greater")
expect_true(tails[2] > 0 && tails[2] < 1e-50)
expect_equal(sum(tails), 1)

## A narrow excluded-side normal approaches a point immediately above zero.
## BF01 = 1 therefore crosses close to z = 0, even when BF01 - 1 is tiny.
for (psd in c(1e-6, 1e-12)) {
    expect_equal(pbf01(k = 1, n = 1, usd = 1, pm = -1, psd = psd,
                        dpm = 0, dpsd = 0, alternative = "greater"),
                 0.5, tolerance = 1e-8)
}

## A retained point alternative is unchanged; an excluded point has zero
## prior mass and cannot define a one-sided Bayes factor.
expect_equal(bf01(0.2, 0.1, pm = 0.5, psd = 0, alternative = "greater"),
             bf01(0.2, 0.1, pm = 0.5, psd = 0))
for (pm in c(-0.5, 0)) {
    expect_error(bf01(0.2, 0.1, pm = pm, psd = 0, alternative = "greater"),
                 "point prior")
    expect_error(pbf01(0.1, 100, 1, pm = pm, psd = 0, alternative = "greater"),
                 "point prior")
    expect_error(nbf01(0.1, 0.8, 1, pm = pm, psd = 0, alternative = "greater"),
                 "point prior")
}

## Plotting must retain the one-sided analysis for both design and null curves.
design <- powerbf01(n = 50, pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
                     alternative = "greater")
expect_equal(design$alternative, "greater")
curves <- plot(design, plot = FALSE)
expect_true(any(grepl("greater", capture.output(print(design)))))
for (target in c("H1", "H0")) {
    curve <- curves[[paste0("powDF", target)]]
    expect_equal(curve$power,
                 pbf01(k = if (target == "H1") design$k else 1/design$k,
                        n = curve$n, usd = sqrt(2), pm = 0, psd = 1,
                        dpm = ifelse(curve$prior == "Design prior", 0.5, 0),
                        dpsd = 0, lower.tail = target == "H1",
                        alternative = "greater"))
}
expect_equal(bf01(0.2, 0.1, pm = 0, psd = 1),
             bf01(0.2, 0.1, pm = 0, psd = 1, alternative = "two.sided"))
