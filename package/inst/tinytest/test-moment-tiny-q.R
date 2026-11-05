library(tinytest)
library(bfpwr)

## Leading-order BF01 - 1 independently resolves narrow-prior effects.
r <- 1e-8
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (z in c(-2, -1, 1, 2, 10)) {
        actual <- expm1(nmbf01(direction*z, se = 1, psd = r,
                               alternative = alternative, log = TRUE))
        reference <- 1.5*r^2 - 4*dnorm(0)*z*r
        expect_equal(actual/reference, 1, tolerance = 1e-6)
    }
    rCritical <- 5e-8
    critical <- bfpwr:::zcrit(1, se = 1, tau = rCritical, type = "moment",
                              alternative = alternative)
    reference <- 1.5*rCritical*sqrt(1 + rCritical^2)/(4*dnorm(0))
    expect_true(abs(direction*critical - reference) <= 2e-10)
}

## Independent quadrature checks both sides of the Taylor transition.
r <- 0.5
for (q in c(-1e-5*(1 + 1e-8), -1e-5*(1 - 1e-8),
             1e-5*(1 - 1e-8), 1e-5*(1 + 1e-8))) {
    integral <- integrate(function(t) t^2*exp(q*t - t^2/2),
                           0, Inf, rel.tol = 1e-12)$value
    reference <- 1.5*log1p(r*r) + log(2*pi)/2 - log(2) - log(integral)
    for (alternative in c("greater", "less")) {
        direction <- if (alternative == "greater") 1 else -1
        actual <- nmbf01(direction*q*sqrt(1 + r*r)/r, se = 1, psd = r,
                         alternative = alternative, log = TRUE)
        expect_true(abs(actual - reference) < 1e-10)
    }
}
