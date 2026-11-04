library(tinytest)
library(bfpwr)

## Finite standardized q survives both small-ratio square/reciprocal limits.
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (r in c(1e-155, 1e-200)) {
        actual <- nmbf01(direction/r, se = 1, psd = r,
                         alternative = alternative, log = TRUE)
        expect_equal(actual, -1.84789662189189, tolerance = 1e-12)
        tiny <- nmbf01(direction, se = 1, psd = r,
                       alternative = alternative, log = TRUE)
        expect_true(tiny < 0)
        expect_equal(tiny/r, -4*dnorm(0), tolerance = 1e-12)
    }
    for (r in c(1e155, 1e200)) {
        actual <- nmbf01(0, se = 1, psd = r,
                         alternative = alternative, log = TRUE)
        expect_true(is.finite(actual))
        expect_equal(actual, 3*log(r), tolerance = 1e-12)
    }
    expect_equal(nmbf01(0, se = 1e-200, psd = 1e200,
        alternative = alternative, log = TRUE), 2763.1021115928552,
        tolerance = 1e-12)

    ## q^2 overflows, although q^2/2 and the log BF are representable.
    q <- 1.5e154
    actual <- nmbf01(direction*q*sqrt(2), se = 1, psd = 1,
                     alternative = alternative, log = TRUE)
    expect_true(is.finite(actual))
    expect_equal(actual/(-1.125e308), 1, tolerance = 1e-14)

    ## Independent recorded q-space quadrature roots; k=1 near zero is
    ## deliberately omitted because its absolute q tolerance amplifies in z.
    for (r in c(1e-160, 1e160)) {
        qref <- if (r < 1) c(1.20763668045903, -0.764891549168737) else
            c(46.8861255416956, 46.8135937888452)
        w <- if (r < 1) r/sqrt(1 + r*r) else 1/sqrt(1 + (1/r)^2)
        for (i in 1:2) {
            k <- c(0.1, 3)[i]
            critical <- bfpwr:::zcrit(k, se = 1, tau = r,
                type = "moment", alternative = alternative)
            expect_true(is.finite(critical))
            expect_equal(direction*critical*w, qref[i], tolerance = 1e-10)
            probability <- pnmbf01(k, n = 1, usd = 1, psd = r,
                dpm = 0, dpsd = 0, alternative = alternative)
            expect_true(is.finite(probability))
            expect_equal(probability, pnorm(qref[i]/w, lower.tail = FALSE),
                         tolerance = 1e-10)
        }
    }
}
