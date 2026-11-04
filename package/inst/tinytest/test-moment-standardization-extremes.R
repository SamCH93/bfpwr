library(tinytest)
library(bfpwr)

for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (scale in c(1, 1e308)) {
        actual <- nmbf01(direction*scale, se = scale, null = -direction*scale,
            psd = scale, alternative = alternative, log = TRUE)
        expect_equal(actual, -1.7425285997906004, tolerance = 1e-12)
    }

    ## z overflows although the marginal-likelihood q remains finite.
    aligned <- nmbf01(direction*1e210, se = 1e-100, psd = 1e-300,
        alternative = alternative, log = TRUE)
    expect_true(is.finite(aligned))
    expect_equal(aligned/(-5e219), 1, tolerance = 1e-10)
    opposite <- nmbf01(-direction*1e210, se = 1e-100, psd = 1e-300,
        alternative = alternative, log = TRUE)
    expect_true(is.finite(opposite))
    expect_true(abs(opposite - 759.38572486011992) <= 1e-9)

    ## Multiplication in this order supplies an independent tiny-q reference.
    for (psd in c(1e-320, 3e-320)) {
        for (sign in c(-1, 1)) {
            estimate <- direction*sign*1e308
            qref <- (direction*estimate*psd)/(1e4*1e4)
            actual <- nmbf01(estimate, se = 1e4, psd = psd,
                alternative = alternative, log = TRUE)
            expect_equal((actual/qref)/(-4*dnorm(0)), 1, tolerance = 1e-10)
        }
    }

    ## Even q is outside double range, but its wrong-tail log BF is finite.
    reference <- 3*(log(1) - log(1e-308)) + .5*log(2*pi) -
        2*log(2) + 3*(log(1e308) - log(1e-308))
    actual <- nmbf01(-direction*1e308, se = 1e-308, psd = 1,
        alternative = alternative, log = TRUE)
    expect_true(is.finite(actual))
    expect_true(abs(actual - reference) <= 1e-9)
}
