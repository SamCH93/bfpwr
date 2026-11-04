library(tinytest)
library(bfpwr)

## Independent q-space roots and corresponding predictive normal probabilities.
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    for (tau in c(1e-300, 1e-308)) {
        for (i in 1:2) {
            k <- c(.1, 3)[i]
            qref <- c(1.2076366804590342, -.76489154916873636)[i]
            pref <- c(.1135935218056585, .77783193672655049)[i]
            critical <- bfpwr:::zcrit(k, se = 1, tau = tau,
                type = "moment", alternative = alternative)
            expect_true(is.finite(critical))
            expect_equal(direction*critical*tau, qref, tolerance = 1e-10)
            for (lower.tail in c(TRUE, FALSE)) {
                actual <- pnmbf01(k, n = 1, usd = 1, psd = tau, dpm = 0,
                    dpsd = 1/tau, alternative = alternative,
                    lower.tail = lower.tail)
                expect_equal(actual, if (lower.tail) pref else 1 - pref,
                    tolerance = 1e-10)
            }
        }
    }
    expect_equal(bfpwr:::zcrit(2^1.5, se = 1, tau = 1,
        type = "moment", alternative = alternative), 0, tolerance = 1e-10)
    expect_error(bfpwr:::zcrit(.1, se = 1, tau = 1e-320,
        type = "moment", alternative = alternative), "finite z range")
}
