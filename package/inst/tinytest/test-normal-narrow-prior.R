library(tinytest)
library(bfpwr)

## Check the finite boundary against a directly bracketed BF root. Its power
## should approach the point-alternative limit as the normal prior narrows.
for (shift in c(-0.5, 0.5)) {
    for (psd in c(1e-4, 1e-8, 1e-12)) {
        for (k in c(0.1, 3)) {
            null <- 2
            pm <- null + shift
            root <- uniroot(function(z) {
                bf01(estimate = null + 0.1*z, se = 0.1, null = null,
                     pm = pm, psd = psd, log = TRUE) - log(k)
            }, interval = if (shift > 0) c(0, 5) else c(-5, 0),
            tol = 1e-10)$root
            critical <- bfpwr:::zcrit(k = k, se = 0.1, null = null,
                                      mu = pm, tau = psd)
            finiteRoot <- critical[which.min(abs(critical))]
            expect_equal(finiteRoot, root, tolerance = 1e-8)

            actual <- pbf01(k = k, n = 100, usd = 1, null = null,
                            pm = pm, psd = psd, dpm = pm, dpsd = 0)
            reference <- pnorm(root, mean = shift/0.1,
                                lower.tail = shift < 0)
            expect_equal(actual, reference, tolerance = 1e-8)
        }
        sequential <- pbf01seq(k1 = 0.1, k0 = 3, se = 0.1, null = null,
                               pm = pm, psd = psd, dpm = pm, dpsd = 0)
        point <- pbf01seq(k1 = 0.1, k0 = 3, se = 0.1, null = null,
                          pm = pm, psd = 0, dpm = pm, dpsd = 0)
        expect_true(abs(sequential$cumpH1 - point$cumpH1) < 1e-6)
        expect_true(abs(sequential$cumpH0 - point$cumpH0) < 1e-6)
    }
}

for (psd in c(0, 1e-12)) {
    found <- nbf01seq(k1 = 0.1, k0 = 3, power = 0.8, usd = 1,
                      pm = 0.5, psd = psd, dpm = 0.5, dpsd = 0,
                      nrange = c(10, 150), minN = 10, by = 10,
                      details = TRUE)
    expect_true(found$reached)
    if (psd == 0) reference <- found$n
    else expect_equal(found$n, reference,
                      info = "narrow-prior search approaches the point prior")
}
