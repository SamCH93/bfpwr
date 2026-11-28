library(tinytest)
library(bfpwr)

## Under a normal prior, the predictive t statistic is a scaled noncentral t.
n <- 30
t <- 2.5
m <- 0.2
s <- 0.5
scale <- sqrt(1 + n*s^2)
expected <- dt(t, df = n - 1)/(dt(t/scale, df = n - 1,
                                   ncp = sqrt(n)*m/scale)/scale)
expect_equal(tbf01(t = t, n = n, plocation = m, pscale = s, pdf = Inf,
                   type = "one.sample"), expected, tolerance = 1e-8)

for (alternative in c("less", "greater")) {
    bf <- tbf01(t = t, n = n, plocation = m, pscale = s, pdf = c(Inf, 1e8),
                type = "one.sample", alternative = alternative)
    expect_equal(bf[1], bf[2], tolerance = 1e-6, info = alternative)
}
for (pdf in c(-Inf, -1, 0, NA, NaN)) {
    expect_error(tbf01(t = t, n = n, pdf = pdf))
}

## Fixed and sequential designs should agree with their large-df limits.
power <- ptbf01(k = 1/6, n = 50, pdf = c(Inf, 1e8), dpm = 0.5, dpsd = 0)
expect_equal(power[1], power[2], tolerance = 1e-6)
expect_equal(powertbf01(k = 1/6, n = 50, pdf = Inf, dpm = 0.5, dpsd = 0)$power,
             power[2], tolerance = 1e-6)

nfixed <- ntbf01(k = 1/6, power = 0.5, pdf = c(Inf, 1e8), dpm = 0.5,
                 dpsd = 0, nrange = c(30, 100))
expect_true(is.finite(nfixed[1]))
expect_equal(nfixed[1], nfixed[2])

seq_normal <- ptbf01seq(k1 = 1/6, k0 = 2, n = c(30, 50), pdf = Inf,
                        dpm = 0.5, dpsd = 0)
seq_t <- ptbf01seq(k1 = 1/6, k0 = 2, n = c(30, 50), pdf = 1e8,
                   dpm = 0.5, dpsd = 0)
expect_equal(seq_normal$cumpH1, seq_t$cumpH1, tolerance = 1e-6)
expect_equal(seq_normal$cumpH0, seq_t$cumpH0, tolerance = 1e-6)

nseq <- ntbf01seq(k1 = 1/6, k0 = 2, power = 0.5, pdf = c(Inf, 1e8),
                  dpm = 0.5, dpsd = 0, nrange = c(30, 100))
expect_true(is.finite(nseq[1]))
expect_equal(nseq[1], nseq[2])
