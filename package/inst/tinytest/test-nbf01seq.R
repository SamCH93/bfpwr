library(tinytest)
library(bfpwr)

k1 <- 1/5
k0 <- 5
pow <- 0.8
usd <- sqrt(2)
pm <- 0
psd <- 1
dpm <- 0.5
dpsd <- 0

search <- nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd, pm = pm,
                   psd = psd, dpm = dpm, dpsd = dpsd, looks = 3,
                   nrange = c(2, 200), details = TRUE)

expect_true(search$reached,
            info = "sequential z search should reach H1 target")
expect_true(search$actualPower >= pow,
            info = "sequential z search achieved power should exceed target")
expect_equal(search$n, max(search$result$n),
             info = "sequential z search n should be the final look")
expect_equal(search$result$solver$n, search$n,
             info = "sequential z design should carry solver metadata")

schedule <- bfpwr:::.bfseq_schedule_spec(looks = 3, nrange = c(2, 200))
if (search$n > search$nrange[1]) {
    nprev <- bfpwr:::.bfseq_schedule_n(search$n - 1, schedule)
    prev <- pbf01seq(k1 = k1, k0 = k0, se = usd/sqrt(nprev), n = nprev,
                     pm = pm, psd = psd, dpm = dpm, dpsd = dpsd)
    expect_true(utils::tail(prev$cumpH1, 1) < pow,
                info = "sequential z search should find first H1 crossing")
}

h0search <- nbf01seq(k1 = k1, k0 = k0, power = 0.7, usd = usd, pm = pm,
                     psd = psd, dpm = 0, dpsd = 0, target = "h0",
                     looks = 2, nrange = c(2, 500), details = TRUE)
expect_true(h0search$reached,
            info = "sequential z search should reach H0 target")
expect_true(h0search$actualPower >= 0.7,
            info = "sequential z H0 search achieved power should exceed target")
expect_equal(h0search$target, "h0",
             info = "sequential z H0 search should retain target label")

powres <- powerbf01seq(power = pow, k1 = k1, k0 = k0, pm = pm, psd = psd,
                       dpm = dpm, dpsd = dpsd, looks = 3,
                       nrange = c(2, 200))
expect_true(inherits(powres, "bfseqdesign"),
            info = "powerbf01seq should return a sequential design")
expect_equal(powres$solver$n, search$n,
             info = "powerbf01seq should use nbf01seq search result")
expect_equal(powres$sample.type, "two.sample",
             info = "powerbf01seq should retain sampling design")

fixed <- powerbf01seq(n = 60, k1 = k1, k0 = k0, pm = pm, psd = psd,
                      dpm = dpm, dpsd = dpsd, looks = 3)
expect_true(inherits(fixed, "bfseqdesign"),
            info = "powerbf01seq fixed-n mode should return a sequential design")
expect_equal(max(fixed$n), 60,
             info = "powerbf01seq fixed-n mode should use requested final n")
