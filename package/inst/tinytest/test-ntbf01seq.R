library(tinytest)
library(bfpwr)

k1 <- 1/2
k0 <- 2
pow <- 0.4

search <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, details = TRUE)
)

expect_true(search$reached,
            info = "sequential t search should reach H1 target")
expect_true(search$actualPower >= pow,
            info = "sequential t search achieved power should exceed target")
expect_true(inherits(search$result, "bfseqdesign"),
            info = "sequential t search details should include design object")
expect_equal(search$result$solver$n, search$n,
             info = "sequential t design should carry solver metadata")

h0search <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0, dpsd = 0,
              alternative = "greater", target = "h0", looks = 1,
              nrange = c(2, 80), strict = FALSE, details = TRUE)
)
expect_true(h0search$reached,
            info = "sequential t search should reach H0 target")
expect_equal(h0search$target, "h0",
             info = "sequential t H0 search should retain target label")

powres <- suppressWarnings(
    powertbf01seq(power = pow, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, nrange = c(2, 80),
                  strict = FALSE)
)
expect_true(inherits(powres, "bfseqdesign"),
            info = "powertbf01seq should return a sequential design")
expect_equal(powres$solver$n, search$n,
             info = "powertbf01seq should use ntbf01seq search result")

ratiores <- suppressWarnings(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, ratio = 2,
                  strict = FALSE)
)
expect_equal(ratiores$n2, ceiling(ratiores$n1*2),
             info = "powertbf01seq should apply two-sample allocation ratio")
expect_equal(ratiores$ratio, 2,
             info = "powertbf01seq should retain allocation ratio")
