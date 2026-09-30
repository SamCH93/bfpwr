library(tinytest)
library(bfpwr)

## A continuous design prior can make H0 power peak inside the search range.
## Reflection and a shifted null describe the same sample-size problem.
for (direction in c(-1, 1)) {
    args <- list(k = 3, usd = 1, null = 2, pm = 2 + direction*0.2,
                 psd = 0.3, dpm = 2 - direction*0.2, dpsd = 0.3,
                 alternative = if (direction == 1) "greater" else "less",
                 lower.tail = FALSE)
    power <- function(n) do.call(pbf01, c(args, list(n = n)))
    expect_true(power(500) < 0.8 && power(250) > 0.8)
    reference <- uniroot(function(n) power(n) - 0.8, c(100, 250), tol = 1e-10)$root
    n <- do.call(nbf01, c(args, list(power = 0.8, nrange = c(10, 500))))
    expect_equal(n, ceiling(reference))
    expect_true(power(n) >= 0.8 && power(n - 1) < 0.8)
    fractional <- do.call(nbf01, c(args,
        list(power = 0.8, nrange = c(10, 500), integer = FALSE)))
    expect_equal(fractional, reference, tolerance = 1e-7)
    expect_true(is.nan(suppressWarnings(do.call(nbf01, c(args,
        list(power = 0.81, nrange = c(10, 500)))))))
    expect_true(is.nan(suppressWarnings(do.call(nbf01, c(args,
        list(power = 0.8, nrange = c(10, 100)))))))
}

## Independent coordinates give an exact rare two-look probability, including
## probabilities far below machine epsilon and reflected upper-tail events.
for (sign in c(-1, 1)) {
    region <- if (sign == -1) rbind(c(-1, -Inf), c(1, -10)) else
        rbind(c(-1, 10), c(1, Inf))
    expected <- (pnorm(1) - pnorm(-1))*pnorm(-10)
    actual <- bfpwr:::.bfseq_intstage(
        list(H0 = list(region)), mean = c(0, 0), sigma = diag(2))[1]
    expect_equal(unname(actual)/expected, 1, tolerance = 1e-10)
}

## Empty intervals remain empty after reordering, including adaptive infinite
## boundaries. Sending [Inf, Inf] to the sampler can instead yield NaN.
for (bound in c(-Inf, 0, Inf)) {
    region <- rbind(c(-1, bound), c(1, bound))
    actual <- bfpwr:::.bfseq_intstage(
        list(H0 = list(region)), mean = c(0, 0), sigma = diag(2))[1]
    expect_equal(unname(actual), 0)
}

## The full corpus case has empty H1 events during adaptive boundary searches.
design <- suppressWarnings(ptbf01seq(
    k1 = 1/30, k0 = 30, n = seq(20, 110, 10),
    plocation = 0.5 - 0.2, pscale = 0.35, pdf = 30,
    dpm = 0 - 0.2, dpsd = 0, alternative = 'greater', type = 'one.sample'))
expect_true(all(is.finite(design$cumpH0)) && all(is.finite(design$cumpH1)))
expect_true(is.finite(design$EN1))

## Conditional-normal quadrature independently checks a correlated rare event.
rho <- 0.5
tail <- pnorm(-5.7)
conditional <- function(u) {
    mean <- rho*qnorm(u*tail)
    pnorm(3, mean, sqrt(1 - rho^2)) - pnorm(-3, mean, sqrt(1 - rho^2))
}
expected <- tail*integrate(conditional, 0, 1, rel.tol = 1e-10)$value
region <- rbind(c(-3, -Inf), c(3, -5.7))
actual <- bfpwr:::.bfseq_intstage(list(H0 = list(region)), c(0, 0),
    matrix(c(1, rho, rho, 1), 2))[1]
expect_true(abs(unname(actual)/expected - 1) < 2e-4)

## The corpus regression must retain a rare H0 event; cumulative first-stop
## probability is bounded by the sum of the marginal fixed-look probabilities.
n <- seq(10, 110, 10)
design <- pbf01seq(k1 = 1/30, k0 = 30, n = n, se = 1/sqrt(n),
                   pm = 0, psd = 0.1, dpm = 0, dpsd = 0.5, alternative = "less")
fixed <- pbf01(k = 30, n = n, usd = 1, pm = 0, psd = 0.1,
               dpm = 0, dpsd = 0.5, alternative = "less", lower.tail = FALSE)
expect_true(tail(design$cumpH0, 1) > 1e-6)
expect_true(tail(design$cumpH0, 1) <= sum(fixed))
fine <- pbf01seq(k1 = 1/30, k0 = 30, n = n, se = 1/sqrt(n),
                 pm = 0, psd = 0.1, dpm = 0, dpsd = 0.5,
                 alternative = "less", ngrid = 100000)
expect_true(abs(tail(design$cumpH0, 1)/tail(fine$cumpH0, 1) - 1) < 0.001)
