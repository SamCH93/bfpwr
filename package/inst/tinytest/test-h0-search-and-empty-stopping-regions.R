library(tinytest)
library(bfpwr)

## Design-prior mass inside H1 can make H0-evidence probability peak inside
## the search range. Reflection and a shifted null describe the same problem.
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

## Probabilities of misleading evidence rise and then fall for any sidedness,
## tail, prior family, or test. Each search returns the first crossing, which
## the bracket places before the interior maximum.
peakCases <- list(
    list(label = "one-sided z, point design inside H1",
         power = function(n) pbf01(k = 3, n = n, usd = 1, pm = -0.2, psd = 0.3,
                                   dpm = -0.2, dpsd = 0, alternative = "less",
                                   lower.tail = FALSE),
         search = function(target) nbf01(k = 3, power = target, usd = 1,
                                         pm = -0.2, psd = 0.3, dpm = -0.2,
                                         dpsd = 0, alternative = "less",
                                         lower.tail = FALSE),
         target = 0.1, bracket = c(10, 25), unattainable = 0.11),
    list(label = "two-sided z, point design",
         power = function(n) pbf01(k = 3, n = n, usd = 1, pm = 0, psd = 1,
                                   dpm = 0.2, dpsd = 0, lower.tail = FALSE),
         search = function(target) nbf01(k = 3, power = target, usd = 1,
                                         pm = 0, psd = 1, dpm = 0.2, dpsd = 0,
                                         lower.tail = FALSE),
         target = 0.45, bracket = c(10, 25), unattainable = 0.51),
    list(label = "two-sided z, H1 evidence under the null",
         power = function(n) pbf01(k = 1/3, n = n, usd = 1, pm = 0, psd = 1,
                                   dpm = 0, dpsd = 0),
         search = function(target) nbf01(k = 1/3, power = target, usd = 1,
                                         pm = 0, psd = 1, dpm = 0, dpsd = 0),
         target = 0.02, bracket = c(1, 5), unattainable = 0.03),
    list(label = "moment prior",
         power = function(n) pnmbf01(k = 3, n = n, usd = 1, psd = 0.5,
                                     dpm = 0.2, dpsd = 0.3, lower.tail = FALSE),
         search = function(target) nnmbf01(k = 3, power = target, usd = 1,
                                           psd = 0.5, dpm = 0.2, dpsd = 0.3,
                                           lower.tail = FALSE),
         target = 0.4, bracket = c(5, 10), unattainable = 0.51)
)
for (case in peakCases) {
    reference <- uniroot(function(n) case$power(n) - case$target,
                         case$bracket, tol = 1e-10)$root
    n <- case$search(case$target)
    expect_equal(n, ceiling(reference), info = case$label)
    expect_true(case$power(n) >= case$target &&
                case$power(n - 1) < case$target, info = case$label)
    expect_true(is.nan(suppressWarnings(case$search(case$unattainable))),
                info = case$label)
}

## The t-test search follows the same rule as the z-test search.
tpower <- function(n, dpm, dpsd) {
    ptbf01(k = 3, n = n, plocation = -0.2, pscale = 0.3, pdf = 30,
           dpm = dpm, dpsd = dpsd, type = "one.sample", alternative = "less",
           lower.tail = FALSE)
}
tsearch <- function(power, dpm, dpsd) {
    ntbf01(k = 3, power = power, plocation = -0.2, pscale = 0.3, pdf = 30,
           dpm = dpm, dpsd = dpsd, type = "one.sample", alternative = "less",
           lower.tail = FALSE, nrange = c(10, 500))
}
expect_true(tpower(500, 0.2, 0.3) < 0.8 && tpower(250, 0.2, 0.3) > 0.8)
n <- tsearch(0.8, 0.2, 0.3)
expect_equal(n, 176)
expect_true(tpower(n, 0.2, 0.3) >= 0.8 && tpower(n - 1, 0.2, 0.3) < 0.8)
expect_true(is.nan(suppressWarnings(tsearch(0.81, 0.2, 0.3))))
n <- tsearch(0.1, -0.2, 0)
expect_true(n < 50)
expect_true(tpower(n, -0.2, 0) >= 0.1 && tpower(n - 1, -0.2, 0) < 0.1)

## powerbf01seq() names the BF family 'bftype', not 'type'.
for (bftype in c("directional", "moment")) {
    expect_error(powerbf01seq(n = 30, pm = 0, psd = 1, dpm = 0.2, dpsd = 0,
                              bftype = bftype, alternative = "greater"),
                 pattern = "requires bftype = \"normal\"", fixed = TRUE)
    expect_error(powerbf01seq(power = 0.8, pm = 0, psd = 1, dpm = 0.2,
                              dpsd = 0, bftype = bftype, alternative = "less",
                              nrange = c(2, 100)),
                 pattern = "requires bftype = \"normal\"", fixed = TRUE)
}

## Empty intervals have zero mass, including adaptive infinite
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

## Shared integration points preserve physical cumulative probabilities.
n <- seq(10, 500, 10)
design <- pbf01seq(k1 = 1/3, k0 = 3, n = n, se = 1/sqrt(n),
                   pm = 0, psd = 0.1, dpm = 0.2, dpsd = 0,
                   alternative = "less")
probabilities <- cbind(design$cumpH0, design$cumpH1, design$cumpInc)
expect_true(all(probabilities >= -8*.Machine$double.eps &
                probabilities <= 1 + 8*.Machine$double.eps))
expect_equal(rowSums(probabilities), rep(1, length(n)), tolerance = 1e-14)
expect_true(all(diff(design$cumpH0) >= -8*.Machine$double.eps) &&
            all(diff(design$cumpH1) >= -8*.Machine$double.eps))
