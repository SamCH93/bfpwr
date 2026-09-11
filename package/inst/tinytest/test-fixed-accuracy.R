library(tinytest)
library(bfpwr)

## Recover a known fractional sample size, so rounding cannot hide root error.
zargs <- list(k = 0.1, usd = 1, pm = 0, psd = 1/sqrt(2),
              dpm = 0.5, dpsd = 0, alternative = "greater")
target <- do.call(pbf01, c(list(n = 50.5), zargs))
zsearch <- c(list(power = target, integer = FALSE, nrange = c(1, 1000)), zargs)
expect_true(abs(do.call(nbf01, zsearch) - 50.5) < 1e-7)
coarse <- do.call(nbf01, c(zsearch, list(tol = 0.01)))
expect_true(abs(coarse - 50.5) > 1e-7,
            info = "callers can trade sample-size accuracy for fewer iterations")
wrapped <- powerbf01(power = target, type = "one.sample", pm = 0,
                      psd = 1/sqrt(2), dpm = 0.5, dpsd = 0,
                      alternative = "greater", nrange = c(1, 1000), tol = 0.01)
expect_equal(wrapped$n, coarse)
expect_equal(wrapped$tol, 0.01)

moment <- powernmbf01(power = 0.8, psd = 1/sqrt(2), dpm = 0.5, dpsd = 0,
                      type = "one.sample", nrange = c(1, 1000), tol = 0.01)
expect_equal(moment$n, nnmbf01(k = 0.1, power = 0.8, usd = 1,
    psd = 1/sqrt(2), dpm = 0.5, dpsd = 0, nrange = c(1, 1000),
    integer = FALSE, tol = 0.01))
expect_equal(moment$tol, 0.01)
binomial <- powerbinbf01(power = 0.8, dp = 0.7, nrange = c(1, 1000), tol = 0.01)
expect_equal(binomial$n, nbinbf01(k = 0.1, power = 0.8, dp = 0.7,
    nrange = c(1, 1000), tol = 0.01))
expect_equal(binomial$tol, 0.01)

## Both adaptive and explicit ranges accept integration and root controls.
## Previously those controls were dropped or leaked into the root function.
targs <- list(k = 0.1, n = 100, dpm = 0.5, dpsd = 0,
              alternative = "greater")
controls <- list(rel.tol = 1e-10, abs.tol = 1e-11,
                 subdivisions = 500, tol = 1e-10)
fine <- do.call(ptbf01, c(targs, controls))
explicit <- do.call(ptbf01, c(targs, controls, list(drange = c(0, 2))))
expect_true(is.finite(explicit))
expect_equal(fine, explicit, tolerance = 1e-9)
expect_equal(do.call(ptbf01, targs), fine, tolerance = 1e-8)
nargs <- targs[names(targs) != "n"]
found <- do.call(ntbf01, c(nargs, controls,
    list(power = fine, type = "one.sample", integer = FALSE, nrange = c(2, 200))))
actual <- do.call(ptbf01, c(nargs, controls, list(n = found, type = "one.sample")))
expect_equal(actual, fine, tolerance = 1e-9)

## Plotting retains the numerical settings chosen for a power result.
wrapped <- do.call(powertbf01, c(targs, controls))
expect_equal(wrapped$numerical[names(controls)], controls)
curves <- suppressWarnings(plot(wrapped, nlim = c(90, 110), ngrid = 3,
                                 plot = FALSE, nullplot = FALSE))$powDFH1
design <- curves[curves$prior == "Design prior", ]
expect_equal(design$power, do.call(ptbf01,
    c(targs[names(targs) != "n"], controls, list(n = design$n))), tolerance = 1e-12)

## Sequential BF controls are separate from the multivariate integral, and
## sample-size searches must use the same boundaries as direct calls.
seqargs <- list(k1 = 0.1, k0 = 3, dpm = 0.5, dpsd = 0,
                alternative = "greater", type = "one.sample",
                bf.control = controls, ngrid = 1000)
seq <- do.call(ptbf01seq, c(seqargs, list(n = c(30, 60))))
expect_equal(seq$integration$bf.control, controls)
expect_equal(seq$cumpH1[1], do.call(ptbf01,
    c(targs[names(targs) != "n"], controls, list(n = 30, type = "one.sample"))),
    tolerance = 1e-8)
search <- do.call(ntbf01seq, c(seqargs, list(power = tail(seq$cumpH1, 1),
    nrange = c(30, 60), minN = 30, by = 30, details = TRUE)))
expect_equal(search$result$zk1, seq$zk1, tolerance = 1e-10)
expect_equal(search$result$integration$bf.control, controls)
expect_error(ptbf01seq(n = c(30, 60), bf.control = list(ngrid = 1000)))

## The stable fallback also converges with its separate node control.
tailargs <- list(t = -10, n = 20, alternative = "greater", log = TRUE)
reference <- do.call(tbf01, c(tailargs, list(tail.nquad = 1024)))
default <- do.call(tbf01, tailargs)
previous <- do.call(tbf01, c(tailargs, list(tail.nquad = 128)))
expect_true(abs(default - reference) < 1e-7)
expect_true(abs(default - reference) < abs(previous - reference)/10)

## A tail event below the old cutoff is retained at the tighter default.
critical <- bfpwr:::tcrit(k = 30, n1 = 41, n2 = 41, plocation = 0,
    pscale = 1/sqrt(2), pdf = 1, type = "two.sample", alternative = "greater",
    trange = c(-20, 0))
se <- sqrt(2/41)
tailpower <- list(k = 30, n = 41, dpm = se*(critical - qnorm(5e-4)),
                  dpsd = 0, alternative = "greater", lower.tail = FALSE)
expect_equal(do.call(ptbf01, tailpower), 5e-4, tolerance = 1e-7)
expect_equal(suppressWarnings(do.call(ptbf01,
    c(tailpower, list(tail.eps = 1e-3)))), 0)
