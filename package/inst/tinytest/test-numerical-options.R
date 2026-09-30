library(tinytest)
library(bfpwr)

## Invalid updates must leave every option intact; saved lists can restore
## the effective settings after an update.
original <- bfpwrOptions()
defaults <- bfpwr:::.bfpwr_defaults
expect_equal(original, defaults)
expect_error(bfpwrOptions(ngrid = 20, rel.tol = -1))
expect_equal(bfpwrOptions(), original)
expect_error(bfpwrOptions(ngird = 20))
expect_error(bfpwrOptions(list(ngrid = 10, ngrid = 20)))
expect_error(bfpwrOptions(tail.nquad = 1))
expect_error(bfpwrOptions(tail.eps = 0.5))
expect_error(bfpwrOptions(ngrid = 1.5))

old <- bfpwrOptions(ngrid = 17, tol = 1e-9, rel.tol = 1e-9,
                    subdivisions = 500, tail.eps = 1e-7, tail.nquad = 128)
expect_equal(old, original)
expect_equal(bfpwrOptions()$ngrid, 17)
expect_equal(bfpwr:::.bfpwr_integrate_dots(list()),
             list(rel.tol = 1e-9, abs.tol = 1e-9, subdivisions = 500))
expect_equal(bfpwr:::.bfpwr_integrate_dots(list(abs.tol = 1e-10))$abs.tol, 1e-10)

zargs <- list(k1 = 0.1, n = c(30, 60), se = 1/sqrt(c(30, 60)),
               pm = 0, psd = 1, dpm = 0, dpsd = 0, alternative = "greater")
z <- do.call(pbf01seq, zargs)
expect_equal(z$integration$ngrid, 17L)
expect_equal(do.call(pbf01seq, c(zargs, list(ngrid = 19)))$integration$ngrid, 19L)
t <- ptbf01seq(k1 = 0.1, n = c(30, 60), dpm = 0, dpsd = 0,
                alternative = "greater", type = "one.sample")
expect_equal(t$tail.eps, 1e-7)
expect_equal(t$tail.nquad, 128)
expect_equal(t$integration$bf.control,
    list(tol = 1e-9, rel.tol = 1e-9, abs.tol = 1e-9, subdivisions = 500))
zplot <- plot(z, plot = FALSE)
tplot <- suppressWarnings(plot(t, plot = FALSE))
bfpwrOptions(original)
expect_equal(plot(z, plot = FALSE)$pDF2, zplot$pDF2)
expect_equal(suppressWarnings(plot(t, plot = FALSE))$pDF2, tplot$pDF2)
expect_equal(bfpwrOptions(), original)

## A misspelled or unnamed numerical control must never be silently ignored.
expect_error(ptbf01(k = 0.1, n = 30, rel.toll = 1e-9))
expect_error(bfpwr:::.bfpwr_validate_controls(list(1e-9)))
expect_error(tbf01(t = 1, n = 30, subdivisons = 500))
expect_error(ptbf01seq(n = c(30, 60), bf.control = list(rel.toll = 1e-9)))

## Public signatures expose their option name and literal package fallback.
for (name in c("powerbf01", "powernmbf01", "powerbinbf01")) {
    expect_equal(eval(formals(get(name))$tol), defaults$tol)
}
expect_equal(eval(formals(tbf01)$tail.nquad), defaults$tail.nquad)
expect_equal(eval(formals(ptbf01)$tail.eps), defaults$tail.eps)
