library(tinytest)
library(bfpwr)

## When the design prior equals the null, both plot panels must be identical,
## including numerical integration accuracy and any backend controls.
for (controls in list(list(ngrid = 17),
                      list(method = "pmvnorm",
                           algorithm = mvtnorm::GenzBretz(maxpts = 10000,
                                                        abseps = 1e-4)))) {
    z <- do.call(pbf01seq, c(list(k1 = 0.1, k0 = 3,
        se = sqrt(2/c(50, 100, 200)), n = c(50, 100, 200),
        pm = 0, psd = 1, dpm = 0, dpsd = 0), controls))
    plotted <- plot(z, plot = FALSE)
    expect_equal(plotted$pDF1, plotted$pDF2,
                 info = "z null plot uses the design integration settings")

    t <- suppressWarnings(do.call(ptbf01seq, c(list(k1 = 0.1, k0 = 3,
        n = c(50, 100, 200), dpm = 0, dpsd = 0,
        alternative = "greater"), controls)))
    plotted <- suppressWarnings(plot(t, plot = FALSE))
    expect_equal(plotted$pDF1, plotted$pDF2,
                 info = "t null plot uses the design integration settings")
}

## Previously saved designs have no integration metadata and use defaults.
z$integration <- NULL
expect_true(is.list(plot(z, plot = FALSE)))
