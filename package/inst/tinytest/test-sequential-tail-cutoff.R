library(tinytest)
library(bfpwr)

## A cutoff can lie beyond either an almost empty or an almost certain H0
## event. Compare to explicitly bracketed t boundaries in both directions.
for (alternative in c("greater", "less")) {
    dpm <- if (alternative == "greater") -1 else 1
    args <- list(k1 = 0.1, k0 = 3, n = c(500, 1000), dpm = dpm,
                 dpsd = 0, alternative = alternative)
    approximate <- suppressWarnings(do.call(ptbf01seq, args))
    reference <- do.call(ptbf01seq, c(args, list(trange = c(-10, 10))))
    expect_equal(approximate$cumpH0, reference$cumpH0, tolerance = 0.004,
                 info = "tail cutoffs preserve almost certain H0 stopping")
    expect_equal(approximate$cumpH1, reference$cumpH1, tolerance = 0.004)
    expect_equal(approximate$EN1, 500,
                 info = "almost certain H0 stopping occurs at the first look")

    ## The sample-size evaluator uses the same approximate boundaries as the
    ## direct probability API, including when both cutoffs equal zero.
    found <- suppressWarnings(ntbf01seq(
        k1 = 0.1, k0 = 3, power = 0.9, target = "H0", dpm = dpm,
        dpsd = 0, alternative = alternative, nrange = c(500, 1000),
        minN = 500, by = 500, details = TRUE
    ))
    expect_true(found$reached)
    expect_equal(found$n, 500)
    expect_equal(found$result$cumpH0, approximate$cumpH0[1])

    local({
        grDevices::pdf(NULL)
        on.exit(grDevices::dev.off())
        plotWarnings <- character()
        plotted <- withCallingHandlers(plot(approximate, zplot = TRUE),
            warning = function(w) {
                plotWarnings <<- c(plotWarnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        expect_equal(plotWarnings, character())
        expect_true(is.list(plotted),
                    info = "unbounded cutoff approximations remain plottable")
    })
}
