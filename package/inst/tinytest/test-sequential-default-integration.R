library(tinytest)
library(bfpwr)

## The public design and search interfaces share the same integration default,
## and an explicit grid size must still reach the returned design.
zargs <- list(k1 = 0.1, pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
              alternative = "greater")
targs <- list(k1 = 0.1, dpm = 0.5, dpsd = 0, alternative = "greater")
calls <- list(
    list(fun = pbf01seq, args = c(zargs, list(n = c(20, 40, 60),
                                             se = 1/sqrt(c(20, 40, 60))))),
    list(fun = ptbf01seq, args = c(targs, list(n = c(20, 40, 60)))),
    list(fun = nbf01seq, args = c(zargs, list(power = 0.1,
        nrange = c(20, 60), minN = 20, by = 20, details = TRUE))),
    list(fun = ntbf01seq, args = c(targs, list(power = 0.1,
        nrange = c(20, 60), minN = 20, by = 20, details = TRUE))),
    list(fun = powerbf01seq, args = c(zargs, list(n = 60, looks = 3))),
    list(fun = powertbf01seq, args = c(targs, list(n = 60, looks = 3)))
)
for (call in calls) {
    for (grid in list(NULL, 10000L, 17L)) {
        args <- call$args
        if (!is.null(grid)) args$ngrid <- grid
        result <- do.call(call$fun, args)
        design <- if ("result" %in% names(result)) result$result else result
        expect_equal(design$integration,
            list(method = "lpmvnorm", ngrid = if (is.null(grid)) 10000L else grid))
        if (is.null(grid)) default <- design
        if (identical(grid, 10000L)) {
            expect_equal(design$cumpH1, default$cumpH1)
            expect_equal(design$cumpH0, default$cumpH0)
        }
    }
}
