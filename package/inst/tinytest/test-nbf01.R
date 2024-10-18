library(tinytest)
library(bfpwr)

## verify that computed sample size leads to desired power
grid <- expand.grid(k = c(1/5), usd = c(0.5, 1.5), null = c(0, 0.1),
                    pm = c(1, 2), psd = c(0, 1), dpm = c(1, 1.5),
                    dpsd = c(1, 2), power = c(0.6, 0.7),
                    analytical = c(FALSE, TRUE))
nupper <- 10^5
nlower <- 1
for (i in seq(1, nrow(grid))) {
    test <- paste(i, paste(colnames(grid), grid[i,], sep = "=", collapse = ", "))
    n <- suppressWarnings(nbf01(k = grid$k[i], power = grid$power[i],
                                usd = grid$usd[i], null = grid$null[i],
                                pm = grid$pm[i], psd = grid$psd[i],
                                dpm = grid$dpm[i], dpsd = grid$dpsd[i],
                                integer = FALSE))
    if (is.nan(n)) {
        pupper <- pbf01(k = grid$k[i], n = nupper, usd = grid$usd[i],
                        null = grid$null[i], pm = grid$pm[i], psd = grid$psd[i],
                        dpm = grid$dpm[i], dpsd = grid$dpsd[i])
        plower <- pbf01(k = grid$k[i], n = nlower, usd = grid$usd[i],
                        null = grid$null[i], pm = grid$pm[i], psd = grid$psd[i],
                        dpm = grid$dpm[i], dpsd = grid$dpsd[i])
        expect_true(pupper < grid$power[i] | plower > grid$power[i], info = test)
    } else {
        p <- pbf01(k = grid$k[i], n = n, usd = grid$usd[i], null = grid$null[i],
                   pm = grid$pm[i], psd = grid$psd[i], dpm = grid$dpm[i],
                   dpsd = grid$dpsd[i])
        expect_equal(p, grid$power[i], info = test, tolerance = 0.0001)
    }
}
