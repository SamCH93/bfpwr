## ## do not run these tests for the moment, because they are there to verify
## ## the power with simulation which takes a long time to run

## set.seed(42)
## library(tinytest)
## library(bfpwr)

## ## compare to simulation-based probabilities
## n <- seq(10, 50, 5) # sample size per stage
## se <- sqrt(2/n) # standard errors per stage
## dpm <- 0.3
## dpsd <- 0.1
## pm <- 0.3
## psd <- 0.5
## type <- "normal"
## k0 <- 6
## k1 <- 1/30
## nsim <- 100000
## results <- replicate(n = nsim, expr = {
##     smd <- rnorm(n = 1, mean = dpm, sd = dpsd)
##     y1 <- rnorm(n = max(n), mean = 0, sd = 1)
##     y2 <- rnorm(n = max(n), mean = smd, sd = 1)
##     smd <- sapply(seq_along(n), FUN = function(i) {
##         (mean(y2[1:n[i]]) - mean(y1[1:n[i]]))
##     })
##     se <- sapply(seq_along(n), FUN = function(i) {
##         sqrt(2/n[i])
##     })
##     bf <- sapply(seq_along(n), FUN = function(i) {
##         if (type == "normal") {
##             bf01(estimate = smd[i], se = se[i], null = 0, pm = pm, psd = psd)
##         } else if (type == "directional") {
##             dirbf01(estimate = smd[i], se = se[i], null = 0, pm = pm, psd = psd)
##         } else {
##             nmbf01(estimate = smd[i], se = se[i], null = 0, psd = psd)
##         }
##     })
##     result <- "inconclusive"
##     for (i in seq_along(bf)) {
##         if (bf[i] >= k0) {
##             result <- "H0"
##             break
##         }
##         if (bf[i] <= k1) {
##         ## if ((bf[i] < k1) & t[i] > 0) {
##             result <- "H1"
##             break
##         }
##     }
##     result
## })

## pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm, psd = psd, dpm = dpm,
##          dpsd = dpsd, type = type, strict = TRUE)
## mean(results == "H1")
## mean(results == "H0")
## mean(results == "inconclusive")
