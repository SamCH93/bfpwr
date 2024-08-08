## ## do not run these tests for the moment, because they are there to verify
## ## the power with simulation which takes a long time to run

## library(tinytest)
## library(bfpwr)

## ## verify that computed power equal to simulated power
## grid <- expand.grid(n = c(50, 500, 5000), null = 0, plocation = 0,
##                     pscale = 0.71, pdf = 1, dpm = c(0.3, 0), dpsd = c(0, 0.1),
##                     alternative = c("two.sided", "greater", "less"),
##                     type = c("one.sample", "two.sample"),
##                     stringsAsFactors = FALSE)

## k <- c(1/4, 3)
## set.seed(42)
## nsim <- 1000
## pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
## pow1 <- numeric(nrow(grid))
## pow2 <- numeric(nrow(grid))
## powSim1 <- numeric(nrow(grid))
## powSim2 <- numeric(nrow(grid))
## for (i in seq(1, nrow(grid))) {
##     setTxtProgressBar(pb, i)
##     test <- paste(i, paste(colnames(grid), grid[i,], sep = "=", collapse = ", "))
##     ## simulate
##     if (grid$type[i] == "one.sample") {
##         se <- sqrt(1/grid$n[i])
##     } else {
##         se <- sqrt(2/grid$n[i])
##     }
##     estimatesim <- rnorm(n = nsim, mean = grid$dpm[i], sd = sqrt(grid$dpsd[i]^2 + se^2))
##     tsim <- (estimatesim - grid$null[i])/se
##     tbfsim <- tbf01(t = tsim, n = grid$n[i], plocation = grid$plocation[i],
##                     pscale = grid$pscale[i], pdf = grid$pdf[i],
##                     alternative = grid$alternative[i], type = grid$type[i])
##     if (grid$alternative[i] == "two.sided") {
##         powSim1[i] <- mean(tbfsim <= k[1])
##         powSim2[i] <- mean(tbfsim <= k[2])
##     } else if (grid$alternative[i] == "less") {
##         powSim1[i] <- mean(tbfsim <= k[1] & tsim <= 0)
##         powSim2[i] <- mean(tbfsim <= k[2] & tsim <= 0)
##     } else {
##         powSim1[i] <- mean(tbfsim <= k[1] & tsim >= 0)
##         powSim2[i] <- mean(tbfsim <= k[2] & tsim >= 0)
##     }

##     ## compute numerically
##     pows <- ptbf01(k = k, n = grid$n[i], null = grid$null[i],
##                    plocation = grid$plocation[i], pscale = grid$pscale[i],
##                    pdf = grid$pdf[i], dpm = grid$dpm[i], dpsd = grid$dpsd[i],
##                    alternative = grid$alternative[i], type = grid$type[i],
##                    lower.tail = TRUE, drange = "adaptive")
##     pow1[i] <- pows[1]
##     pow2[i] <- pows[2]
##     cat("\n")
##     print(expect_equal(pow1[i], powSim1[i], info = test, tolerance = 0.04))
##     print(expect_equal(pow2[i], powSim2[i], info = test, tolerance = 0.04))
## }
