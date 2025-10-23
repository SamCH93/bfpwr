## ## do not run these tests for the moment, because they are there to verify
## ## the power with simulation which takes a long time to run

## library(tinytest)
## library(bfpwr)

## ## verify that computed power equal to simulated power
## grid <- expand.grid(plocation = 0,
##                     pscale = 0.71,
##                     pdf = 1,
##                     dpm = c(0.5, 0),
##                     dpsd = c(0, 0.1),
##                     alternative = c("two.sided", "greater", "less"),
##                     type = c("one.sample", "two.sample"),
##                     stringsAsFactors = FALSE)
## n <- seq(40, 100, 10)
## k0 <- 6
## k1 <- 1/30
## set.seed(42)
## nsim <- 10000
## pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
## pH1sim <- numeric(nrow(grid))
## pH0sim <- numeric(nrow(grid))
## pIncsim <- numeric(nrow(grid))
## pH1 <- numeric(nrow(grid))
## pH0 <- numeric(nrow(grid))
## pInc <- numeric(nrow(grid))
## for (i in seq(1, nrow(grid))) {
##     setTxtProgressBar(pb, i)
##     test <- paste(i, paste(colnames(grid), grid[i,], sep = "=", collapse = ", "))
##     ## simulate
##     results <- replicate(n = nsim, expr = {
##         smd <- rnorm(n = 1, mean = grid$dpm[i], sd = grid$dpsd[i])
##         if (grid$type[i] == "one.sample") {
##             y <- rnorm(n = max(n), mean = smd, sd = 1)
##             t <- sapply(seq_along(n), FUN = function(j) {
##                 ttest <- t.test(y[1:n[j]], alternative = "two.sided")$statistic
##             })
##         } else {
##             y1 <- rnorm(n = max(n), mean = 0, sd = 1)
##             y2 <- rnorm(n = max(n), mean = smd, sd = 1)
##             t <- sapply(seq_along(n), FUN = function(j) {
##                 ttest <- t.test(y2[1:n[j]], y1[1:n[j]], var.equal = TRUE,
##                                 alternative = "two.sided")$statistic
##             })
##         }
##         bf <- sapply(seq_along(n), FUN = function(j) {
##             tbf01(t = t[j], n = n[j], plocation = grid$plocation[i],
##                   pscale = grid$pscale[i], pdf = grid$pdf[i],
##                   type = grid$type[i], alternative = grid$alternative[i])
##         })
##         result <- "inconclusive"
##         for (j in seq_along(bf)) {
##             if (bf[j] >= k0) {
##                 result <- "H0"
##                 break
##             }
##             if (bf[j] <= k1) {
##                 result <- "H1"
##                 break
##             }
##         }
##         result
##     })
##     pH1sim[i] <- mean(results == "H1")
##     pH0sim[i] <- mean(results == "H0")
##     pIncsim[i] <- mean(results == "inconclusive")
##     res <- ptbf01seq(k1 = k1, k0 = k0, n1 = n, n2 = n,
##                      plocation = grid$plocation[i], pscale = grid$pscale[i],
##                      pdf = grid$pdf[i], type = grid$type[i], dpm = grid$dpm[i],
##                      dpsd = grid$dpsd[i], alternative = grid$alternative[i],
##                      strict = FALSE)
##     pH1[i] <- tail(res$cumpH1, n = 1)
##     pH0[i] <- tail(res$cumpH0, n = 1)
##     pInc[i] <- tail(res$cumpInc, n = 1)
##     cat("\n")
##     print(expect_equal(pH1[i], pH1sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pH0[i], pH0sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pInc[i], pIncsim[i], info = test, tolerance = 0.02))
## }


## for (i in seq(1, nrow(grid))) {
##     cat("\n")
##     print(expect_equal(pH1[i], pH1sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pH0[i], pH0sim[i], info = test, tolerance = 0.02))
##     print(expect_equal(pInc[i], pIncsim[i], info = test, tolerance = 0.02))
## }

## ## percentage scale
## round(cbind(pH1, pH1sim, pH0, pH0sim, pInc, pIncsim)*100, 2)

## ## percentage differences
## round(cbind(pH1 - pH1sim, pH0 - pH0sim, pInc - pIncsim)*100, 2)
