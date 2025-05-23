## do not run these tests for the moment, because they are there to verify
## the power with simulation which takes a long time to run

## library(tinytest)
## library(bfpwr)

## ## verify that computed power equal to simulated power
## grid1 <- expand.grid(n = c(10, 100, 1000),
##                     p0 = 0.4, type = "point",
##                     a = c(1, 2), b = c(1, 2),
##                     dp = c(NA, 0.4, 0.5),
##                     da = 2, db = 3,
##                     dl = 0, du = 1,
##                     stringsAsFactors = FALSE)
## grid2 <- expand.grid(n = c(10, 100, 1000),
##                     p0 = 0.4, type = "direction",
##                     a = c(1, 2), b = c(1, 2),
##                     dp = c(NA, 0.5),
##                     da = 2, db = 3,
##                     dl = 0.2, du = 0.8,
##                     stringsAsFactors = FALSE)
## grid <- rbind(grid1, grid2)
## k <- c(1/4, 3)
## set.seed(42)
## nsim <- 100000
## pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
## pow1 <- numeric(nrow(grid))
## pow2 <- numeric(nrow(grid))
## powSim1 <- numeric(nrow(grid))
## powSim2 <- numeric(nrow(grid))
## for (i in seq(1, nrow(grid))) {
##     setTxtProgressBar(pb, i)
##     test <- paste(i, paste(colnames(grid), grid[i,], sep = "=", collapse = ", "))
##     ## simulate
##     if (!is.na(grid$dp[i])) {
##         p <- grid$dp[i]
##     } else {
##         if (grid$type[i] == "point") {
##             p <- rbeta(n = nsim, shape1 = grid$da[i], shape2 = grid$db[i])
##         } else {
##             pprop <- rbeta(n = 100*nsim, shape1 = grid$da[i], shape2 = grid$db[i])
##             p <- pprop[(pprop >= grid$dl[i]) & (pprop <= grid$du[i])][1:nsim]
##         }
##     }
##     x <- rbinom(n = nsim, size = grid$n[i], prob = p)
##     bfsim <- binbf01(x = x, n = grid$n[i], p0 = grid$p0[i], type = grid$type[i],
##                      a = grid$a[i], b = grid$b[i])
##     powSim1[i] <- mean(bfsim <= k[1])
##     powSim2[i] <- mean(bfsim <= k[2])

##     ## compute numerically
##     pows <- pbinbf01(k = k, n = grid$n[i], p0 = grid$p0[i], type = grid$type[i],
##                      a = grid$a[i], b = grid$b[i], dp = grid$dp[i],
##                      da = grid$da[i], db = grid$db[i],
##                      dl = grid$dl[i], du = grid$du[i])
##     pow1[i] <- pows[1]
##     pow2[i] <- pows[2]
##     cat("\n")
##     print(expect_equal(pow1[i], powSim1[i], info = test, tolerance = 0.01))
##     print(expect_equal(pow2[i], powSim2[i], info = test, tolerance = 0.01))
## }

## round(cbind(powSim1, pow1, powSim1 - pow1,
##             powSim2, pow2, powSim2 - pow2), 4)

#>       powSim1   pow1         powSim2   pow2
#>  [1,]  0.1624 0.1608  0.0015  1.0000 1.0000  0.0000
#>  [2,]  0.5774 0.5762  0.0011  0.7648 0.7651 -0.0003
#>  [3,]  0.8384 0.8387 -0.0003  0.8883 0.8882  0.0001
#>  [4,]  0.0944 0.0949 -0.0005  0.7159 0.7163 -0.0004
#>  [5,]  0.5737 0.5732  0.0005  0.7313 0.7322 -0.0009
#>  [6,]  0.8343 0.8352 -0.0009  0.8829 0.8830 -0.0001
#>  [7,]  0.1605 0.1608 -0.0003  1.0000 1.0000  0.0000
#>  [8,]  0.5763 0.5762  0.0000  0.7825 0.7824  0.0001
#>  [9,]  0.8438 0.8421  0.0017  0.8950 0.8934  0.0016
#> [10,]  0.1568 0.1608 -0.0040  1.0000 1.0000  0.0000
#> [11,]  0.5904 0.5901  0.0004  0.8145 0.8155 -0.0010
#> [12,]  0.8438 0.8437  0.0001  0.8981 0.8985 -0.0004
#> [13,]  0.0186 0.0183  0.0003  1.0000 1.0000  0.0000
#> [14,]  0.0080 0.0078  0.0002  0.1570 0.1545  0.0025
#> [15,]  0.0024 0.0024  0.0000  0.0371 0.0359  0.0013
#> [16,]  0.0122 0.0123 -0.0001  0.5359 0.5342  0.0017
#> [17,]  0.0082 0.0082  0.0000  0.1046 0.1036  0.0009
#> [18,]  0.0022 0.0020  0.0002  0.0291 0.0282  0.0009
#> [19,]  0.0187 0.0183  0.0003  1.0000 1.0000  0.0000
#> [20,]  0.0080 0.0078  0.0002  0.1839 0.1842 -0.0004
#> [21,]  0.0029 0.0030 -0.0001  0.0464 0.0454  0.0010
#> [22,]  0.0184 0.0183  0.0000  1.0000 1.0000  0.0000
#> [23,]  0.0103 0.0104  0.0000  0.2637 0.2614  0.0022
#> [24,]  0.0036 0.0033  0.0003  0.0575 0.0568  0.0007
#> [25,]  0.0563 0.0557  0.0006  1.0000 1.0000  0.0000
#> [26,]  0.2424 0.2421  0.0003  0.7595 0.7581  0.0013
#> [27,]  0.9996 0.9996  0.0001  1.0000 1.0000  0.0000
#> [28,]  0.0544 0.0547 -0.0003  0.6781 0.6777  0.0004
#> [29,]  0.3093 0.3087  0.0007  0.6908 0.6914 -0.0007
#> [30,]  0.9996 0.9996  0.0001  1.0000 1.0000  0.0000
#> [31,]  0.0555 0.0557 -0.0002  1.0000 1.0000  0.0000
#> [32,]  0.2419 0.2421 -0.0001  0.7586 0.7584  0.0002
#> [33,]  0.9996 0.9996 -0.0001  1.0000 1.0000  0.0000
#> [34,]  0.0559 0.0557  0.0003  1.0000 1.0000  0.0000
#> [35,]  0.3099 0.3087  0.0012  0.8194 0.8168  0.0026
#> [36,]  0.9997 0.9997  0.0000  1.0000 1.0000  0.0000
#> [37,]  0.3194 0.3189  0.0005  0.6423 0.6444 -0.0021
#> [38,]  0.4508 0.4518 -0.0010  0.6177 0.6193 -0.0017
#> [39,]  0.5304 0.5300  0.0004  0.5809 0.5798  0.0010
#> [40,]  0.1873 0.1883 -0.0011  0.6460 0.6444  0.0016
#> [41,]  0.4114 0.4123 -0.0009  0.5552 0.5550  0.0002
#> [42,]  0.5112 0.5108  0.0004  0.5559 0.5559  0.0000
#> [43,]  0.3200 0.3189  0.0010  0.8018 0.8025 -0.0007
#> [44,]  0.5116 0.5130 -0.0014  0.6610 0.6628 -0.0019
#> [45,]  0.5488 0.5472  0.0016  0.6018 0.5996  0.0022
#> [46,]  0.3193 0.3189  0.0003  0.6466 0.6444  0.0022
#> [47,]  0.4516 0.4518 -0.0002  0.5985 0.5978  0.0007
#> [48,]  0.5252 0.5279 -0.0027  0.5751 0.5777 -0.0026
#> [49,]  0.3773 0.3770  0.0003  0.8295 0.8281  0.0014
#> [50,]  0.8161 0.8159  0.0002  0.9940 0.9940  0.0001
#> [51,]  1.0000 1.0000  0.0000  1.0000 1.0000  0.0000
#> [52,]  0.1742 0.1719  0.0023  0.8272 0.8281 -0.0009
#> [53,]  0.6877 0.6914 -0.0037  0.9714 0.9716 -0.0002
#> [54,]  1.0000 1.0000  0.0000  1.0000 1.0000  0.0000
#> [55,]  0.3818 0.3770  0.0048  0.9452 0.9453 -0.0001
#> [56,]  0.9328 0.9334 -0.0006  0.9984 0.9982  0.0001
#> [57,]  1.0000 1.0000  0.0000  1.0000 1.0000  0.0000
#> [58,]  0.3785 0.3770  0.0016  0.8272 0.8281 -0.0009
#> [59,]  0.8161 0.8159  0.0002  0.9895 0.9895  0.0000
#> [60,]  1.0000 1.0000  0.0000  1.0000 1.0000  0.0000
