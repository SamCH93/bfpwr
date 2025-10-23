## library(tinytest)
## library(bfpwr)

## ## check with simulation that correct probabilities calculated
## simbenchmark <- function(nsim, k1, k0, usd, n, pm, psd, dpm, dpsd, type) {
##     ## simulate stage-wise BFs
##     smd <- rnorm(n = nsim, mean = dpm, sd = dpsd)
##     bfmat <- sapply(X = smd, FUN = function(smdi) {
##         y1 <- rnorm(n = max(n), mean = 0, sd = usd)
##         y2 <- rnorm(n = max(n), mean = smdi*usd, sd = usd)
##         est <- sapply(seq_along(n), FUN = function(i) {
##             (mean(y2[1:n[i]]) - mean(y1[1:n[i]]))/usd
##         })
##         se <- sapply(seq_along(n), FUN = function(i) {
##             sqrt(2/n[i])
##         })
##         bf <- sapply(seq_along(n), FUN = function(i) {
##             if (type == "normal") {
##                 bf01(estimate = est[i], se = se[i], null = 0, pm = pm,
##                      psd = psd)
##             } else if (type == "directional") {
##                 dirbf01(estimate = est[i], se = se[i], null = 0, pm = pm,
##                         psd = psd)
##             } else {
##                 nmbf01(estimate = est[i], se = se[i], null = 0, psd = psd)
##             }
##         })
##         return(bf)
##     })

##     ## estimate probabilities
##     stop <- apply(X = bfmat, MARGIN = 2, FUN = function(x) {
##         result <- "inconclusive"
##         for (xi in x) {
##             if (xi >= k0) {
##                 result <- "H0"
##                 break
##             }
##             if (xi <= k1) {
##                 result <- "H1"
##                 break
##             }
##         }
##         return(result)
##     })
##     pH1sim <- mean(stop == "H1")
##     pH0sim <- mean(stop == "H0")
##     pIncsim <- mean(stop == "inconclusive")

##     ## compute probabilities numerically
##     se <- sqrt(2/n)
##     res <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm,
##                            psd = psd, dpm = dpm, dpsd = dpsd, type = type,
##                            strict = TRUE)

##     ## put everything together
##     out <- data.frame(method = c("simulation", "numerical"),
##                       pH1 = c(pH1sim, tail(res$cumpH1, n = 1)),
##                       pH0 = c(pH0sim, tail(res$cumpH0, n = 1)),
##                       pInc = c(pIncsim, tail(res$cumpInc, n = 1)))
##     return(out)
## }

## simgrid <- expand.grid(k1 = c(1/10, 0.3),
##                        k0 = c(3, 10),
##                        usd = c(sqrt(2)),
##                        pm = c(-0.25, 0.5),
##                        dpm = c(0, 0.5),
##                        psd = c(0.1, 1),
##                        dpsd = c(0, 0.1),
##                        type = c("normal", "moment", "directional"),
##                        stringsAsFactors = FALSE)
## n <- seq(25, 250, 25)
## set.seed(42)
## nsim <- 10000
## simres <- do.call("rbind", lapply(X = seq(1, nrow(simgrid)), FUN = function(i) {
##     res <- simbenchmark(nsim = nsim, k1 = simgrid$k1[i], k0 = simgrid$k0[i],
##                         usd = simgrid$usd[i], n = n, pm = simgrid$pm[i],
##                         psd = simgrid$psd[i], dpm = simgrid$dpm[i],
##                         dpsd = simgrid$dpsd[i], type = simgrid$type[i])
##     condition <- simgrid[i,]
##     data.frame(i = i, condition, res)
## }))

## for (i in seq(1, nrow(simgrid))) {
##     pH1sim <- simres[simres$i == i,]$pH1[1]
##     pH1num <- simres[simres$i == i,]$pH1[2]
##     print(expect_equal(pH1sim, pH1num, info = i, tolerance = 0.01))
## }

## library(ggplot2)
## simplt <- ggplot(data = simres, aes(x = method, y = pH0, color = method,
##                                     group = interaction(pm, dpm, psd, type, k1,
##                                                         k0, dpsd))) +
##     facet_grid(pm + dpm + psd + type ~ k1 + k0 + dpsd,
##                labeller = label_both) +
##     geom_line(color = 1, alpha = 0.8) +
##     geom_point(size = 1) +
##     theme_bw()
## ggsave(plot = simplt, filename = "simbenchmark.pdf", width = 20, height = 20)
