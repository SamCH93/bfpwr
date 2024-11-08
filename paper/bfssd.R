## ----"main-setup", include = FALSE--------------------------------------------
## knitr options
library(knitr)
opts_chunk$set(fig.height = 4,
               echo = FALSE,
               warning = FALSE,
               message = FALSE,
               cache = FALSE,
               eval = TRUE)

## should sessionInfo be printed at the end?
Reproducibility <- TRUE

## packages
library(bfpwr)
library(xtable)
library(lamW)
library(BayesRep)
library(BFDA)
library(dplyr)
library(ggplot2)


## ----"plot-power", fig.height = 4.5-------------------------------------------
## BF parameters
k <- 1/10
null <- 0
sd <- sqrt(2) # unit SD for SMD effect size
pm <- dpm <- 0.3 # large SMD
psd1 <- dpsd1 <- 0 # point prior
psd2 <- dpsd2 <- 0.2 # normal prior

## plot power curves
par(mfrow = c(1, 2), mar = c(5.1, 4.2, 2, 1))
nseq <- exp(seq(log(10), log(10^5), length.out = 500))
yticks <- seq(0, 100, 20)
xticks <- c(1, 10, 10^2, 10^3, 10^4, 10^5)
lwd <- 1.5
xlabs <- as.expression(c(bquote(1), bquote(10),
                         sapply(seq(2, 5), FUN = function(x) bquote(10^.(x)))))
cols <- adjustcolor(col = c(2, 4), alpha = 0.8)
## point prior under the alternative
plot(nseq, pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd1,
                 dpm = dpm, dpsd = dpsd1)*100,
     xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote({"Pr(BF"["01"] <= 1/.(1/k)} ~ "|" ~ tau == 0 *")"),
     type = "l", xaxt = "n", yaxt = "n", ylim = c(0, 100), col = cols[1],
     lwd = lwd, log = "x", main = "Point analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
axis(side = 1, at = xticks, labels = xlabs)
axis(side = 2, at = yticks, labels = paste0(yticks, "%"), las = 1)
lines(nseq,
      pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd1,
            dpm = dpm, dpsd = dpsd2)*100,
      lwd = lwd, col = cols[2])
zlim <- (null - pm)/(2*dpsd2)
plim <- 1 - pnorm(q = zlim)
abline(h = c(100, plim*100), lty = 2, col = adjustcolor(col = 1, alpha = 0.6))
legend("bottomright", bg = "white", title = "Design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(pm)} * ","
                         ~ tau[italic("d")] == .(psd1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(pm)} * ","
                         ~ tau[italic("d")] == .(psd2) * ")")),
       lwd = lwd, lty = 1, col = cols, cex = 0.7)

## normal prior under the alternative
plot(nseq, pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd2,
                 dpm = dpm, dpsd = dpsd1)*100,
     xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote({"Pr(BF"["01"] <= 1/.(1/k)} ~ "|" ~ tau == .(psd2) *")"),
     type = "l", xaxt = "n", yaxt = "n", ylim = c(0, 100), col = cols[1],
     lwd = lwd, log = "x", main = "Normal analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
axis(side = 1, at = xticks, labels = xlabs)
axis(side = 2, at = yticks, labels = paste0(yticks, "%"), las = 1)
lines(nseq,
      pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd2,
            dpm = dpm, dpsd = dpsd2)*100, lwd = lwd,
      col = cols[2])
abline(h = 100, lty = 2, col = adjustcolor(col = 1, alpha = 0.6))


## ----"verify-equations", eval = FALSE-----------------------------------------
## ## ## verify closed-form sample size formulas
## ## null <- 0.1
## ## pm <- 0.4
## ## psd <- 0
## ## dpm1 <- pm
## ## dpm2 <- 0.35
## ## dpsd1 <- 0
## ## dpsd2 <- 0.3
## ## sd <- sqrt(2)
## 
## ## k <- 1/10
## ## power <- 0.63
## ## zb <- qnorm(p = power)
## 
## ## ## formula (9)
## ## nbf01(k = k, power = power, usd = sd, null = null, pm = pm, psd = 0, dpm = pm,
## ##       dpsd = 0, analytical = c(FALSE, TRUE), integer = FALSE)
## ## sd^2*(zb + sqrt(zb^2 - log(k^2)))^2/(pm - null)^2
## 
## ## ## formula (8)
## ## nbf01(k = k, power = power, usd = sd, null = null, pm = pm, psd = 0, dpm = dpm2,
## ##       dpsd = 0, analytical = c(FALSE, TRUE), integer = FALSE)
## ## sd^2*(zb + sqrt(zb^2 - log(k^2)*(null + pm - 2*dpm2)/(null - pm)))^2/(null + pm - 2*dpm2)^2
## 
## ## ## formula (7)
## ## dpsd2 <- 0.1
## ## nbf01(k = k, power = power, usd = sd, null = null, pm = pm, psd = 0, dpm = dpm2,
## ##       dpsd = dpsd2, analytical = c(FALSE, TRUE), integer = FALSE)
## ## A <- sqrt(zb^2 - (2*dpm2 - null - pm)/(pm- null)*log(k^2) + (dpsd2*log(k^2)/(pm - null))^2)
## ## ((zb + A)^2 - (dpsd2*log(k^2)/(pm - null))^2)/(((2*dpm2- pm - null)^2 - 4*zb^2*dpsd2^2)/sd^2)


## ----"nTable1", results = "asis"----------------------------------------------
## produce a table with sample sizes based on closed-form sample size for LRs
kseq <- rev(c(1/1000, 1/300, 1/100, 1/30, 1/10, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3))
powseq <- seq(0.5, 0.95, 0.05)
smd <- 1
results <- sapply(X = kseq, FUN = function(k) {
    beta <- 1 - powseq
    zb <- qnorm(p = 1 - beta)
    (n <- (zb + sqrt(zb^2 - log(k^2)))^2/(smd^2/2))
})
colnames(results) <- formatBF(kseq)
rownames(results) <- paste0(powseq*100, "\\%")
xtab <- xtable(ceiling(results), digits = 0)
addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c(paste0("& \\multicolumn{", length(kseq), "}{c}{$k$} \\\\\n"),
                      paste(paste0("\\cmidrule{2-", length(kseq) + 1, "} \n"),
                            "$1 - \\beta$ &",
                            paste(formatBF(kseq), collapse = " & "), "\\\\"))
print(xtab, booktabs = TRUE, floating = FALSE,
      sanitize.rownames.function = function(x) x, add.to.row = addtorow,
      include.colnames = FALSE)


## ----"plot-n", fig.height = 4.5-----------------------------------------------
## BF parameters
k <- 1/10
null <- 0
sd <- sqrt(2) # unit SD for SMD effect size
pm <- dpm <- 1 # large SMD
psd1 <- dpsd1 <- 0 # point prior
psd2 <- dpsd2 <- 0.5 # normal prior

## plot sample sizes as a function of power
par(mfrow = c(1, 2), mar = c(5.1, 4.2, 2, 1))
powseq <- seq(0.25, 0.999999, length.out = 1000)
lwd <- 1.5
xticks <- seq(0, 100, 25)
yticks <- c(1, 10, 10^2, 10^3, 10^4, 10^5)
ylabs <- as.expression(c(bquote(1), bquote(10),
                         sapply(seq(2, 5), FUN = function(x) bquote(10^.(x)))))
cols <- adjustcolor(col = c(2, 4), alpha = 0.8)
## point prior under the alternative
plot(powseq*100,
     nbf01(k = k, power = powseq, usd = sd, null = null, pm = pm, psd = psd1,
           dpm = dpm, dpsd = dpsd1, analytical = TRUE, integer = TRUE),
     xlab = "Target power", ylab = bquote("Sample size per group" ~ italic(n)),
     type = "s", xaxt = "n", yaxt = "n", ylim = c(1, 10^4), col = cols[1],
     lwd = lwd, log = "y", main = "Point analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
axis(side = 1, at = xticks, labels = paste0(xticks, "%"))
axis(side = 2, at = yticks, labels = ylabs, las = 1)
lines(powseq*100,
      nbf01(k = k, power = powseq, usd = sd, null = null, pm = pm, psd = psd1,
           dpm = dpm, dpsd = dpsd2, analytical = TRUE, integer = TRUE),
      col = cols[2], type = "s")
zlim <- (null - pm)/(2*dpsd2)
plim <- 1 - pnorm(q = zlim)
abline(v = c(100, plim*100), lty = 2, col = adjustcolor(col = 1, alpha = 0.6))
legend("bottomleft", bg = "white", title = "Design prior",
       inset = c(0.15, 0),
       legend = c(bquote("N(" * {mu[italic("d")] == .(pm)} * ", "
                         * tau[italic("d")] == .(psd1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(pm)} * ", "
                         * tau[italic("d")] == .(psd2) * ")")),
       lwd = 1.5, lty = 1, col = cols, cex = 0.7)

## normal prior under the alternative
cols2 <- palette.colors(n = 2, palette = "Dark2", alpha = 0.9)
psdlocal1 <- 1
psdlocal2 <- sqrt(2)
nnormal1 <- ceiling(k^2*exp(-lambertWm1(x = -k^2*qnorm(p = powseq/2)^2))*2/psdlocal1^2)
nnormal2 <- ceiling(k^2*exp(-lambertWm1(x = -k^2*qnorm(p = powseq/2)^2))*2/psdlocal2^2)
plot(powseq*100, nnormal1, xlab = "Target power",
     ylab = bquote("Sample size per group" ~ italic(n)), type = "s", xaxt = "n",
     yaxt = "n", ylim = c(1, 10^4), col = cols2[1], lwd = lwd, log = "y",
     main = "Local normal analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(powseq*100, nnormal2, type = "s", col = cols2[2], lwd = lwd)
axis(side = 1, at = xticks, labels = paste0(xticks, "%"))
axis(side = 2, at = yticks, labels = ylabs, las = 1)
abline(v = 100, lty = 2, col = adjustcolor(col = 1, alpha = 0.6))
legend("bottom", bg = "white", title = "Analysis/design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(null)} * ", "
                         * tau[italic("d")] == .(psdlocal1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(null)} * ", "
                         * tau[italic("d")] == sqrt(.(psdlocal2^2)) * ")")),
       lwd = 1.5, lty = 1, col = cols2, cex = 0.7)


## ----"nTable2", results = "asis"----------------------------------------------
## produce a table with sample sizes computed with formula
kseq <- rev(c(1/1000, 1/300, 1/100, 1/30, 1/10, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3))
powseq <- seq(0.5, 0.95, 0.05)
results <- sapply(X = kseq, FUN = function(k) {
    beta <- 1 - powseq
    k^2*exp(-lambertWm1(x = -k^2*qnorm(p = (1 - beta)/2)^2))
})
colnames(results) <- formatBF(kseq)
rownames(results) <- paste0(powseq*100, "\\%")
xtab <- xtable(ceiling(results), digits = 0)
addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c(paste0("& \\multicolumn{", length(kseq), "}{c}{$k$} \\\\\n"),
                      paste(paste0("\\cmidrule{2-", length(kseq) + 1, "} \n"),
                            "$1 - \\beta$ &",
                            paste(formatBF(kseq), collapse = " & "), "\\\\"))
print(xtab, booktabs = TRUE, floating = FALSE,
      sanitize.rownames.function = function(x) x, add.to.row = addtorow,
      include.colnames = FALSE)


## ----"mirtazapine-example"----------------------------------------------------
## sample size calculation as in Banerjee et al. (2021, p. 1490-1491)
## https://doi.org/10.1016/S0140-6736(21)01210-1
pm <- -6
pow <- 0.8
sd <- 15
alpha <- 0.05
attrition <- 0.1
n <- 2*sd^2*(qnorm(p = 1 - alpha/2) + qnorm(p = pow))^2/pm^2
## power.t.test(sd = sd, power = pow, delta = pm)

## study results
est <- -1.74
ci <- c(-7.17, 3.69)
se <- (ci[2] - ci[1])/(2*qnorm(p = 0.975))
p <- 2*pnorm(q = abs(est/se), lower.tail = FALSE)
lr01 <- bf01(estimate = est, se = se, null = 0, pm = pm, psd = 0)


## ----"mirtazapine-example-design", fig.height = 6-----------------------------
## BF parameter
null <- 0
usd <- sqrt(2)*sd # unit sd for mean difference so that n is the group size
psd <- 0 # point analysis prior
k <- 1/10

## design parameters
power <- 0.8
dpm <- pm
dpsd1 <- 0
dpsd2 <- 2

## compute required sample size to achieve 80% power for LR
nnum <- nbf01(k = k, power = power, usd = usd, null = null, pm = pm, psd = psd,
              dpm = dpm, dpsd = dpsd1)
zb <- qnorm(p = power)
nanalyt <- ceiling((zb + sqrt(zb^2 - log(k^2)*(pm + null - 2*dpm)/(null - pm)))^2/
                   (pm + null - 2*dpm)^2*usd^2)
a <- ((null + pm)/2 - dpm)^2 - zb^2*dpsd2^2
b <- usd^2*((null + pm - 2*dpm)*log(k)/(null - pm) - zb^2)
c <- (usd^2*log(k)/(null - pm))^2
nanalyt2 <- ceiling((-b + sqrt(b^2 - 4*a*c))/(2*a))

## plot power curves
par(mfrow = c(2, 1), mar = c(2.5, 5, 2.5, 2.5))
cols <- rev(palette.colors(n = 4, alpha = 0.95)[2:4])
transpblack <- adjustcolor(col = 1, alpha = 0.2)
nseq <- seq(from = 5, to = 350, by = 1)
pow <- pbf01(k = k, n = nseq, usd = usd, null = null, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd1)
pow2 <- pbf01(k = k, n = nseq, usd = usd, null = null, pm = pm, psd = psd,
              dpm = dpm, dpsd = dpsd2)
powNull <- pbf01(k = k, n = nseq, usd = usd, null = null, pm = pm, psd = psd,
                 dpm = null, dpsd = 0)
plot(nseq, pow*100, xlab = "",
     ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * " )"), type = "l",
     ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = cols[1],
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(nseq, pow2*100, type = "l", col = cols[2], lwd = 1.5)
lines(nseq, powNull*100, type = "l", col = cols[3], lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
axis(side = 4, at = c(power, 0.05)*100,
     labels = paste0(c(power, 0.05)*100, "%"), las = 1, col = transpblack,
     cex.axis = 0.8)
axis(side = 3, at = c(nanalyt, nanalyt2), col = transpblack, cex.axis = 0.8)
abline(h = c(power, 0.05)*100, col = transpblack)
abline(v = c(nanalyt, nanalyt2), col = transpblack)
legend("right", title = "Design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(dpm)} * ", "
                         * tau[italic("d")] == .(dpsd1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpm)} * ", "
                         * tau[italic("d")] == .(dpsd2) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(null)} * ", "
                         * tau[italic("d")] == 0 * ")")),
       lty = 1, lwd = 1.5, col = cols, bg = "white", cex = 0.7)
par(mar = c(4, 5, 1, 2.5))
powH0 <- pbf01(k = 1/k, n = nseq, usd = usd, null = null, pm = pm, psd = psd,
               dpm = dpm, dpsd = dpsd1, lower.tail = FALSE)
pow2H0 <- pbf01(k = 1/k, n = nseq, usd = usd, null = null, pm = pm, psd = psd,
                dpm = dpm, dpsd = dpsd2, lower.tail = FALSE)
powNullH0 <- pbf01(k = 1/k, n = nseq, usd = usd, null = null, pm = pm, psd = psd,
                   dpm = null, dpsd = 0, lower.tail = FALSE)
plot(nseq, powH0*100, xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote("Pr(BF"["01"] > .(1/k) * " )"), type = "l",
     ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = cols[1],
     panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(nseq, pow2H0*100, type = "l", col = cols[2], lwd = 1.5)
lines(nseq, powNullH0*100, type = "l", col = cols[3], lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
abline(v = nanalyt, col = transpblack)
abline(h = c(power, 0.05)*100, col = transpblack)
axis(side = 4, at = c(power, 0.05)*100,
     labels = paste0(c(power, 0.05)*100, "%"), las = 1, col = transpblack,
     cex.axis = 0.8)


## ----"BFDA-comparison", fig.height = 6----------------------------------------
## "1) In the example given below, we used two populations with normal
## distributions and a fixed standardized mean difference of delta = 0.5
## 3) In the example given below, we analyzed simulated data with a Cauchy
## prior" (scale parameter = 1/sqrt(2))"


## ## compare Normal to Cauchy prior
## xseq <- seq(-3, 3, 0.01)
## matplot(xseq, cbind(dnorm(xseq, mean = 0, sd = 1/sqrt(2)),
##                     dcauchy(xseq, location = 0, scale = 1/sqrt(2))),
##         type = "l", lty = c(1, 2), col = 1, xlab = "x", ylab = "density")
## legend("topright", legend = c("Normal(0, variance = 1/2)", "Cauchy(0, scale = 1/sqrt(2))"),
##        lty = c(1, 2))

## BF parameter
null <- 0
sd <- sqrt(2) # unit sd for an SMD effect size so that n is the group size
pm <- 0
psd <- 1/sqrt(2)
k <- 1/6

## design parameters
power <- 0.95
dpm <- 0.5
dpsd <- 0
dpsd2 <- 0.1

## compute required sample size to achieve target power
n <- nbf01(k = k, power = power, usd = sd, null = null, pm = pm, psd = psd,
           dpm = dpm, dpsd = dpsd)
n2 <- nbf01(k = k, power = power, usd = sd, null = null, pm = pm, psd = psd,
            dpm = dpm, dpsd = dpsd2)
nH0 <- nbf01(k = 1/k, power = power, usd = sd, null = null, pm = pm, psd = psd,
             dpm = null, dpsd = 0, lower.tail = FALSE)

## compute power curve
nseq <- seq(from = 2, to = 800, by = 1)
pow <- pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd)
powH0 <- pbf01(k = 1/k, n = nseq, usd = sd, null = null, pm = pm, psd = psd,
               dpm = dpm, dpsd = dpsd, lower.tail = FALSE)
powNull <- pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd,
                 dpm = null, dpsd = 0)
powNullH0 <- pbf01(k = 1/k, n = nseq, usd = sd, null = null, pm = pm, psd = psd,
                   dpm = null, dpsd = 0, lower.tail = FALSE)
pow2 <- pbf01(k = k, n = nseq, usd = sd, null = null, pm = pm, psd = psd,
              dpm = dpm, dpsd = dpsd2)
pow2H0 <- pbf01(k = 1/k, n = nseq, usd = sd, null = null, pm = pm, psd = psd,
                dpm = dpm, dpsd = dpsd2, lower.tail = FALSE)

## compute power curve with BFDA simulation package
nMC <- 1000
## library(BFDA)
## BFDAsim <- BFDA.sim(expected.ES = dpm, type = "t.between",
##                     alternative = "two.sided",
##                     prior = list("normal", list(prior.mean = pm, prior.variance = psd^2)),
##                     design = "sequential", stepsize = 5, n.min = min(nseq),
##                     n.max = max(nseq), B = nMC, seed = 4242, cores = 10)
## ##> Duration: Time difference of 11.64552 mins
## set.seed(43)
## BFDAsim2 <- BFDA.sim(expected.ES = rnorm(n = nMC, mean = dpm, sd = dpsd2),
##                      type = "t.between", alternative = "two.sided",
##                      prior = list("normal", list(prior.mean = pm, prior.variance = psd^2)),
##                      design = "sequential", stepsize = 5, n.min = min(nseq),
##                      n.max = max(nseq), B = nMC, seed = 42462, cores = 10)
## ##> Duration:  Time difference of 11.88571 mins
## BFDAsimNull <- BFDA.sim(expected.ES = null, type = "t.between",
##                         alternative = "two.sided",
##                         prior = list("normal", list(prior.mean = pm, prior.variance = psd^2)),
##                         design = "sequential", stepsize = 5, n.min = min(nseq),
##                         n.max = max(nseq), B = nMC, seed = 4245, cores = 10)
## ##> Duration: Time difference of 8.955649 mins
## save(BFDAsim, BFDAsim2, BFDAsimNull, file = "./BFDAsim.RData")
load(file = "./BFDAsim.RData")
nseqBFDA <- unique(BFDAsim$sim$n)
powBFDA <- lapply(X= list(BFDAsim, BFDAsim2, BFDAsimNull), FUN = function(df) {
    t(sapply(X = nseqBFDA, FUN = function(ni) {
        subdf <- subset(df$sim, n == ni)
        c(mean(subdf$logBF > log(1/k)),
          mean(subdf$logBF < log(k)))
    }))
})

## plot power curves
plotSimCurve <- function(nseq, powe, nMC, col) {
    lines(nseq, powe*100, type = "l", col = col, lwd = 1, lty = 2)
    mcse <- sqrt(powe*(1 - powe)/nMC)
    polygon(x = c(nseq, rev(nseq)), y = c(powe + mcse, rev(powe - mcse))*100,
            border = FALSE, col = adjustcolor(col = col, alpha.f = 0.3),
            lty = 2)
}
par(mfrow = c(2, 1), mar = c(2.5, 5, 2.5, 2.5))
transpblack <- adjustcolor(col = 1, alpha = 0.2)
plot(nseq, pow*100, xlab = "",
     ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * " )"), type = "l",
     ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = cols[1],
     panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(nseq, pow2*100, type = "l", col = cols[2], lwd = 1.5)
lines(nseq, powNull*100, type = "l", col = cols[3], lwd = 1.5)
for (i in seq(1, length(powBFDA))) {
    plotSimCurve(nseq = nseqBFDA, pow = powBFDA[[i]][,1], nMC = nMC,
                 col = cols[i])
}
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
     col = transpblack, cex.axis = 0.8)
axis(side = 3, at = c(n, n2), col = transpblack, cex.axis = 0.8)
abline(h = power*100, col = transpblack)
abline(v = c(n, n2), col = transpblack)
legend("topright", inset = c(0, 1/6), title = "Design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(dpm)} * ", "
                         * tau[italic("d")] == .(dpsd) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpm)} * ", "
                         * tau[italic("d")] == .(dpsd2) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(null)} * ", "
                         * tau[italic("d")] == 0 * ")")),
       lty = 1, lwd = 1.5, col = cols, bg = "white", cex = 0.7)
legend("bottomright", inset = c(0, 1/6), title = "Computational method",
       legend = c("Closed-form", "Simulation"), lty = c(1, 2), lwd = c(1.5, 1),
       col = 1, bg = "white", cex = 0.7)
par(mar = c(4, 5, 1, 2.5))
plot(nseq, powH0*100, xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote("Pr(BF"["01"] > .(1/k) * " )"), type = "l",
     ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = cols[1],
     panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(nseq, pow2H0*100, type = "l", col = cols[2], lwd = 1.5)
lines(nseq, powNullH0*100, type = "l", col = cols[3], lwd = 1.5)
for (i in seq(1, length(powBFDA))) {
    plotSimCurve(nseq = nseqBFDA, pow = powBFDA[[i]][,2], nMC = nMC,
                 col = cols[i])
}
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)


## ----"extensions-t-example", fig.height = 4.5---------------------------------
null <- 0
plocation <- 0
pscale <- 1/sqrt(2)
pdf <- 1
dpm <- 0.5
dpsd1 <- 0
dpsd2 <- 0.1
power <- 0.95
k <- 1/6
alternative <- "greater"
type <- "two.sample"
nex <- ntbf01(k = k, power = power, null = null, plocation = plocation,
              pscale = pscale, pdf = pdf, alternative = alternative,
              type = type, dpm = dpm, dpsd = c(dpsd1, dpsd2))


## library(microbenchmark)
## microbenchmark({
##     ntbf01(k = k, power = power, null = null, plocation = plocation,
##               pscale = pscale, pdf = pdf, alternative = alternative,
##               type = type, dpm = dpm, dpsd = dpsd1)
## })
## #> Unit: milliseconds
## #>       min       lq     mean   median       uq      max neval
## #>  94.27834 96.52802 99.69758 97.41885 98.79241 177.8769   100

## ## plot power curves
## nseq <- seq(5, 300, 1)
## pow <- ptbf01(k = k, n = nseq, null = null, plocation = plocation,
##               pscale = pscale, pdf = pdf, alternative = alternative,
##               type = type, dpm = dpm, dpsd = dpsd1)
## pow2 <- ptbf01(k = k, n = nseq, null = null, plocation = plocation,
##                pscale = pscale, pdf = pdf, alternative = alternative,
##                type = type, dpm = dpm, dpsd = dpsd2)
## powNull <- ptbf01(k = k, n = nseq, null = null, plocation = plocation,
##                   pscale = pscale, pdf = pdf, alternative = alternative,
##                   type = type, dpm = null, dpsd = 0)
## par(mar = c(4, 5, 4, 2.5))
## matplot(nseq, cbind(pow, pow2, powNull)*100, type = "s", ylim = c(0, 100),
##         lty = 1, col = cols, ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * " )"),
##         xlab = bquote("Sample size per group" ~ italic(n)), lwd = 1.5,
##         yaxt = "n",
##         panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
## axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
## abline(v = nex, col = transpblack)
## axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
##      col = transpblack, cex.axis = 0.8)
## axis(side = 3, at = nex, col = transpblack, cex.axis = 0.8)
## abline(h = power*100, col = transpblack)
## legend("right",
##        legend = c("Point design prior", "Normal design prior", "Null hypothesis"),
##        lty = 1, lwd = 1.5, col = cols, bg = "white", cex = 0.7)


## -----------------------------------------------------------------------------
## microbenchmark::microbenchmark({
##     powertbf01(k = 1/6, power = 0.95, null = 0, plocation = 0, pscale = 1/sqrt(2),
##                pdf = 1, alternative = "greater", dpm = 0.5, dpsd = 0)
## })
## ##> + Unit: milliseconds
## ##>       min       lq     mean   median       uq      max neval
## ##>  83.67747 85.16477 92.08307 88.00866 98.46287 117.2702   100


## ----"moment-priors-illustration-plot", fig.height = 3.25---------------------
dnmoment <- function(x, location = 0, spread = 1) {
    stats::dnorm(x = x, mean = location, sd = spread)*(x - location)^2/spread^2
}
xseq <- seq(-4, 4, length.out = 500)
taus <- c(0.5, 1, 2)
null <- 0
dens <- sapply(X = taus, FUN = function(tau) dnmoment(x = xseq, location = 0, spread = tau))
cols <- hcl.colors(n = length(taus), alpha = 0.8)
par(mar = c(4, 5, 1, 1))
matplot(xseq, dens, type = "l", lty = 1, col = cols, lwd = 1.5,
        xlab = bquote("Parameter" ~ theta), ylab = "Density", las = 1,
        panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
leg <- sapply(taus, function(tau) as.expression(bquote({"NM("* theta[0] == "0,"} * tau == .(tau) * ")")))
legend("topright", legend = leg, col = cols, lty = 1, lwd = 1.5, cex = 0.9,
       bg = "white")


## ----"normal-moment-example", fig.height = 6----------------------------------
## plot power curves for normal moment prior example
nseq <- seq(1, 1100, 1)
dpm <- 0.5
null <- 0
k <- 1/6
sd <- 2
psd <- 0.5/sqrt(2) # mode at 0.5
dpriors <- cbind(dpm = c(0.5, 0.5, 0), dpsd = c(0, 0.1, 0))
power <- 0.95
powH1 <- sapply(X = seq(1, nrow(dpriors)), FUN = function(i) {
    pnmbf01(k = k, n = nseq, usd = sd, null = null, psd = psd, dpm = dpriors[i,1],
            dpsd = dpriors[i,2])
})
powH0 <- sapply(X = seq(1, nrow(dpriors)), FUN = function(i) {
    pnmbf01(k = 1/k, n = nseq, usd = sd, null = null, psd = psd,
            dpm = dpriors[i,1], dpsd = dpriors[i,2], lower.tail = FALSE)
})
nH1 <- ceiling(sapply(X = c(1, 2), FUN = function(i) {
    nnmbf01(k = k, power = power,, usd = sd, null = null, psd = psd,
            dpm = dpriors[i,1], dpsd = dpriors[i,2])
}))
nH0 <- ceiling(sapply(X = 3, FUN = function(i) {
    nnmbf01(k = 1/k, power = power, usd = sd, null = null, psd = psd,
            dpm = dpriors[i,1], dpsd = dpriors[i,2], lower.tail = FALSE)
}))

par(mfrow = c(2, 1), mar = c(3, 5, 2.5, 2.5))
cols <- rev(palette.colors(n = 4, alpha = 0.95)[2:4])
transpblack <- adjustcolor(col = 1, alpha = 0.2)
matplot(nseq, powH1*100, lty = 1, type = "l", las = 1, ylim = c(0, 100),
        panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)),
        yaxt = "n", xlab = "",
        ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * " )"), col = cols, lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
     col = transpblack, cex.axis = 0.8)
axis(side = 3, at = nH1, col = transpblack, cex.axis = 0.8)
abline(v = nH1, col = transpblack)
abline(h = power*100, col = transpblack)
legend("right", bg = "white", title = "Design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(dpriors[1,1])} * ", "
                         * tau[italic("d")] == .(dpriors[1,2]) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpriors[2,1])} * ", "
                         * tau[italic("d")] == .(dpriors[2,2]) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpriors[3,1])} * ", "
                         * tau[italic("d")] == .(dpriors[3,2]) * ")")),
       lty = 1, lwd = 1.5, col = cols, cex = 0.7)
par(mar = c(4, 5, 1.5, 2.5))
matplot(nseq, powH0*100, lty = 1, type = "l", las = 1, ylim = c(0, 100),
        panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)),
        yaxt = "n", xlab = bquote("Sample size per group" ~ italic(n)),
        ylab = bquote("Pr(BF"["01"] > .(1/k) * " )"), col = cols, lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
     col = transpblack, cex.axis = 0.8)
axis(side = 3, at = nH0, col = transpblack, cex.axis = 0.8)
abline(v = nH0, col = transpblack)
abline(h = power*100, col = transpblack)

## nH0 for normal analysis prior for comparison
nH0normal <- nbf01(k = 6, power = 0.95, usd = sqrt(2), null = 0, pm = 0,
                   psd = 1/sqrt(2), dpm = 0, dpsd = 0, lower.tail = FALSE)


## ----"package-illustration", echo = TRUE, fig.height = 6----------------------
library(bfpwr)

## BF parameters
k <- 1/6 # BF threshold
null <- 0 # null value is zero
sd <- 1 # standard deviation of one observation
pm <- 0 # analysis prior mean set to zero
psd <- sqrt(2) # analysis prior SD set to sqrt(2)
type <- "two.sample" # two-sample test

## design prior
dpm <- 0.5 # design prior mean equal to medium SMD effect size
dpsd <- 0.1 # positive design prior SD to incorporate parameter uncertainty

## determine sample size to achieve 85% power
power <- 0.85
ssd <- powerbf01(k = k, power = power, sd = sd, null = null, pm = pm, psd = psd,
                 dpm = dpm, dpsd = dpsd, type = type)
ssd

## plot power curve
plot(ssd, nlim = c(1, 400))


## ----child = "appendix.Rnw"---------------------------------------------------

## ----"sim-evaluation-params"--------------------------------------------------
pow <- 0.8
k <- 1/10
usd <- sqrt(2)
null <- 0
nsim <- 50000
pargrid <- expand.grid(k = k,
                       power = pow,
                       usd = usd,
                       null = null,
                       pm = c(0, 0.2, 0.5, 0.8),
                       psd = c(0, 0.5, 1, 2),
                       dpm = c(0, 0.2, 0.5, 0.8),
                       dpsd = c(0, 0.5, 1))


## ----"additional-simulations1", cache = TRUE----------------------------------
library(bfpwr)

## function to simulate data and see whether closed-form / root-finding solution
## aligns with Monte Carlo power
simbenchmark <- function(nsim = 10000, k = 1/10, power = 0.8, usd, null = 0, pm,
                         psd, dpm, dpsd, lower.tail = TRUE, analytical = TRUE,
                         integer = FALSE) {

    ## compute n with closed-form solution
    n <- try(nbf01(k = k, power = power, usd = usd, null = null, pm = pm,
                   psd = psd, dpm = dpm, dpsd = dpsd, lower.tail = lower.tail,
                   analytical = analytical, integer = integer))
    if (is.nan(n)) {
        power_closed <- NaN
        power_MC <- NaN
    } else {

        ## simulate parameter estimates
        se <- usd/sqrt(n)
        est <- rnorm(n = nsim, mean = dpm, sd = sqrt(se^2 + dpsd^2))

        ## compute Bayes factors
        bf <- bf01(estimate = est, se = se, null = null, pm = pm, psd = psd)

        ## recompute power
        power_closed <- pbf01(k = k, n = n, usd = usd, null = null, pm = pm,
                              psd = psd, dpm = dpm, dpsd = dpsd,
                              lower.tail = lower.tail)
        if (lower.tail == TRUE) {
            power_MC <- mean(bf <= k)

        } else {
            power_MC <- mean(bf > k)
        }
    }
    res <- data.frame(nsim, k, usd, null, pm, psd, dpm, dpsd, lower.tail, n,
                      power, power_closed, power_MC)
    return(res)
}

## simulation benchmarking
set.seed(424242)
res <- do.call("rbind", lapply(X = seq(1, nrow(pargrid)), FUN = function(i) {
    simbenchmark(nsim = nsim, k = pargrid$k[i], power = pargrid$power[i],
                 usd = pargrid$usd[i], null = pargrid$null[i],
                 pm = pargrid$pm[i], psd = pargrid$psd[i], dpm = pargrid$dpm[i],
                 dpsd = pargrid$dpsd[i])
}))
res$mcse <- sqrt(res$power_MC*(1 - res$power_MC)/nsim)

## ----"additional-simulations2", fig.height = 8, fig.width = 9-----------------
library(dplyr)
library(ggplot2)
nDF <- res |>
    group_by(k, power, usd, null, pm, psd, dpm, dpsd) |>
    summarise(n = ceiling(unique(n))) |>
    ungroup() |>
    mutate(nFormat = ifelse(is.nan(n), "x", n))

ggplot(data = res,
       aes(x = factor(dpm, ordered = TRUE), y = power_MC,
           color = factor(dpsd, ordered = TRUE))) +
    geom_vline(xintercept = seq(0.5, 10.5, 1), alpha = 0.1, lty = 3) +
    geom_hline(yintercept = pow, lty = 2, alpha = 0.7) +
    facet_grid(psd ~pm,
               labeller = label_bquote(cols = "Analysis prior mean" ~ mu == .(pm),
                                       rows = "Analysis prior SD" ~ tau == .(psd))) +
    geom_pointrange(aes(ymin = power_MC - mcse, ymax = power_MC + mcse),
                    position = position_dodge2(width = 0.5),
                    fatten = 1.5) +
    geom_text(data = nDF, aes(y = pow + 0.0075, label = nFormat), size = 2.5,
              position = position_dodge2(width = 1), angle = 30,
              fontface = "bold") +
    labs(x = bquote("Design prior mean" ~ mu[italic(d)]),
         y = bquote("Estimated power" %+-% "MCSE"),
         color =  bquote("Design prior SD" ~ tau[italic(d)])) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = pow + c(-1, 1)*0.008) +
    scale_color_viridis_d(end = 0.7, option = "A") +
    theme_bw() +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = alpha("black", 0.01)))




## ----"sessionInfo2", echo = Reproducibility, results = Reproducibility--------
cat(paste(Sys.time(), Sys.timezone(), "\n"))
sessionInfo()

