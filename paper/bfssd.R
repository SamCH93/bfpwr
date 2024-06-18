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


## ----"plot-power", fig.height = 4.5-------------------------------------------
## BF parameters
k <- 1/10
null <- 0
sd <- sqrt(2) # unit SD for SMD effect size
pm <- dpm <- 0.3 # large SMD
psd1 <- dpsd1 <- 0 # point prior
psd2 <- dpsd2 <- 0.2 # normal prior

nseq <- exp(seq(log(10), log(10^5), length.out = 500))

par(mfrow = c(1, 2), mar = c(5.1, 4.2, 2, 1))
yticks <- seq(0, 100, 20)
xticks <- c(1, 10, 10^2, 10^3, 10^4, 10^5)
lwd <- 1.5
xlabs <- as.expression(c(bquote(1), bquote(10),
                         sapply(seq(2, 5), FUN = function(x) bquote(10^.(x)))))

cols <- adjustcolor(col = c(2, 4), alpha = 0.8)
## point prior under the alternative
plot(nseq, pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd1,
                 dpm = dpm, dpsd = dpsd1)*100,
     xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote({"Pr(BF"["01"] <= 1/.(1/k)} ~ "|" ~ tau == 0 *")"),
     type = "l", xaxt = "n", yaxt = "n", ylim = c(0, 100), col = cols[1],
     lwd = lwd, log = "x", main = "Point analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
axis(side = 1, at = xticks, labels = xlabs)
axis(side = 2, at = yticks, labels = paste0(yticks, "%"), las = 1)
lines(nseq,
      pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd1,
            dpm = dpm, dpsd = dpsd2)*100, lwd = lwd,
      col = cols[2])
zlim <- (null - pm)/(2*dpsd2)
plim <- 1 - pnorm(q = zlim)
abline(h = c(100, plim*100), lty = 2, col = adjustcolor(col = 1, alpha = 0.6))
legend("bottomright", bg = "white", title = "Design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(pm)} * "," ~ tau[italic("d")] == .(psd1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(pm)} * "," ~ tau[italic("d")] == .(psd2) * ")")),
       lwd = lwd, lty = 1, col = cols, cex = 0.7)

## normal prior under the alternative
plot(nseq, pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd2,
                 dpm = dpm, dpsd = dpsd1)*100,
     xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote({"Pr(BF"["01"] <= 1/.(1/k)} ~ "|" ~ tau == .(psd2) *")"),
     type = "l", xaxt = "n", yaxt = "n", ylim = c(0, 100), col = cols[1],
     lwd = lwd, log = "x", main = "Normal analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
axis(side = 1, at = xticks, labels = xlabs)
axis(side = 2, at = yticks, labels = paste0(yticks, "%"), las = 1)
lines(nseq,
      pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd2,
            dpm = dpm, dpsd = dpsd2)*100, lwd = lwd,
      col = cols[2])
abline(h = 100, lty = 2, col = adjustcolor(col = 1, alpha = 0.6))


## ----"asymptotics", eval = FALSE----------------------------------------------
## ## check asymptotics
## library(bfpwr)
## nseq <- exp(seq(log(1), log(10^10), length.out = 500))
## k <- 1/10
## null <- 0
## sd <- sqrt(2)
## dpm <- 0.6
## dpsd1 <- 0
## dpsd2 <- 1.5
## 
## ## BF with point analysis prior
## pm <- 1
## psd <- 0
## pow1 <- pbf01(k = k, n = nseq, null = null, sd = sd, pm = pm, psd = psd, dpm = dpm,
##               dpsd = dpsd1)
## pow2 <- pbf01(k = k, n = nseq, null = null, sd = sd, pm = pm, psd = psd, dpm = dpm,
##               dpsd = dpsd2)
## zlim2 <- (null + pm - 2*dpm)*0.5/dpsd2
## matplot(nseq, cbind(pow1, pow2), type = "l", lty = c(1, 2), col = c(1, 2), lwd = 1.5,
##         xlab = "n", ylab = "Power", las = 1, log = "x", ylim = c(0, 1),
##         main = "Bayes factor with point analysis prior")
## abline(h = 1 - pnorm(zlim2), lty = 3)
## legend("bottomright", c("Point design prior", "Normal design prior"), lty = c(1, 2),
##        lwd = 1.5, col = c(1, 2))
## 
## 
## ## BF with normal alternative
## pm <- 1
## psd <- 2
## pow1 <- pbf01(k = k, n = nseq, null = null, sd = sd, pm = pm, psd = psd, dpm = dpm,
##               dpsd = dpsd1)
## pow2 <- pbf01(k = k, n = nseq, null = null, sd = sd, pm = pm, psd = psd, dpm = dpm,
##               dpsd = dpsd2)
## matplot(nseq, cbind(pow1, pow2), type = "l", lty = c(1, 2), col = c(1, 2), lwd = 1.5,
##         xlab = "n", ylab = "Power", las = 1, log = "x", ylim = c(0, 1),
##         main = "Bayes factor with normal analysis prior")
## abline(h = 1, lty = 3)
## legend("bottomright", c("Point design prior", "Normal design prior"), lty = c(1, 2),
##        lwd = 1.5, col = c(1, 2))


## ----"nTable1", results = "asis"----------------------------------------------
## compute a table with sample sizes
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

powseq <- seq(0.25, 0.999999, length.out = 1000)

par(mfrow = c(1, 2), mar = c(5.1, 4.2, 2, 1))
lwd <- 1.5
xticks <- seq(0, 100, 25)
yticks <- c(1, 10, 10^2, 10^3, 10^4, 10^5)
ylabs <- as.expression(c(bquote(1), bquote(10),
                         sapply(seq(2, 5), FUN = function(x) bquote(10^.(x)))))

cols <- adjustcolor(col = c(2, 4), alpha = 0.8)
## point prior under the alternative
plot(powseq*100,
     nbf01(k = k, power = powseq, sd = sd, null = null, pm = pm, psd = psd1,
           dpm = dpm, dpsd = dpsd1, analytical = TRUE, integer = TRUE),
     xlab = "Target power", ylab = bquote("Sample size per group" ~ italic(n)),
     type = "s", xaxt = "n", yaxt = "n", ylim = c(1, 10^4), col = cols[1],
     lwd = lwd, log = "y", main = "Point analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
axis(side = 1, at = xticks, labels = paste0(xticks, "%"))
axis(side = 2, at = yticks, labels = ylabs, las = 1)
lines(powseq*100,
      nbf01(k = k, power = powseq, sd = sd, null = null, pm = pm, psd = psd1,
           dpm = dpm, dpsd = dpsd2, analytical = TRUE, integer = TRUE),
      col = cols[2], type = "s")
zlim <- (null - pm)/(2*dpsd2)
plim <- 1 - pnorm(q = zlim)
abline(v = c(100, plim*100), lty = 2, col = adjustcolor(col = 1, alpha = 0.6))
legend("bottomleft", bg = "white", title = "Design prior",
       inset = c(0.15, 0),
       legend = c(bquote("N(" * {mu[italic("d")] == .(pm)} * ", " * tau[italic("d")] == .(psd1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(pm)} * ", " * tau[italic("d")] == .(psd2) * ")")),
       lwd = 1.5, lty = 1, col = cols, cex = 0.7)

## normal prior under the alternative
cols2 <- palette.colors(n = 2, palette = "Dark2", alpha = 0.9)
psdlocal1 <- 1
psdlocal2 <- sqrt(2)
nnormal1 <- ceiling(k^2*exp(-lambertWm1(x = -k^2*qnorm(p = powseq/2)^2))*2/psdlocal1^2)
nnormal2 <- ceiling(k^2*exp(-lambertWm1(x = -k^2*qnorm(p = powseq/2)^2))*2/psdlocal2^2)
plot(powseq*100,
     nnormal1,
     ## nbf01(k = k, power = powseq, sd = sd, null = null, pm = null, psd = psdlocal1,
     ##       dpm = null, dpsd = psdlocal1),
     xlab = "Target power", ylab = bquote("Sample size per group" ~ italic(n)),
     type = "s", xaxt = "n", yaxt = "n", ylim = c(1, 10^4), col = cols2[1],
     lwd = lwd, log = "y", main = "Local normal analysis prior",
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(powseq*100, nnormal2, type = "s", col = cols2[2], lwd = lwd)
axis(side = 1, at = xticks, labels = paste0(xticks, "%"))
axis(side = 2, at = yticks, labels = ylabs, las = 1)
abline(v = 100, lty = 2, col = adjustcolor(col = 1, alpha = 0.6))
legend("bottom", bg = "white", title = "Analysis/design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(null)} * ", " * tau[italic("d")] == .(psdlocal1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(null)} * ", " * tau[italic("d")] == sqrt(.(psdlocal2^2)) * ")")),
       lwd = 1.5, lty = 1, col = cols2, cex = 0.7)


## ----"nTable2", results = "asis"----------------------------------------------
## compute a table with sample sizes
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


## ----fig.height = 6-----------------------------------------------------------
## A randomized trial to assess the effectiveness of zanamivir, a new treatment for
## influenza, compared a group randomly allocated to the new treatment with a group
## randomly allocated to a sham (or placebo) treatment. When planning the trial the
## investigators decided that the primary variable would be the number of days to
## the alleviation of symptoms (with alleviation and symptoms defined precisely
## elsewhere in the plan of the study).

## A previous study suggested that a sensible value for sd was 2.75 d and the
## minimal clinically relevant difference that the trial should have good power
## to detect was taken to be 1 d.

## BF parameter
null <- 0
sd <- sqrt(2)*2.75 # unit sd for an MD effect size so that n is the group size
pm <- 1 # the effect on MD scale
psd <- 0 # point analysis prior
k <- 1/10

## design parameters
power <- 0.9
dpm <- pm
dpsd1 <- 0
dpsd2 <- 1/4

## compute required sample size to achieve 80% power for LR
nnum <- nbf01(k = k, power = power, sd = sd, null = null, pm = pm, psd = psd,
              dpm = dpm, dpsd = dpsd1)
zb <- qnorm(p = power)
nanalyt <- ceiling((zb + sqrt(zb^2 - log(k^2)*(pm + null - 2*dpm)/(null - pm)))^2/
                   (pm + null - 2*dpm)^2*sd^2)
## a <- zb^2*dpsd2^2 - ((null + pm)/2 - dpm)^2
## b <- sd^2*(zb^2 - (null + pm - 2*dpm)*log(k)/(null - pm))
## c <- (sd^2*log(k)/(null - pm))^2
a <- ((null + pm)/2 - dpm)^2 - zb^2*dpsd2^2
b <- sd^2*((null + pm - 2*dpm)*log(k)/(null - pm) - zb^2)
c <- (sd^2*log(k)/(null - pm))^2
nanalyt2 <- ceiling((-b + sqrt(b^2 - 4*a*c))/(2*a))

## ## sample size to achieve 80% power under the null
## nbf01(k = 1/k, power = power, sd = sd, null = null, pm = pm, psd = psd,
##       dpm = null, dpsd = 0, lower.tail = FALSE)

## plot power curve
par(mfrow = c(2, 1), mar = c(2.5, 5, 2.5, 2.5))
cols <- rev(palette.colors(n = 4, alpha = 0.95)[2:4])
transpblack <- adjustcolor(col = 1, alpha = 0.2)
nseq <- seq(from = 5, to = 550, by = 1)
pow <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd1)
pow2 <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd2)
powNull <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
                 dpm = null, dpsd = 0)
plot(nseq, pow*100, xlab = "",
     ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * " )"), type = "l",
     ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = cols[1],
     panel.first = graphics::grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(nseq, pow2*100, type = "l", col = cols[2], lwd = 1.5)
lines(nseq, powNull*100, type = "l", col = cols[3], lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
     col = transpblack, cex.axis = 0.8)
axis(side = 3, at = c(nanalyt, nanalyt2), col = transpblack, cex.axis = 0.8)
abline(h = power*100, col = transpblack)
abline(v = c(nanalyt, nanalyt2), col = transpblack)
legend("right", title = "Design prior",
       ## legend = c("Point design prior", "Normal design prior", "Null hypothesis"),
       legend = c(bquote("N(" * {mu[italic("d")] == .(dpm)} * ", " * tau[italic("d")] == .(dpsd1) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpm)} * ", " * tau[italic("d")] == .(dpsd2) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(null)} * ", " * tau[italic("d")] == 0 * ")")),
       lty = 1, lwd = 1.5, col = cols, bg = "white", cex = 0.7)

par(mar = c(4, 5, 1, 2.5))
powH0 <- pbf01(k = 1/k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
               dpm = dpm, dpsd = dpsd1, lower.tail = FALSE)
pow2H0 <- pbf01(k = 1/k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
                dpm = dpm, dpsd = dpsd2, lower.tail = FALSE)
powNullH0 <- pbf01(k = 1/k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
                   dpm = null, dpsd = 0, lower.tail = FALSE)
plot(nseq, powH0*100, xlab = bquote("Sample size per group" ~ italic(n)),
     ylab = bquote("Pr(BF"["01"] > .(1/k) * " )"), type = "l",
     ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = cols[1],
     panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
lines(nseq, pow2H0*100, type = "l", col = cols[2], lwd = 1.5)
lines(nseq, powNullH0*100, type = "l", col = cols[3], lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
abline(v = nanalyt, col = transpblack)
abline(h = power*100, col = transpblack)
axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
     col = transpblack, cex.axis = 0.8)


## ----fig.height = 6-----------------------------------------------------------

## 1) In the example given below, we used two populations with normal
## distributions and a fixed standardized mean difference of delta =0.5 3) In
## the example given below, we analyzed simulated data with a Cauchy prior
## (scale parameter = 1/sqrt(2))


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
n <- nbf01(k = k, power = power, sd = sd, null = null, pm = pm, psd = psd,
           dpm = dpm, dpsd = dpsd)
n2 <- nbf01(k = k, power = power, sd = sd, null = null, pm = pm, psd = psd,
            dpm = dpm, dpsd = dpsd2)
nH0 <- nbf01(k = 1/k, power = power, sd = sd, null = null, pm = pm, psd = psd,
            dpm = null, dpsd = 0, lower.tail = FALSE)

## ## under the null
## nbf01(k = 1/k, power = power, sd = sd, null = null, pm = pm, psd = psd, dpm = null,
##       dpsd = 0, lower.tail = FALSE)

## compute power curve
nseq <- seq(from = 2, to = 800, by = 1)
pow <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd)
powH0 <- pbf01(k = 1/k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd, lower.tail = FALSE)
powNull <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
                 dpm = null, dpsd = 0)
powNullH0 <- pbf01(k = 1/k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
                   dpm = null, dpsd = 0, lower.tail = FALSE)
pow2 <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
              dpm = dpm, dpsd = dpsd2)
pow2H0 <- pbf01(k = 1/k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
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
       ## legend = c("Point design prior", "Normal design prior", "Null hypothesis"),
       legend = c(bquote("N(" * {mu[italic("d")] == .(dpm)} * ", " * tau[italic("d")] == .(dpsd) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpm)} * ", " * tau[italic("d")] == .(dpsd2) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(null)} * ", " * tau[italic("d")] == 0 * ")")),
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
##abline(v = c(n, n2), col = transpblack)


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
## legend("right", title = "Data distribution",
##        legend = c("Point design prior", "Normal design prior", "Null hypothesis"),
##        lty = 1, lwd = 1.5, col = cols, bg = "white", cex = 0.7)


## ----fig.height = 3.5---------------------------------------------------------
dnmoment <- function(x, location = 0, spread = 1) {
    stats::dnorm(x = x, mean = location, sd = spread)*(x - location)^2/spread^2
}
xseq <- seq(-4, 4, length.out = 500)
taus <- c(0.5, 1, 2)
null <- 0
dens <- sapply(X = taus, FUN = function(tau) dnmoment(x = xseq, location = 0, spread = tau))
cols <- hcl.colors(n = length(taus), alpha = 0.8)
par(mar = c(4, 5, 2, 2))
matplot(xseq, dens, type = "l", lty = 1, col = cols, lwd = 1.5, xlab = bquote(theta),
        ylab = "Density", las = 1,
        panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)))
leg <- sapply(taus, function(tau) as.expression(bquote({"NM("* theta[0] == "0,"} * tau == .(tau) * ")")))
legend("topright", legend = leg, col = cols, lty = 1, lwd = 1.5, cex = 0.7,
       bg = "white")


## ----"normal-moment-example", fig.height = 6----------------------------------
nseq <- seq(1, 1100, 1)
dpm <- 0.5
null <- 0
k <- 1/6
sd <- 2
psd <- 0.5/sqrt(2) # mode at 0.5
dpriors <- cbind(dpm = c(0.5, 0.5, 0), dpsd = c(0, 0.1, 0))
power <- 0.95
powH1 <- sapply(X = seq(1, nrow(dpriors)), FUN = function(i) {
    pnmbf01(k = k, n = nseq, sd = sd, null = null, psd = psd, dpm = dpriors[i,1],
            dpsd = dpriors[i,2])
})
powH0 <- sapply(X = seq(1, nrow(dpriors)), FUN = function(i) {
    pnmbf01(k = 1/k, n = nseq, sd = sd, null = null, psd = psd,
            dpm = dpriors[i,1], dpsd = dpriors[i,2], lower.tail = FALSE)
})
nH1 <- ceiling(sapply(X = c(1, 2), FUN = function(i) {
    nnmbf01(k = k, power = power,, sd = sd, null = null, psd = psd,
            dpm = dpriors[i,1], dpsd = dpriors[i,2])
}))
nH0 <- ceiling(sapply(X = 3, FUN = function(i) {
    nnmbf01(k = 1/k, power = power, sd = sd, null = null, psd = psd,
            dpm = dpriors[i,1], dpsd = dpriors[i,2], lower.tail = FALSE)
}))

par(mfrow = c(2, 1), mar = c(3, 5, 2.5, 2.5))
cols <- rev(palette.colors(n = 4, alpha = 0.95)[2:4])
transpblack <- adjustcolor(col = 1, alpha = 0.2)
matplot(nseq, powH1*100, lty = 1, type = "l", las = 1, ylim = c(0, 100),
        panel.first = grid(lty = 3, col = adjustcolor(col = 1, alpha = 0.1)),
        yaxt = "n", xlab = bquote("Sample size per group" ~ italic(n)),
        ylab = bquote("Pr(BF"["01"] < 1/.(1/k) * " )"), col = cols, lwd = 1.5)
axis(side = 2, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"), las = 1)
axis(side = 4, at = power*100, labels = paste0(power*100, "%"), las = 1,
     col = transpblack, cex.axis = 0.8)
axis(side = 3, at = nH1, col = transpblack, cex.axis = 0.8)
abline(v = nH1, col = transpblack)
abline(h = power*100, col = transpblack)
legend("right", bg = "white", title = "Design prior",
       legend = c(bquote("N(" * {mu[italic("d")] == .(dpriors[1,1])} * ", " * tau[italic("d")] == .(dpriors[1,2]) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpriors[2,1])} * ", " * tau[italic("d")] == .(dpriors[2,2]) * ")"),
                  bquote("N(" * {mu[italic("d")] == .(dpriors[3,1])} * ", " * tau[italic("d")] == .(dpriors[3,2]) * ")")),
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
nH0normal <- nbf01(k = 6, power = 0.95, sd = sqrt(2), null = 0, pm = 0,
                   psd = 1/sqrt(2), dpm = 0, dpsd = 0, lower.tail = FALSE)


## ----"package-illustration", echo = TRUE, fig.height = 6----------------------
## install from CRAN or GitHub (the latter requires "remotes" package)
## install.packages("bfpwr") # not yet on CRAN
## remotes::install_github(repo = "SamCH93/bfpwr", subdir = "package")

## load package
library(bfpwr)

## BF parameters
k <- 1/6 # BF threshold
null <- 0 # null value
sd <- 1 # standard deviation of one observation
pm <- null # analysis prior centered around null value
psd <- sqrt(2) # unit information sd for a standardized mean difference
type <- "two.sample" # two-sample test

## design prior
dpm <- 0.5 # design prior mean equal to large SMD effect size
dpsd <- 0.1 # design prior sd to incorporate parameter uncertainty

## determine sample size to achieve 85% power
power <- 0.85
ssd <- powerbf01(k = k, power = power, sd = sd, null = null, pm = pm, psd = psd,
                 dpm = dpm, dpsd = dpsd, type = type)
ssd

## plot power curve
plot(ssd, nlim = c(1, 400))


## ----"simulation-verification", eval = FALSE----------------------------------
## nseq <- seq(from = 1, to = 1000, by = 1)
## pow <- pbf01(k = k, n = nseq, sd = sd, null = null, pm = pm, psd = psd,
##              dpm = dpm, dpsd = dpsd)
## plot(nseq, pow, xlab = "n", ylab = "Power", type = "s", ylim = c(0, 1), las = 1)
## set.seed(123)
## powsim <- sapply(X = nseq, FUN = function(n) {
##     ## verify power with simulation
##     se <- sd/sqrt(n)
##     y <- rnorm(n = 1000, mean = dpm, sd = sqrt(dpsd^2 + se^2))
##     bf <- bf01(estimate = y, se = se, null = null, pm = pm, psd = psd)
##     mean(bf < k)
## })
## lines(nseq, powsim, lty = 2, col = 2)



## ----"sessionInfo2", echo = Reproducibility, results = Reproducibility--------
cat(paste(Sys.time(), Sys.timezone(), "\n"))
sessionInfo()

