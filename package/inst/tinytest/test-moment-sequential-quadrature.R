library(tinytest)
library(bfpwr)

## Conditional-normal integration independently checks two-look regions.
n <- c(15, 40)
usd <- sqrt(2)
se <- usd/sqrt(n)
null <- 2
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    dpm <- null + direction*0.4
    dpsd <- 0.3
    design <- pbf01seq(k1 = 0.1, k0 = 3, se = se, n = n, null = null,
        psd = 0.7, dpm = dpm, dpsd = dpsd, type = "moment",
        alternative = alternative, ngrid = 50000)
    mean <- direction*(dpm - null)/se
    variance <- 1 + (dpsd/se)^2
    covariance <- se[2]/se[1] + dpsd^2/(se[1]*se[2])
    a <- direction*design$zk0
    b <- direction*design$zk1
    firstH1 <- pnorm(b[1], mean[1], sqrt(variance[1]), lower.tail = FALSE)
    firstH0 <- pnorm(a[1], mean[1], sqrt(variance[1]))
    conditionalSD <- sqrt(variance[2] - covariance^2/variance[1])
    additionalH1 <- integrate(function(u1) {
        conditionalMean <- mean[2] + covariance/variance[1]*(u1 - mean[1])
        dnorm(u1, mean[1], sqrt(variance[1]))*
            pnorm(b[2], conditionalMean, conditionalSD, lower.tail = FALSE)
    }, a[1], b[1], rel.tol = 1e-10)$value
    additionalH0 <- integrate(function(u1) {
        conditionalMean <- mean[2] + covariance/variance[1]*(u1 - mean[1])
        dnorm(u1, mean[1], sqrt(variance[1]))*
            pnorm(a[2], conditionalMean, conditionalSD)
    }, a[1], b[1], rel.tol = 1e-10)$value
    expect_equal(design$cumpH1, c(firstH1, firstH1 + additionalH1), tolerance = 1e-5)
    expect_equal(design$cumpH0, c(firstH0, firstH0 + additionalH0), tolerance = 1e-5)
}
