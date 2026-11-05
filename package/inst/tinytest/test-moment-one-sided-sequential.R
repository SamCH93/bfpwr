library(tinytest)
library(bfpwr)

## One look equals the fixed design, including a nonzero null and a design
## prior with mass on both sides. Explicit two-sided calls are unchanged.
for (alternative in c("less", "greater", "two.sided")) {
    direction <- if (alternative == "less") -1 else 1
    for (dpsd in c(0, 0.3)) {
        sequential <- pbf01seq(k1 = 0.1, k0 = 3, se = 0.2, null = 2,
            psd = 0.7, dpm = 2 + direction*0.5, dpsd = dpsd,
            type = "moment", alternative = alternative)
        fixed <- pnmbf01(k = c(0.1, 3), n = 25, usd = 1, null = 2,
            psd = 0.7, dpm = 2 + direction*0.5, dpsd = dpsd,
            alternative = alternative, lower.tail = c(TRUE, FALSE))
        expect_equal(c(sequential$cumpH1, sequential$cumpH0), fixed,
                     tolerance = 1e-10)
    }
}

## Independent normal increments, evaluated with the posterior half-line
## second-moment formula rather than package Bayes factors or boundaries.
set.seed(42251)
nsim <- 100000
n <- c(15, 40, 80)
usd <- sqrt(2)
psd <- 0.7
theta <- rnorm(nsim, 0.4, 0.3)
sumY <- numeric(nsim)
active <- rep(TRUE, nsim)
stopN <- rep(max(n), nsim)
stopH1 <- stopH0 <- integer(nsim)
for (i in seq_along(n)) {
    dn <- if (i == 1) n[i] else n[i] - n[i - 1]
    sumY <- sumY + rnorm(nsim, dn*theta, usd*sqrt(dn))
    estimate <- sumY/n[i]
    se <- usd/sqrt(n[i])
    postSD <- 1/sqrt(1/se^2 + 1/psd^2)
    postMean <- estimate*postSD^2/se^2
    q <- postMean/postSD
    halfSecondMoment <- (postMean^2 + postSD^2)*pnorm(q) +
        postMean*postSD*dnorm(q)
    logbf <- dnorm(estimate, 0, se, log = TRUE) -
        dnorm(estimate, 0, sqrt(se^2 + psd^2), log = TRUE) -
        log(2*halfSecondMoment/psd^2)
    h1 <- active & logbf <= log(0.1)
    h0 <- active & logbf >= log(3)
    stopH1[h1] <- i
    stopH0[h0] <- i
    stopN[h1 | h0] <- n[i]
    active[h1 | h0] <- FALSE
}
for (alternative in c("greater", "less")) {
    direction <- if (alternative == "greater") 1 else -1
    design <- pbf01seq(k1 = 0.1, k0 = 3, se = usd/sqrt(n), n = n,
        null = 2, psd = psd, dpm = 2 + direction*0.4, dpsd = 0.3,
        type = "moment", alternative = alternative, ngrid = 20000)
    for (i in seq_along(n)) {
        for (target in c("H1", "H0")) {
            stopped <- if (target == "H1") stopH1 else stopH0
            reference <- mean(stopped > 0 & stopped <= i)
            mcse <- sqrt(reference*(1 - reference)/nsim)
            expect_true(abs(design[[paste0("cump", target)]][i] - reference) <
                        5*mcse + 1e-4)
        }
    }
    expect_true(abs(design$EN - mean(stopN)) < 5*sd(stopN)/sqrt(nsim))
    expect_true(abs(design$VarN - var(stopN)) <
                5*sd((stopN - mean(stopN))^2)/sqrt(nsim))
    expect_equal(length(design$zk1), length(n))
    for (target in c("0", "1")) {
        expect_equal(nmbf01(2 + design[[paste0("zk", target)]]*design$se,
            se = design$se, null = 2, psd = psd, alternative = alternative),
            rep(design[[paste0("k", target)]], length(n)), tolerance = 1e-8)
    }
    printed <- capture.output(print(design))
    expect_true(any(grepl(if (direction == 1) "_+" else "_-", printed, fixed = TRUE)))

    nullDesign <- pbf01seq(k1 = 0.1, k0 = 3, se = usd/sqrt(n), n = n,
        null = 2, psd = psd, dpm = 2, dpsd = 0, type = "moment",
        alternative = alternative, ngrid = 17)
    curves <- plot(nullDesign, plot = FALSE)
    expect_true(!is.null(curves$pDF1))
    expect_true(!is.null(curves$pDF2))
    expect_equal(curves$pDF1, curves$pDF2)
}

## Both cached schedule paths reach the first supported sample size for H1
## and H0; the wrapper and vectorized interface preserve the alternative.
for (target in c("H1", "H0")) {
    targetPower <- if (target == "H1") 0.6 else 0.95
    for (alternative in c("greater", "less")) {
        direction <- if (alternative == "greater") 1 else -1
        dpm <- if (target == "H1") direction*0.5 else 0
        for (schedule in list(list(looks = 3, search = "exhaustive"),
                              list(minN = 10, by = 5))) {
            found <- do.call(nbf01seq, c(list(k1 = 0.1, k0 = 3, power = targetPower,
                usd = 1, psd = 0.7, dpm = dpm, dpsd = 0, type = "moment",
                alternative = alternative, target = target,
                nrange = c(10, 200), details = TRUE), schedule))
            expect_true(found$reached)
            expect_true(found$actualPower >= targetPower)
            expect_equal(found$result$alternative, alternative)
            previousN <- found$n - if (is.null(schedule$by)) 1 else schedule$by
            previous <- do.call(powerbf01seq, c(list(n = previousN, k1 = 0.1,
                k0 = 3, psd = 0.7, dpm = dpm, dpsd = 0, type = "one.sample",
                bftype = "moment", alternative = alternative,
                target = target, nrange = c(10, 200)), schedule))
            expect_true(tail(previous[[paste0("cump", target)]], 1) < targetPower)
        }
    }
}
expect_equal(nbf01seq(k1 = 0.1, power = 0.8, usd = 1, psd = 0.7,
    dpm = c(-0.5, 0.5), dpsd = 0, type = "moment",
    alternative = c("less", "greater"), looks = 1, nrange = c(2, 200)),
    nnmbf01(k = 0.1, power = 0.8, usd = 1, psd = 0.7,
    dpm = c(-0.5, 0.5), dpsd = 0, alternative = c("less", "greater")))
wrapped <- powerbf01seq(power = 0.7, k0 = 3, psd = 0.7,
    dpm = 0.5, dpsd = 0, bftype = "moment", type = "one.sample",
    alternative = "greater", looks = 3, nrange = c(10, 200))
expect_true(wrapped$solver$reached)
expect_equal(wrapped$alternative, "greater")
expect_error(pbf01seq(k1 = 0.1, se = 0.2, psd = 0, dpm = 0.5,
                      type = "moment", alternative = "greater"), "psd > 0")

## Moment priors are centered at the null; pm is ignored, as documented.
expect_equal(pbf01seq(k1 = 0.1, se = c(0.3, 0.2), psd = 0.7,
    dpm = 0.5, dpsd = 0, pm = "ignored", type = "moment",
    alternative = "greater")$cumpH1,
    pbf01seq(k1 = 0.1, se = c(0.3, 0.2), psd = 0.7,
    dpm = 0.5, dpsd = 0, type = "moment", alternative = "greater")$cumpH1)
