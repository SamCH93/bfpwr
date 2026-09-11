library(tinytest)
library(bfpwr)

## A single look must match fixed-design probabilities on either side of a
## nonzero point null, for point and continuous design and analysis priors.
for (alternative in c("less", "greater")) {
    direction <- if (alternative == "greater") 1 else -1
    for (psd in c(0, 0.7)) {
        for (dpsd in c(0, 0.3)) {
            args <- list(null = 2, pm = 2 + direction*0.3, psd = psd,
                         dpm = 2 + direction*0.5, dpsd = dpsd,
                         alternative = alternative)
            sequential <- do.call(pbf01seq, c(list(k1 = 0.1, k0 = 3,
                                                   se = 0.2), args))
            fixed <- do.call(pbf01, c(list(k = c(0.1, 3), n = 50,
                usd = sqrt(2), lower.tail = c(TRUE, FALSE)), args))
            expect_equal(c(sequential$cumpH1, sequential$cumpH0), fixed,
                         tolerance = 1e-10)
        }
    }
}

## Simulate independent data increments and apply the Table 1 BF formula
## directly (not package critical values). Compare first-stopping events and
## sample-size moments with the sequential integration, including a design
## prior with mass on both sides of the null.
set.seed(64365)
nsim <- 100000
n <- c(15, 40, 80)
usd <- sqrt(2)
theta <- rnorm(nsim, mean = 0.4, sd = 0.3)
sumY <- numeric(nsim)
active <- rep(TRUE, nsim)
stopN <- rep(max(n), nsim)
stopH1 <- stopH0 <- integer(nsim)
for (i in seq_along(n)) {
    dn <- if (i == 1) n[i] else n[i] - n[i - 1]
    sumY <- sumY + rnorm(nsim, mean = dn*theta, sd = usd*sqrt(dn))
    estimate <- sumY/n[i]
    se <- usd/sqrt(n[i])
    postSD <- 1/sqrt(1/se^2 + 1/0.7^2)
    postMean <- (estimate/se^2 - 0.2/0.7^2)*postSD^2
    logbf <- dnorm(estimate, 0, se, log = TRUE) -
        dnorm(estimate, -0.2, sqrt(se^2 + 0.7^2), log = TRUE) +
        pnorm(-0.2/0.7, log.p = TRUE) - pnorm(postMean/postSD, log.p = TRUE)
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
                       pm = -direction*0.2, psd = 0.7,
                       dpm = direction*0.4, dpsd = 0.3,
                       alternative = alternative, ngrid = 20000)
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
    expect_equal(design$alternative, alternative)
    expect_equal(length(design$zk1), length(n))
    printed <- capture.output(print(design))
    expect_true(any(grepl(paste0("parameter ", if (direction == 1) ">" else "<"),
                         printed, fixed = TRUE)))
}

## Both cached search paths must use one-sided boundaries and return a design
## that reaches the target. Exhaustive search checks the first integer crossing.
for (target in c("H1", "H0")) {
    dpm <- if (target == "H1") 0.5 else 0
    args <- list(k1 = 0.1, k0 = 3, usd = 1, pm = 0, psd = 1,
                 dpm = dpm, dpsd = 0, alternative = "greater")
    for (schedule in list(list(looks = 3, search = "exhaustive"),
                          list(minN = 10, by = 5))) {
        found <- do.call(nbf01seq, c(list(power = 0.7, target = target,
                         nrange = c(10, 150), details = TRUE), args, schedule))
        expect_true(found$reached)
        expect_true(found$actualPower >= 0.7)
        expect_equal(found$result$alternative, "greater")
        previousN <- found$n - if (is.null(schedule$by)) 1 else schedule$by
        previous <- do.call(powerbf01seq, c(list(n = previousN,
            type = "one.sample", nrange = c(10, 150), target = target),
            args[setdiff(names(args), "usd")], schedule))
        expect_true(tail(previous[[paste0("cump", target)]], 1) < 0.7)
    }
}
expect_equal(nbf01seq(k1 = 0.1, power = 0.8, pm = 0, psd = 1,
    dpm = c(-0.5, 0.5), dpsd = 0, alternative = c("less", "greater"),
    looks = 1, nrange = c(2, 150)), rep(101, 2))
wrapped <- powerbf01seq(power = 0.8, pm = 0, psd = 1,
                        dpm = 0.5, dpsd = 0, alternative = "greater",
                        looks = 3, nrange = c(2, 150))
expect_true(wrapped$solver$reached)
expect_equal(wrapped$alternative, "greater")

## Null-plot reconstruction retains the alternative and integration settings.
for (alternative in c("less", "greater")) {
    design <- pbf01seq(k1 = 0.1, k0 = 3, se = sqrt(2/n), n = n,
                       null = 2, pm = 2, psd = 1, dpm = 2, dpsd = 0,
                       alternative = alternative, ngrid = 17)
    curves <- plot(design, plot = FALSE)
    expect_equal(curves$pDF1, curves$pDF2)
    expect_true(is.list(plot(design, plot = FALSE, zplot = TRUE)))
}

## Directional composite hypotheses and moment priors keep their existing
## meaning; these are not point-null one-sided normal designs.
for (type in c("directional", "moment")) {
    expect_error(pbf01seq(k1 = 0.1, se = 0.2, pm = 0, psd = 1,
        dpm = 0.5, type = type, alternative = "greater"), "normal analysis prior")
    expect_error(nbf01seq(k1 = 0.1, power = 0.8, pm = 0, psd = 1,
        dpm = 0.5, type = type, alternative = "greater"), "normal analysis prior")
}
expect_error(pbf01seq(k1 = 0.1, se = 0.2, pm = 0, psd = 0,
                      alternative = "greater"), "point prior")
