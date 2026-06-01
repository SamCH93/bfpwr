library(tinytest)
library(bfpwr)

## Dense information fractions should advance the lower bound to the first
## feasible maximum N instead of wasting evaluations on duplicate schedules.
dense <- nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2), pm = 0,
                  psd = 1, dpm = 0.5, dpsd = 0, timing = c(0.9, 1),
                  nrange = c(2, 80), details = TRUE)
expect_equal(dense$nrange[1], 10,
             info = "dense timing search should start at first feasible max N")
expect_true(all(diff(dense$result$n) > 0),
            info = "dense timing search should generate increasing looks")
expect_true(min(dense$result$n) >= 2,
            info = "dense timing search should keep look sizes feasible")

increment <- powerbf01seq(n = 23, k1 = 1/2, k0 = 2, pm = 0, psd = 1,
                          dpm = 0.5, dpsd = 0, by = 10, minN = 5)
expect_equal(increment$n, c(5, 15, 23),
             info = "z increment schedule should append requested final n")

unreachableWarning <- NULL
unreachable <- withCallingHandlers(
    nbf01seq(k1 = 1/10, k0 = 10, power = 0.99, usd = sqrt(2), pm = 0,
             psd = 1, dpm = 0, dpsd = 0, nrange = c(2, 5),
             details = TRUE),
    warning = function(w) {
        unreachableWarning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    }
)
expect_true(is.nan(unreachable$n),
            info = "unreachable z target should return NaN sample size")
expect_false(unreachable$reached,
             info = "unreachable z target should be flagged as unreached")
expect_true(grepl("upper bound", unreachableWarning, fixed = TRUE),
            info = "unreachable z target should warn about upper nrange bound")

vecn <- nbf01seq(k1 = 1/2, k0 = 2, power = c(0.3, 0.4), usd = sqrt(2),
                 pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
                 nrange = c(2, 80))
expect_equal(length(vecn), 2,
             info = "z search should retain vectorized numeric output")
expect_true(all(is.finite(vecn)),
            info = "z vectorized search should return finite values here")

vecpm <- nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2),
                  pm = c(0, 0.1), psd = 1, dpm = 0.5, dpsd = 0,
                  nrange = c(2, 100))
expect_equal(length(vecpm), 2,
             info = "z search should vectorize over analysis prior mean")
expect_true(all(is.finite(vecpm)),
            info = "z pm-vectorized search should return finite values here")

oneLookH1 <- nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2),
                      pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
                      looks = 1, nrange = c(2, 100))
fixedH1 <- nbf01(k = 1/2, power = 0.4, usd = sqrt(2), pm = 0,
                 psd = 1, dpm = 0.5, dpsd = 0, nrange = c(2, 100),
                 analytical = FALSE)
oneLookH0 <- nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2),
                      pm = 0, psd = 1, dpm = 0, dpsd = 0, target = "h0",
                      looks = 1, nrange = c(2, 100))
fixedH0 <- nbf01(k = 2, power = 0.4, usd = sqrt(2), pm = 0, psd = 1,
                 dpm = 0, dpsd = 0, lower.tail = FALSE,
                 nrange = c(2, 100), analytical = FALSE)
expect_equal(oneLookH1, fixedH1,
             info = "one-look z H1 search should match fixed-design search")
expect_equal(oneLookH0, fixedH0,
             info = "one-look z H0 search should match fixed-design search")

early <- nbf01seq(k1 = 1/3, k0 = 5, power = 0.9826,
                  pm = 0.5, psd = 2, dpm = 1, dpsd = 0, looks = 4,
                  nrange = c(2, 120), search = "exhaustive",
                  details = TRUE)
expect_equal(early$n, 44,
             info = "timing search should return the first crossing")
expect_true(early$actualPower >= 0.9826,
            info = "timing first crossing should reach the target")
expect_true(early$firstCrossingCertified,
            info = "exhaustive timing search should certify the first crossing")

timingSchedule <- bfpwr:::.bfseq_schedule_spec(looks = 5,
                                               nrange = c(2, 1000))
stepEvaluator <- function(maxN) {
    power <- if (maxN >= 800) 0.9 else 0.5
    list(result = list(cumpH1 = power, cumpH0 = 0, n = maxN),
         power = power)
}
adaptiveStep <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 1000),
    schedule = timingSchedule,
    evaluate = stepEvaluator,
    search = "adaptive"
)
exhaustiveStep <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 1000),
    schedule = timingSchedule,
    evaluate = stepEvaluator,
    search = "exhaustive"
)
expect_equal(adaptiveStep$n, exhaustiveStep$n,
             info = "adaptive timing search should find the bracketing solution in monotone cases")
expect_true(adaptiveStep$evaluations < exhaustiveStep$evaluations/10,
            info = "adaptive timing search should avoid exhaustive first-crossing scans")
expect_false(adaptiveStep$firstCrossingCertified,
             info = "adaptive multi-look timing search should report uncertified first crossing")
expect_true(exhaustiveStep$firstCrossingCertified,
            info = "exhaustive multi-look timing search should report certified first crossing")

detailsVector <- try(
    nbf01seq(k1 = c(1/2, 1/3), k0 = 2, power = 0.4, usd = sqrt(2),
             pm = 0, psd = 1, dpm = 0.5, dpsd = 0, details = TRUE),
    silent = TRUE
)
expect_true(inherits(detailsVector, "try-error"),
            info = "details mode should remain scalar for z search")

detailsTargetVector <- try(
    nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2),
             pm = 0, psd = 1, dpm = 0.5, dpsd = 0,
             target = c("h1", "h0"), details = TRUE),
    silent = TRUE
)
expect_true(inherits(detailsTargetVector, "try-error"),
            info = "details mode should reject vectorized z targets")

moment <- powerbf01seq(power = 0.4, k1 = 1/2, k0 = 2, psd = 1,
                       bftype = "moment", dpm = 0.5, dpsd = 0,
                       nrange = c(2, 80))
expect_true(inherits(moment, "bfseqdesign"),
            info = "z moment-prior wrapper should not require pm")
expect_equal(moment$type, "moment",
             info = "z moment-prior wrapper should pass through BF type")

missingMomentDpm <- try(
    powerbf01seq(power = 0.4, k1 = 1/2, k0 = 2, psd = 1,
                 bftype = "moment", nrange = c(2, 80)),
    silent = TRUE
)
expect_true(
    inherits(missingMomentDpm, "try-error") &&
        grepl("dpm", conditionMessage(attr(missingMomentDpm, "condition")),
              fixed = TRUE),
    info = "z moment-prior wrapper should require explicit dpm"
)

badTimingFlat <- try(
    nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2), pm = 0,
             psd = 1, dpm = 0.5, dpsd = 0, timing = c(0.5, 0.5, 1)),
    silent = TRUE
)
badTimingEnd <- try(
    nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2), pm = 0,
             psd = 1, dpm = 0.5, dpsd = 0, timing = c(0.5, 0.9)),
    silent = TRUE
)
badScheduleMix <- try(
    nbf01seq(k1 = 1/2, k0 = 2, power = 0.4, usd = sqrt(2), pm = 0,
             psd = 1, dpm = 0.5, dpsd = 0, timing = c(0.5, 1), by = 5),
    silent = TRUE
)
expect_true(inherits(badTimingFlat, "try-error"),
            info = "z timing should reject duplicate information fractions")
expect_true(inherits(badTimingEnd, "try-error"),
            info = "z timing should reject schedules not ending at one")
expect_true(inherits(badScheduleMix, "try-error"),
            info = "z schedule should reject timing and increment together")

badK1 <- try(
    nbf01seq(k1 = 0, k0 = 2, power = 0.4, usd = sqrt(2), pm = 0,
             psd = 1, dpm = 0.5, dpsd = 0),
    silent = TRUE
)
badFixedNextend <- try(
    powerbf01seq(n = 20, k1 = 1/2, k0 = 2, pm = 0, psd = 1,
                 dpm = 0.5, dpsd = 0, nextend = NA_real_),
    silent = TRUE
)
expect_true(inherits(badK1, "try-error"),
            info = "z search should reject non-positive H1 BF thresholds")
expect_true(inherits(badFixedNextend, "try-error"),
            info = "z fixed-n wrapper should validate nextend")

lowerSchedule <- bfpwr:::.bfseq_schedule_spec(looks = 1,
                                              nrange = c(2, 5))
lowerDip <- suppressWarnings(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 5),
        schedule = lowerSchedule,
        nextend = 1,
        evaluate = function(maxN) {
            power <- if (maxN == 3) 0.5 else 0.9
            list(result = list(cumpH1 = power, cumpH0 = 0, n = maxN),
                 power = power, schedule = maxN)
        }
    )
)
expect_equal(lowerDip$n, 4,
             info = "nextend should be checked when lower bound already reaches")
expect_true(lowerDip$reached,
            info = "lower-bound nextend certification should advance to stable n")

syntheticSchedule <- bfpwr:::.bfseq_schedule_spec(looks = 1,
                                                  nrange = c(2, 10))
synthetic <- suppressWarnings(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 10),
        schedule = syntheticSchedule,
        nextend = 3,
        evaluate = function(maxN) {
            list(result = list(cumpH1 = if (maxN >= 9) 0.9 else 0.5,
                               cumpH0 = 0, n = maxN),
                 power = if (maxN >= 9) 0.9 else 0.5,
                 schedule = maxN)
        }
    )
)
expect_false(synthetic$reached,
             info = "nextend should fail closed when stability cannot be certified")
expect_true(is.nan(synthetic$n),
            info = "uncertified nextend search should return NaN sample size")
