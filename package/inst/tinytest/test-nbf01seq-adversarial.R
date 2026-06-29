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
expect_true(inherits(badK1, "try-error"),
            info = "z search should reject non-positive H1 BF thresholds")

lowerSchedule <- bfpwr:::.bfseq_schedule_spec(looks = 1,
                                              nrange = c(2, 5))
lowerDip <- suppressWarnings(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 5),
        schedule = lowerSchedule,
        evaluate = function(maxN) {
            power <- if (maxN == 3) 0.5 else 0.9
            list(result = list(cumpH1 = power, cumpH0 = 0, n = maxN),
                 power = power, schedule = maxN)
        }
    )
)
expect_equal(lowerDip$n, 2,
             info = "search should return the lower bound when it already reaches")
expect_true(lowerDip$reached,
            info = "lower-bound first crossing should be considered reached")

syntheticSchedule <- bfpwr:::.bfseq_schedule_spec(looks = 1,
                                                  nrange = c(2, 10))
synthetic <- suppressWarnings(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 10),
        schedule = syntheticSchedule,
        evaluate = function(maxN) {
            list(result = list(cumpH1 = if (maxN >= 9) 0.9 else 0.5,
                               cumpH0 = 0, n = maxN),
                 power = if (maxN >= 9) 0.9 else 0.5,
                 schedule = maxN)
        }
    )
)
expect_equal(synthetic$n, 9,
             info = "search should return the first candidate that reaches")
expect_true(synthetic$reached,
            info = "first-crossing search should not require post-crossing stability")

searchPolicySchedule <- bfpwr:::.bfseq_schedule_spec(looks = 1,
                                                     nrange = c(2, 20))
searchPolicyResult <- function(maxN, power) {
    list(result = list(cumpH1 = power, cumpH0 = 0, n = maxN),
         power = power)
}

ordinarySearchError <- try(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 20),
        schedule = searchPolicySchedule,
        evaluate = function(maxN) stop("ordinary evaluator bug")
    ),
    silent = TRUE
)
expect_true(
    inherits(ordinarySearchError, "try-error") &&
        grepl("ordinary evaluator bug",
              conditionMessage(attr(ordinarySearchError, "condition")),
              fixed = TRUE),
    info = "ordinary sequential evaluator errors should propagate"
)

malformedSearchResult <- try(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 20),
        schedule = searchPolicySchedule,
        evaluate = function(maxN) list(power = 0.1)
    ),
    silent = TRUE
)
expect_true(
    inherits(malformedSearchResult, "try-error") &&
        grepl("missing 'result'",
              conditionMessage(attr(malformedSearchResult, "condition")),
              fixed = TRUE),
    info = "malformed sequential evaluator results should remain structural errors"
)

terminalInvalidSearch <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 20),
    schedule = searchPolicySchedule,
    evaluate = function(maxN) {
        if (maxN >= 4) {
            bfpwr:::.bfseq_candidate_invalid(
                "synthetic terminal invalid",
                reason = "synthetic_terminal",
                terminal = TRUE
            )
        }
        searchPolicyResult(maxN, 0.2)
    }
)
expect_false(terminalInvalidSearch$reached,
             info = "terminal invalid candidates should stop adaptive search")
expect_true(is.nan(terminalInvalidSearch$n),
            info = "terminal invalid candidates should return unresolved n")
expect_equal(terminalInvalidSearch$reason, "synthetic_terminal",
             info = "terminal invalid candidate reason should be preserved")
expect_true(terminalInvalidSearch$terminal,
            info = "terminal invalid candidate status should be preserved")

transientAdaptiveSearch <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 20),
    schedule = searchPolicySchedule,
    evaluate = function(maxN) {
        if (maxN == 4) {
            bfpwr:::.bfseq_candidate_invalid(
                "synthetic transient invalid",
                reason = "synthetic_transient",
                terminal = FALSE
            )
        }
        searchPolicyResult(maxN, if (maxN >= 8) 0.9 else 0.2)
    }
)
expect_false(transientAdaptiveSearch$reached,
             info = "adaptive search should stop at transient invalid candidates")
expect_true(
    grepl("search = \"exhaustive\"", transientAdaptiveSearch$error,
          fixed = TRUE),
    info = "adaptive transient-invalid diagnostics should recommend exhaustive search"
)
expect_false(transientAdaptiveSearch$terminal,
             info = "transient invalid candidate status should be preserved")

transientExhaustiveSearch <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 20),
    schedule = searchPolicySchedule,
    search = "exhaustive",
    evaluate = function(maxN) {
        if (maxN == 4) {
            bfpwr:::.bfseq_candidate_invalid(
                "synthetic transient invalid",
                reason = "synthetic_transient",
                terminal = FALSE
            )
        }
        searchPolicyResult(maxN, if (maxN == 8) 0.9 else 0.2)
    }
)
expect_equal(transientExhaustiveSearch$n, 8,
             info = "exhaustive search should scan past transient invalid candidates")
expect_true(transientExhaustiveSearch$reached,
            info = "exhaustive search should return a later finite crossing")
expect_false(
    transientExhaustiveSearch$firstCrossingCertified,
    info = "skipped invalid candidates should prevent absolute first-crossing certification"
)

terminalExhaustiveSearch <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 20),
    schedule = searchPolicySchedule,
    search = "exhaustive",
    evaluate = function(maxN) {
        if (maxN == 4) {
            bfpwr:::.bfseq_candidate_invalid(
                "synthetic terminal invalid",
                reason = "synthetic_terminal",
                terminal = TRUE
            )
        }
        searchPolicyResult(maxN, if (maxN >= 8) 0.9 else 0.2)
    }
)
expect_false(terminalExhaustiveSearch$reached,
             info = "exhaustive search should stop at terminal invalid candidates")
expect_true(terminalExhaustiveSearch$terminal,
            info = "exhaustive terminal-invalid diagnostics should retain terminal status")
expect_true(
    grepl("synthetic terminal invalid", terminalExhaustiveSearch$error,
          fixed = TRUE),
    info = "exhaustive terminal-invalid diagnostics should retain the invalid message"
)

nonfinitePowerSearch <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 20),
    schedule = searchPolicySchedule,
    evaluate = function(maxN) {
        searchPolicyResult(maxN, NaN)
    }
)
expect_false(nonfinitePowerSearch$reached,
             info = "non-finite evaluator power should be classified as invalid")
expect_equal(nonfinitePowerSearch$reason, "nonfinite_power",
             info = "non-finite evaluator power should retain a structured reason")
expect_true(nonfinitePowerSearch$terminal,
            info = "non-finite evaluator power should be terminal by default")

islandAdaptiveSearch <- suppressWarnings(
    bfpwr:::.bfseq_search(
        power = 0.8,
        target = "h1",
        nrange = c(2, 20),
        schedule = searchPolicySchedule,
        evaluate = function(maxN) {
            searchPolicyResult(maxN, if (maxN == 5) 0.9 else 0.2)
        }
    )
)
islandExhaustiveSearch <- bfpwr:::.bfseq_search(
    power = 0.8,
    target = "h1",
    nrange = c(2, 20),
    schedule = searchPolicySchedule,
    search = "exhaustive",
    evaluate = function(maxN) {
        searchPolicyResult(maxN, if (maxN == 5) 0.9 else 0.2)
    }
)
expect_true(is.nan(islandAdaptiveSearch$n),
            info = "adaptive bracketing can miss an isolated early crossing")
expect_equal(islandExhaustiveSearch$n, 5,
             info = "exhaustive search should scan the full candidate range")
expect_true(islandExhaustiveSearch$firstCrossingCertified,
            info = "full-range exhaustive search should certify isolated first crossings")
