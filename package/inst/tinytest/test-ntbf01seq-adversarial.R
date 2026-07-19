library(tinytest)
library(bfpwr)

k1 <- 1/2
k0 <- 2

dense <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", timing = c(0.9, 1),
              nrange = c(2, 80), strict = FALSE, details = TRUE)
)
expect_equal(dense$nrange[1], 10,
             info = "dense t timing search should start at first feasible max N")
expect_true(all(diff(dense$result$n1) > 0),
            info = "dense t timing search should generate increasing n1 looks")

defaultSearchLower <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 1, nrange = c(2, 80),
              strict = FALSE, details = TRUE)
)
expect_equal(defaultSearchLower$nrange[1], 2,
             info = "t search should respect the requested lower bound")

raisedSearchLower <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 1, nrange = c(12, 80),
              strict = FALSE, details = TRUE)
)
expect_equal(raisedSearchLower$nrange[1], 12,
             info = "t non-increment search should respect a higher nrange lower bound")

customIncrement <- suppressWarnings(
    powertbf01seq(n = 23, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", by = 10, minN = 5,
                  strict = FALSE)
)
expect_equal(customIncrement$n1, c(5, 15, 23),
             info = "t increment schedule should respect custom first look")

defaultIncrement <- suppressWarnings(
    powertbf01seq(n = 23, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", by = 10, strict = FALSE)
)
expect_equal(defaultIncrement$n1, c(2, 12, 22, 23),
             info = "t increment schedule without minN should start at lower bound")

shortRangeIncrement <- suppressWarnings(
    powertbf01seq(n = 5, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", by = 10, nrange = c(2, 5),
                  strict = FALSE)
)
expect_equal(shortRangeIncrement$n1, c(2, 5),
             info = "t increment schedule should use lower bound when nrange is below 10")

highRangeIncrement <- suppressWarnings(
    powertbf01seq(n = 23, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", by = 10, nrange = c(12, 100),
                  strict = FALSE)
)
expect_equal(highRangeIncrement$n1, c(12, 22, 23),
             info = "t increment schedule should use lower bound when nrange starts above 10")

badRatio <- try(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", ratio = 0, strict = FALSE),
    silent = TRUE
)
expect_true(inherits(badRatio, "try-error"),
            info = "t wrapper should reject non-positive allocation ratios")

tinyRatio <- try(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", ratio = 0.01, strict = FALSE),
    silent = TRUE
)
expect_true(inherits(tinyRatio, "try-error"),
            info = "t wrapper should reject ratios producing n2 below two")

badK1 <- try(
    ntbf01seq(k1 = 0, k0 = 2, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", strict = FALSE),
    silent = TRUE
)
expect_true(inherits(badK1, "try-error"),
            info = "t search should reject non-positive H1 BF thresholds")

expect_error(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", nrange = c(2.2, 2.8),
              strict = FALSE),
    "no integer candidate",
    info = "t search should not evaluate an integer outside nrange"
)

transientBoundarySearch <- suppressWarnings(
    ntbf01seq(
        k1 = 0.05, k0 = 2, power = 0.1, dpm = 0.5, dpsd = 0,
        alternative = "greater", target = "H1", looks = 1,
        nrange = c(2, 80), trange = c(-100, 3),
        search = "exhaustive", details = TRUE
    )
)
expect_equal(transientBoundarySearch$n, 33,
             info = "exhaustive t search should scan past transient boundary failures")
expect_true(transientBoundarySearch$reached,
            info = "exhaustive t search should find a later valid target crossing")
expect_false(transientBoundarySearch$firstCrossingCertified,
             info = "skipped t boundary failures should prevent certification")

smallRatio <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", ratio = 0.1, nrange = c(2, 200),
              strict = FALSE, details = TRUE)
)
expect_true(smallRatio$reached,
            info = "valid small t allocation ratio should still be searchable")
expect_equal(smallRatio$nrange[1], 11,
             info = "small t allocation ratio should raise the lower search bound")
expect_true(all(smallRatio$result$n2 >= 2),
            info = "small t allocation ratio should keep all n2 looks feasible")

unreachableWarning <- NULL
unreachable <- withCallingHandlers(
    suppressWarnings(
        ntbf01seq(k1 = 1/10, k0 = 10, power = 0.99, dpm = 0, dpsd = 0,
                  alternative = "greater", nrange = c(2, 5),
                  strict = FALSE, details = TRUE)
    ),
    warning = function(w) {
        unreachableWarning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
    }
)
expect_true(is.nan(unreachable$n),
            info = "unreachable t target should return NaN sample size")
expect_false(unreachable$reached,
             info = "unreachable t target should be flagged as unreached")

vecn <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = c(0.3, 0.4), dpm = 0.5,
              dpsd = 0, alternative = "greater", nrange = c(2, 80),
              strict = FALSE)
)
expect_equal(length(vecn), 2,
             info = "t search should retain vectorized numeric output")
expect_true(all(is.finite(vecn)),
            info = "t vectorized search should return finite values here")

oneLookH1 <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 1, nrange = c(2, 80),
              strict = FALSE)
)
fixedH1 <- suppressWarnings(
    ntbf01(k = k1, power = 0.4, dpm = 0.5, dpsd = 0,
           alternative = "greater", nrange = c(2, 80))
)
oneLookH0 <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0, dpsd = 0,
              alternative = "greater", target = "H0", looks = 1,
              nrange = c(2, 80), strict = FALSE)
)
expect_equal(oneLookH1, fixedH1,
             info = "one-look t H1 search should match fixed-design search")
expect_equal(oneLookH0, 3,
             info = "one-look t H0 search should return the first requested-range crossing")

detailsVector <- try(
    ntbf01seq(k1 = c(1/2, 1/3), k0 = 2, power = 0.4, dpm = 0.5,
              dpsd = 0, alternative = "greater", strict = FALSE,
              details = TRUE),
    silent = TRUE
)
expect_true(inherits(detailsVector, "try-error"),
            info = "details mode should remain scalar for t search")

detailsTargetVector <- try(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", target = c("H1", "H0"),
              strict = FALSE, details = TRUE),
    silent = TRUE
)
expect_true(inherits(detailsTargetVector, "try-error"),
            info = "details mode should reject vectorized t targets")

badTiming <- try(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", timing = c(0.5, 0.5, 1),
              strict = FALSE),
    silent = TRUE
)
badScheduleMix <- try(
    ntbf01seq(k1 = k1, k0 = k0, power = 0.4, dpm = 0.5, dpsd = 0,
              alternative = "greater", timing = c(0.5, 1), by = 5,
              strict = FALSE),
    silent = TRUE
)
expect_true(inherits(badTiming, "try-error"),
            info = "t timing should reject duplicate information fractions")
expect_true(inherits(badScheduleMix, "try-error"),
            info = "t schedule should reject timing and increment together")

greaterTimingEvents <- list()
greaterTiming <- suppressWarnings(
    ntbf01seq(k1 = 1/10, k0 = 10, power = 0.8, dpm = 0.5, dpsd = 0,
              type = "two.sample", alternative = "greater", ratio = 1,
              timing = seq(0.2, 1, length.out = 5), nrange = c(20, 200),
              strict = TRUE, trange = "adaptive", details = TRUE,
              progress = function(info) {
                  greaterTimingEvents[[length(greaterTimingEvents) + 1L]] <<- info
              })
)
expect_equal(greaterTiming$n, 96,
             info = "one-sided strict t sequential search should keep accelerated adaptive result")
expect_true(greaterTiming$actualPower >= 0.8,
            info = "one-sided strict t sequential search should reach the target")
expect_equal(length(greaterTimingEvents), greaterTiming$evaluations,
             info = "one-sided strict t sequential search should still report progress")
