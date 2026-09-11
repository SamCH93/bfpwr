library(tinytest)
library(bfpwr)

k1 <- 1/2
k0 <- 2
pow <- 0.4

search <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, details = TRUE)
)

expect_true(search$reached,
            info = "sequential t search should reach H1 target")
expect_true(search$actualPower >= pow,
            info = "sequential t search achieved power should exceed target")
expect_true(inherits(search$result, "bfseqdesign"),
            info = "sequential t search details should include design object")
expect_equal(search$result$solver$n, search$n,
             info = "sequential t design should carry solver metadata")
expect_equal(search$result$tail.eps, 1e-6,
             info = "ntbf01seq search result should store default tail.eps")
expect_equal(search$result$tail.nquad, 512,
             info = "ntbf01seq search result should store default tail.nquad")

tailSearch <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, details = TRUE, tail.eps = 1e-2,
              tail.nquad = 64)
)
expect_equal(tailSearch$result$tail.eps, 1e-2,
             info = "ntbf01seq should pass custom tail.eps into computed designs")
expect_equal(tailSearch$result$tail.nquad, 64,
             info = "ntbf01seq should pass custom tail.nquad into computed designs")

narrowRangeSearch <- suppressWarnings(
    ntbf01seq(k1 = 1/10, k0 = 10, power = 0.8, dpm = 0.5, dpsd = 0.1,
              type = "two.sample", alternative = "greater", target = "H1",
              nrange = c(20, 40), strict = FALSE, trange = c(-1, 1),
              details = TRUE)
)
expect_true(is.nan(narrowRangeSearch$n),
            info = "ntbf01seq should reject searches with failed numeric t boundaries")
expect_true(
    grepl("Failed to compute", narrowRangeSearch$error, fixed = TRUE),
    info = "ntbf01seq should report the failed boundary search in solver diagnostics"
)
narrowLookSearch <- suppressWarnings(
    ntbf01seq(k1 = 0.9, k0 = 2, power = 0.8, dpm = 0.5, dpsd = 0,
              alternative = "greater", timing = c(0.25, 0.5, 1),
              nrange = c(80, 100), strict = FALSE, trange = c(-2, 1.5),
              details = TRUE)
)
expect_true(
    grepl("look 3", narrowLookSearch$error, fixed = TRUE),
    info = "cached sequential t searches should report the failed look index"
)
narrowRangeWarning <- character()
narrowRangeN <- withCallingHandlers(
    ntbf01seq(k1 = 1/10, k0 = 10, power = 0.8, dpm = 0.5, dpsd = 0.1,
              type = "two.sample", alternative = "greater", target = "H1",
              nrange = c(20, 40), strict = FALSE, trange = c(-1, 1)),
    warning = function(w) {
        narrowRangeWarning <<- c(narrowRangeWarning, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)
expect_true(is.nan(narrowRangeN) &&
                any(grepl("Failed to compute", narrowRangeWarning,
                          fixed = TRUE)),
            info = "ntbf01seq should warn with the boundary failure when details = FALSE")
expect_false(any(grepl("Power = NaN", narrowRangeWarning, fixed = TRUE)),
             info = "ntbf01seq should not warn with generic NaN power text for boundary failures")
powerBoundaryFailure <- suppressWarnings(
    try(
        powertbf01seq(k1 = 1/10, k0 = 10, power = 0.8,
                      dpm = 0.5, dpsd = 0.1, type = "two.sample",
                      alternative = "greater", target = "H1",
                      nrange = c(20, 40), strict = FALSE,
                      trange = c(-1, 1)),
        silent = TRUE
    )
)
expect_true(
    inherits(powerBoundaryFailure, "try-error") &&
        grepl("Failed to compute",
              conditionMessage(attr(powerBoundaryFailure, "condition")),
              fixed = TRUE),
    info = "powertbf01seq should preserve boundary failures from sample-size search"
)

fakeBoundaryEval <- function(n) {
    if (n < 8) {
        return(list(n = n, criterion = -0.1, power = 0.4,
                    error = NULL, result = NULL))
    }
    list(n = n, criterion = NA_real_, power = NA_real_,
         error = "boundary failure", result = NULL)
}
hiddenBoundaryFailure <- bfpwr:::.bfseq_search_before_invalid(
    evalN = fakeBoundaryEval, validN = 4, invalidN = 16
)
expect_true(is.null(hiddenBoundaryFailure$upper) &&
                identical(hiddenBoundaryFailure$limit$error,
                          "boundary failure"),
            info = "sample-size search should not hide boundary failures behind the last finite candidate")

vectorTail <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, integer = FALSE,
              tail.eps = c(1e-2, 1e-3))
)
tailOne <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, integer = FALSE, tail.eps = 1e-2)
)
tailTwo <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, integer = FALSE, tail.eps = 1e-3)
)
expect_equal(as.numeric(vectorTail), c(tailOne, tailTwo),
             info = "ntbf01seq should vectorize over tail.eps like other scalar controls")

progressEvents <- list()
progressSearch <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, details = TRUE,
              progress = function(info) {
                  progressEvents[[length(progressEvents) + 1L]] <<- info
              })
)
expect_equal(length(progressEvents), progressSearch$evaluations,
             info = "sequential t progress callback should tick once per new evaluation")
expect_equal(progressEvents[[length(progressEvents)]]$evaluations,
             progressSearch$evaluations,
             info = "sequential t progress callback should report evaluation count")
expect_equal(progressEvents[[1]]$event, "evaluate",
             info = "sequential t progress callback should report evaluate events")
expect_true(all(vapply(progressEvents, function(x) x$target == "H1", logical(1))),
            info = "sequential t progress callback should retain target label")

ticks <- 0L
tickSearch <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, details = TRUE,
              progress = function() {
                  ticks <<- ticks + 1L
              })
)
expect_equal(ticks, tickSearch$evaluations,
             info = "sequential t progress callback should support zero-argument tick functions")

h0search <- suppressWarnings(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0, dpsd = 0,
              alternative = "greater", target = "H0", looks = 1,
              nrange = c(2, 80), strict = FALSE, details = TRUE)
)
expect_true(h0search$reached,
            info = "sequential t search should reach H0 target")
expect_equal(h0search$target, "H0",
             info = "sequential t H0 search should retain target label")

powres <- suppressWarnings(
    powertbf01seq(power = pow, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, nrange = c(2, 80),
                  strict = FALSE)
)
expect_true(inherits(powres, "bfseqdesign"),
            info = "powertbf01seq should return a sequential design")
expect_equal(powres$solver$n, search$n,
             info = "powertbf01seq should use ntbf01seq search result")

powTail <- suppressWarnings(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, strict = FALSE,
                  tail.eps = 1e-2, tail.nquad = 64)
)
expect_equal(powTail$tail.eps, 1e-2,
             info = "powertbf01seq should pass custom tail.eps in fixed-n mode")
expect_equal(powTail$tail.nquad, 64,
             info = "powertbf01seq should pass custom tail.nquad in fixed-n mode")

powerProgressEvents <- list()
powerProgress <- suppressWarnings(
    powertbf01seq(power = pow, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, nrange = c(2, 80),
                  strict = FALSE,
                  progress = function(info) {
                      powerProgressEvents[[length(powerProgressEvents) + 1L]] <<- info
                  })
)
expect_true(inherits(powerProgress, "bfseqdesign"),
            info = "powertbf01seq should return a design when progress callback is supplied")
expect_equal(length(powerProgressEvents), powerProgress$solver$evaluations,
             info = "powertbf01seq should pass progress callback to sample-size search")
expect_equal(powerProgress$solver$search, "adaptive",
             info = "powertbf01seq should record adaptive search strategy")
expect_false(powerProgress$solver$firstCrossingCertified,
             info = "adaptive multi-look timing t search should mark first crossing uncertified")

exhaustivePower <- suppressWarnings(
    powertbf01seq(power = pow, k1 = k1, k0 = k0, dpm = 0.5,
                  dpsd = 0, alternative = "greater", looks = 2,
                  nrange = c(2, 80), strict = FALSE,
                  search = "exhaustive")
)
expect_equal(exhaustivePower$solver$search, "exhaustive",
             info = "powertbf01seq should record exhaustive search strategy")
expect_true(exhaustivePower$solver$firstCrossingCertified,
            info = "exhaustive timing t search should mark first crossing certified")

increaseProgressEvents <- list()
increaseSearch <- suppressWarnings(
    ntbf01seq(k1 = 1/10, k0 = 10, power = 0.8, dpm = 0, dpsd = 0,
              type = "two.sample", alternative = "two.sided",
              target = "H0", minN = 20, by = 20,
              nrange = c(20, 800), strict = FALSE, details = TRUE,
              progress = function(info) {
                  increaseProgressEvents[[length(increaseProgressEvents) + 1L]] <<- info
              })
)
increaseFull <- suppressWarnings(
    ptbf01seq(k1 = 1/10, k0 = 10, n1 = increaseSearch$result$n1,
              n2 = increaseSearch$result$n2, dpm = 0, dpsd = 0,
              type = "two.sample", alternative = "two.sided",
              strict = FALSE, trange = "adaptive")
)
expect_equal(increaseSearch$actualPower, utils::tail(increaseFull$cumpH0, 1),
             tolerance = 1e-10,
             info = "incremental t increase search should match full ptbf01seq at selected maximum")
expect_true((increaseSearch$n - 20) %% 20 == 0,
            info = "incremental t increase search should return a scheduled maximum")
expect_true(increaseSearch$firstCrossingCertified,
            info = "incremental t increase search should certify first scheduled crossing")
expect_equal(length(increaseProgressEvents), increaseSearch$evaluations,
             info = "incremental t increase search should keep progress ticks per candidate maximum")
if (increaseSearch$n > 20) {
    prevN <- increaseSearch$n - 20
    prevFull <- suppressWarnings(
        ptbf01seq(k1 = 1/10, k0 = 10, n1 = seq(20, prevN, by = 20),
                  n2 = seq(20, prevN, by = 20), dpm = 0, dpsd = 0,
                  type = "two.sample", alternative = "two.sided",
                  strict = FALSE, trange = "adaptive")
    )
    expect_true(utils::tail(prevFull$cumpH0, 1) < 0.8,
                info = "incremental t increase search should find first scheduled H0 crossing")
}

fixedIncrease <- suppressWarnings(
    powertbf01seq(n = increaseSearch$n, k1 = 1/10, k0 = 10,
                  dpm = 0, dpsd = 0, type = "two.sample",
                  alternative = "two.sided", target = "H0",
                  minN = 20, by = 20, nrange = c(20, 800),
                  strict = FALSE)
)
expect_equal(fixedIncrease$n1, increaseSearch$result$n1,
             info = "powertbf01seq fixed-n mode should preserve searched increment n1 schedule")
expect_equal(fixedIncrease$n2, increaseSearch$result$n2,
             info = "powertbf01seq fixed-n mode should preserve searched increment n2 schedule")
expect_equal(utils::tail(fixedIncrease$cumpH0, 1), increaseSearch$actualPower,
             tolerance = 1e-10,
             info = "powertbf01seq fixed-n mode should round-trip searched increment power")

tTimingSchedule <- bfpwr:::.bfseq_schedule_spec(
    looks = 2, timing = c(0.4, 1), nrange = c(2, 80)
)
tTimingEval <- bfpwr:::.bfseq_t_schedule_evaluator(
    k1 = k1, k0 = k0, plocation = 0, pscale = 1/sqrt(2), pdf = 1,
    dpm = 0.5, dpsd = 0, type = "two.sample",
    alternative = "two.sided", target = "H1", ratio = 1,
    schedule = tTimingSchedule, strict = FALSE, trange = "adaptive",
    tail.eps = 1e-2, tail.nquad = 64, dots = list()
)
for (maxN in c(30, 43)) {
    tTimingN <- bfpwr:::.bfseq_schedule_n(maxN = maxN,
                                          schedule = tTimingSchedule)
    tTimingDirect <- suppressWarnings(
        ptbf01seq(k1 = k1, k0 = k0, n1 = tTimingN, n2 = tTimingN,
                  plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                  dpm = 0.5, dpsd = 0, type = "two.sample",
                  alternative = "two.sided", strict = FALSE,
                  trange = "adaptive", tail.eps = 1e-2,
                  tail.nquad = 64)
    )
    tTimingCached <- suppressWarnings(tTimingEval(maxN))
    expect_equal(tTimingCached$result$cumpH1, tTimingDirect$cumpH1,
                 tolerance = 1e-10,
                 info = paste("cached t timing evaluator should match ptbf01seq at maxN =", maxN))
    expect_equal(tTimingCached$result$cumpH0, tTimingDirect$cumpH0,
                 tolerance = 1e-10,
                 info = paste("cached t timing evaluator should match H0 probabilities at maxN =", maxN))
    expect_equal(tTimingCached$power, utils::tail(tTimingDirect$cumpH1, 1),
                 tolerance = 1e-10,
                 info = paste("cached t timing evaluator should return target probability at maxN =", maxN))
    expect_equal(tTimingCached$result$tail.eps, 1e-2,
                 info = paste("cached t timing evaluator should retain tail.eps at maxN =", maxN))
    expect_equal(tTimingCached$result$tail.nquad, 64,
                 info = paste("cached t timing evaluator should retain tail.nquad at maxN =", maxN))
}

ratiores <- suppressWarnings(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, ratio = 2,
                  strict = FALSE)
)
expect_equal(ratiores$n2, ceiling(ratiores$n1*2),
             info = "powertbf01seq should apply two-sample allocation ratio")
expect_equal(ratiores$ratio, 2,
             info = "powertbf01seq should retain allocation ratio")

smallRatio <- suppressWarnings(
    powertbf01seq(n = 25, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", by = 1, minN = 6,
                  ratio = 0.2, strict = FALSE)
)
expect_true(inherits(smallRatio, "bfseqdesign"),
            info = "powertbf01seq fixed-n mode should allow non-decreasing n2 schedules")
expect_true(any(diff(smallRatio$n2) == 0),
            info = "small allocation ratios can validly repeat group-2 sample sizes")

fixedTicks <- 0L
suppressWarnings(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, strict = FALSE,
                  progress = function() fixedTicks <<- fixedTicks + 1L)
)
expect_equal(fixedTicks, 0L,
             info = "powertbf01seq fixed-n mode should not call search progress callback")

expect_error(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, progress = 1),
    "progress",
    info = "sequential t progress callback should be NULL or a function"
)

expect_error(
    ntbf01seq(k1 = k1, k0 = k0, power = pow, dpm = 0.5, dpsd = 0,
              alternative = "greater", looks = 2, nrange = c(2, 80),
              strict = FALSE, search = "bad"),
    "should be one of",
    info = "sequential t search strategy should be validated"
)

expect_error(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, strict = FALSE,
                  progress = 1),
    "progress",
    info = "powertbf01seq should validate progress in fixed-n mode"
)

expect_error(
    powertbf01seq(n = 20, k1 = k1, k0 = k0, dpm = 0.5, dpsd = 0,
                  alternative = "greater", looks = 2, strict = FALSE,
                  search = "bad"),
    "should be one of",
    info = "powertbf01seq should validate search strategy in fixed-n mode"
)
