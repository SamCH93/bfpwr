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
expect_true(all(vapply(progressEvents, function(x) x$target == "h1", logical(1))),
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
              alternative = "greater", target = "h0", looks = 1,
              nrange = c(2, 80), strict = FALSE, details = TRUE)
)
expect_true(h0search$reached,
            info = "sequential t search should reach H0 target")
expect_equal(h0search$target, "h0",
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
              target = "h0", minN = 20, by = 20,
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

tTimingSchedule <- bfpwr:::.bfseq_schedule_spec(
    looks = 2, timing = c(0.4, 1), nrange = c(2, 80)
)
tTimingEval <- bfpwr:::.bfseq_t_schedule_evaluator(
    k1 = k1, k0 = k0, plocation = 0, pscale = 1/sqrt(2), pdf = 1,
    dpm = 0.5, dpsd = 0, type = "two.sample",
    alternative = "two.sided", target = "h1", ratio = 1,
    schedule = tTimingSchedule, strict = FALSE, trange = "adaptive",
    dots = list()
)
for (maxN in c(30, 43)) {
    tTimingN <- bfpwr:::.bfseq_schedule_n(maxN = maxN,
                                          schedule = tTimingSchedule)
    tTimingDirect <- suppressWarnings(
        ptbf01seq(k1 = k1, k0 = k0, n1 = tTimingN, n2 = tTimingN,
                  plocation = 0, pscale = 1/sqrt(2), pdf = 1,
                  dpm = 0.5, dpsd = 0, type = "two.sample",
                  alternative = "two.sided", strict = FALSE,
                  trange = "adaptive")
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
