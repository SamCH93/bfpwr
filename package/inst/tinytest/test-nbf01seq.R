library(tinytest)
library(bfpwr)

k1 <- 1/5
k0 <- 5
pow <- 0.8
usd <- sqrt(2)
pm <- 0
psd <- 1
dpm <- 0.5
dpsd <- 0

helperNull <- 0.2
helperSe <- 0.1
helperPm <- 0.6
helperPsd <- 0.8
helperK <- 3
normalZ <- bfpwr:::zcrit(k = helperK, se = helperSe, null = helperNull,
                         mu = helperPm, tau = helperPsd, type = "normal")
expect_equal(bf01(estimate = helperNull + normalZ*helperSe, se = helperSe,
                  null = helperNull, pm = helperPm, psd = helperPsd),
             rep(helperK, length(normalZ)),
             tolerance = 1e-10,
             info = "zcrit normal boundaries should use original-scale null and prior mean")
pointZ <- bfpwr:::zcrit(k = helperK, se = helperSe, null = helperNull,
                        mu = helperPm, tau = 0, type = "normal")
expect_equal(bf01(estimate = helperNull + pointZ*helperSe, se = helperSe,
                  null = helperNull, pm = helperPm, psd = 0),
             helperK,
             tolerance = 1e-10,
             info = "zcrit point-prior boundary should use original-scale null and prior mean")
directionalZ <- bfpwr:::zcrit(k = helperK, se = helperSe, null = helperNull,
                              mu = helperPm, tau = helperPsd,
                              type = "directional")
expect_equal(dirbf01(estimate = helperNull + directionalZ*helperSe,
                     se = helperSe, null = helperNull, pm = helperPm,
                     psd = helperPsd),
             helperK,
             tolerance = 1e-10,
             info = "zcrit directional boundary should use original-scale null and prior mean")
momentZ <- bfpwr:::zcrit(k = helperK, se = helperSe, null = helperNull,
                         tau = 0.5, type = "moment")
expect_equal(nmbf01(estimate = helperNull + momentZ*helperSe, se = helperSe,
                    null = helperNull, psd = 0.5),
             rep(helperK, length(momentZ)),
             tolerance = 1e-10,
             info = "zcrit moment boundaries should use original-scale null")
helperPars <- bfpwr:::predpars(se = c(helperSe, helperSe/2),
                               null = helperNull, dpm = 0.45, dpsd = 0.1)
expect_equal(helperPars$mean, (0.45 - helperNull)/c(helperSe, helperSe/2),
             tolerance = 1e-10,
             info = "predpars should use original-scale null and design prior mean")

search <- nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd, pm = pm,
                   psd = psd, dpm = dpm, dpsd = dpsd, looks = 3,
                   nrange = c(2, 200), details = TRUE)

expect_true(search$reached,
            info = "sequential z search should reach H1 target")
expect_true(search$actualPower >= pow,
            info = "sequential z search achieved power should exceed target")
expect_equal(search$n, max(search$result$n),
             info = "sequential z search n should be the final look")
expect_equal(search$result$solver$n, search$n,
             info = "sequential z design should carry solver metadata")

progressEvents <- list()
progressSearch <- nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd,
                           pm = pm, psd = psd, dpm = dpm, dpsd = dpsd,
                           looks = 3, nrange = c(2, 200), details = TRUE,
                           progress = function(info) {
                               progressEvents[[length(progressEvents) + 1L]] <<- info
                           })
expect_equal(length(progressEvents), progressSearch$evaluations,
             info = "sequential z progress callback should tick once per new evaluation")
expect_equal(progressEvents[[length(progressEvents)]]$evaluations,
             progressSearch$evaluations,
             info = "sequential z progress callback should report evaluation count")
expect_equal(progressEvents[[1]]$event, "evaluate",
             info = "sequential z progress callback should report evaluate events")
expect_true(all(vapply(progressEvents, function(x) x$target == "H1", logical(1))),
            info = "sequential z progress callback should retain target label")

ticks <- 0L
tickSearch <- nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd,
                       pm = pm, psd = psd, dpm = dpm, dpsd = dpsd,
                       looks = 3, nrange = c(2, 200), details = TRUE,
                       progress = function() {
                           ticks <<- ticks + 1L
                       })
expect_equal(ticks, tickSearch$evaluations,
             info = "sequential z progress callback should support zero-argument tick functions")

powerProgressEvents <- list()
powerProgress <- powerbf01seq(power = pow, k1 = k1, k0 = k0,
                              pm = pm, psd = psd, dpm = dpm, dpsd = dpsd,
                              looks = 3, nrange = c(2, 200),
                              progress = function(info) {
                                  powerProgressEvents[[length(powerProgressEvents) + 1L]] <<- info
                              })
expect_true(inherits(powerProgress, "bfseqdesign"),
            info = "powerbf01seq should return a design when progress callback is supplied")
expect_equal(length(powerProgressEvents), powerProgress$solver$evaluations,
             info = "powerbf01seq should pass progress callback to sample-size search")
expect_equal(powerProgress$solver$search, "adaptive",
             info = "powerbf01seq should record adaptive search strategy")
expect_false(powerProgress$solver$firstCrossingCertified,
             info = "adaptive multi-look timing search should mark first crossing uncertified")

exhaustivePower <- powerbf01seq(power = pow, k1 = k1, k0 = k0,
                                pm = pm, psd = psd, dpm = dpm,
                                dpsd = dpsd, looks = 3,
                                nrange = c(2, 200),
                                search = "exhaustive")
expect_equal(exhaustivePower$solver$search, "exhaustive",
             info = "powerbf01seq should record exhaustive search strategy")
expect_true(exhaustivePower$solver$firstCrossingCertified,
            info = "exhaustive timing search should mark first crossing certified")

increaseProgressEvents <- list()
increaseSearch <- nbf01seq(k1 = 1/10, k0 = 10, power = 0.8, usd = usd,
                           pm = 0, psd = 1/sqrt(2), dpm = 0, dpsd = 0,
                           type = "normal", target = "H0",
                           minN = 20, by = 20, nrange = c(20, 1200),
                           strict = FALSE, details = TRUE,
                           progress = function(info) {
                               increaseProgressEvents[[length(increaseProgressEvents) + 1L]] <<- info
                           })
increaseFull <- pbf01seq(k1 = 1/10, k0 = 10,
                         se = usd/sqrt(increaseSearch$result$n),
                         n = increaseSearch$result$n, pm = 0,
                         psd = 1/sqrt(2), dpm = 0, dpsd = 0,
                         type = "normal", strict = FALSE)
expect_equal(increaseSearch$actualPower, utils::tail(increaseFull$cumpH0, 1),
             tolerance = 1e-10,
             info = "incremental z increase search should match full pbf01seq at selected maximum")
expect_true((increaseSearch$n - 20) %% 20 == 0,
            info = "incremental z increase search should return a scheduled maximum")
expect_true(increaseSearch$firstCrossingCertified,
            info = "incremental z increase search should certify first scheduled crossing")
expect_equal(length(increaseProgressEvents), increaseSearch$evaluations,
             info = "incremental z increase search should keep progress ticks per candidate maximum")
if (increaseSearch$n > 20) {
    prevN <- increaseSearch$n - 20
    prevFull <- pbf01seq(k1 = 1/10, k0 = 10,
                         se = usd/sqrt(seq(20, prevN, by = 20)),
                         n = seq(20, prevN, by = 20), pm = 0,
                         psd = 1/sqrt(2), dpm = 0, dpsd = 0,
                         type = "normal", strict = FALSE)
    expect_true(utils::tail(prevFull$cumpH0, 1) < 0.8,
                info = "incremental z increase search should find first scheduled H0 crossing")
}

timingSchedule <- bfpwr:::.bfseq_schedule_spec(
    looks = 3, timing = c(0.25, 0.55, 1), nrange = c(2, 200)
)
timingEval <- bfpwr:::.bfseq_z_schedule_evaluator(
    k1 = k1, k0 = k0, usd = usd, null = 0, pm = pm, psd = psd,
    dpm = dpm, dpsd = dpsd, type = "normal", target = "H1",
    schedule = timingSchedule, strict = TRUE, dots = list()
)
for (maxN in c(80, 123)) {
    timingN <- bfpwr:::.bfseq_schedule_n(maxN = maxN,
                                         schedule = timingSchedule)
    timingDirect <- pbf01seq(k1 = k1, k0 = k0, se = usd/sqrt(timingN),
                             n = timingN, pm = pm, psd = psd,
                             dpm = dpm, dpsd = dpsd, type = "normal",
                             strict = TRUE)
    timingCached <- timingEval(maxN)
    expect_equal(timingCached$result$cumpH1, timingDirect$cumpH1,
                 tolerance = 1e-10,
                 info = paste("cached z timing evaluator should match pbf01seq at maxN =", maxN))
    expect_equal(timingCached$result$cumpH0, timingDirect$cumpH0,
                 tolerance = 1e-10,
                 info = paste("cached z timing evaluator should match H0 probabilities at maxN =", maxN))
    expect_equal(timingCached$power, utils::tail(timingDirect$cumpH1, 1),
                 tolerance = 1e-10,
                 info = paste("cached z timing evaluator should return target probability at maxN =", maxN))
}

expect_error(
    nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd, pm = pm,
             psd = psd, dpm = dpm, dpsd = dpsd, looks = 3,
             nrange = c(2, 200), progress = function(info) stop("progress stopped")),
    "progress stopped",
    info = "sequential z progress callback errors should propagate"
)

schedule <- bfpwr:::.bfseq_schedule_spec(looks = 3, nrange = c(2, 200))
if (search$n > search$nrange[1]) {
    nprev <- bfpwr:::.bfseq_schedule_n(search$n - 1, schedule)
    prev <- pbf01seq(k1 = k1, k0 = k0, se = usd/sqrt(nprev), n = nprev,
                     pm = pm, psd = psd, dpm = dpm, dpsd = dpsd)
    expect_true(utils::tail(prev$cumpH1, 1) < pow,
                info = "sequential z search should find first H1 crossing")
}

h0search <- nbf01seq(k1 = k1, k0 = k0, power = 0.7, usd = usd, pm = pm,
                     psd = psd, dpm = 0, dpsd = 0, target = "H0",
                     looks = 2, nrange = c(2, 500), details = TRUE)
expect_true(h0search$reached,
            info = "sequential z search should reach H0 target")
expect_true(h0search$actualPower >= 0.7,
            info = "sequential z H0 search achieved power should exceed target")
expect_equal(h0search$target, "H0",
             info = "sequential z H0 search should retain target label")

powres <- powerbf01seq(power = pow, k1 = k1, k0 = k0, pm = pm, psd = psd,
                       dpm = dpm, dpsd = dpsd, looks = 3,
                       nrange = c(2, 200))
expect_true(inherits(powres, "bfseqdesign"),
            info = "powerbf01seq should return a sequential design")
expect_equal(powres$solver$n, search$n,
             info = "powerbf01seq should use nbf01seq search result")
expect_equal(powres$sample.type, "two.sample",
             info = "powerbf01seq should retain sampling design")

fixed <- powerbf01seq(n = 60, k1 = k1, k0 = k0, pm = pm, psd = psd,
                      dpm = dpm, dpsd = dpsd, looks = 3)
expect_true(inherits(fixed, "bfseqdesign"),
            info = "powerbf01seq fixed-n mode should return a sequential design")
expect_equal(max(fixed$n), 60,
             info = "powerbf01seq fixed-n mode should use requested final n")

theta0 <- 0.2
shiftN <- 50
shiftPm <- 0.6
shiftDpm <- 0.45
shiftPsd <- 0.8
shiftDpsd <- 0.1
shiftSeq <- pbf01seq(k1 = 1/3, k0 = 3, se = usd/sqrt(shiftN), n = shiftN,
                     null = theta0, pm = shiftPm, psd = shiftPsd,
                     dpm = shiftDpm, dpsd = shiftDpsd, type = "normal")
expect_equal(shiftSeq$cumpH1,
             pbf01(k = 1/3, n = shiftN, usd = usd, null = theta0,
                   pm = shiftPm, psd = shiftPsd, dpm = shiftDpm,
                   dpsd = shiftDpsd),
             tolerance = 1e-10,
             info = "one-look sequential z H1 probability should match fixed pbf01 with nonzero null")
expect_equal(shiftSeq$cumpH0,
             pbf01(k = 3, n = shiftN, usd = usd, null = theta0,
                   pm = shiftPm, psd = shiftPsd, dpm = shiftDpm,
                   dpsd = shiftDpsd, lower.tail = FALSE),
             tolerance = 1e-10,
             info = "one-look sequential z H0 probability should match fixed pbf01 with nonzero null")
expect_equal(shiftSeq$null, theta0,
             info = "sequential z design should retain user-facing null value")
expect_equal(shiftSeq$pm, shiftPm,
             info = "sequential z design should retain user-facing analysis prior mean")
expect_equal(shiftSeq$dpm, shiftDpm,
             info = "sequential z design should retain user-facing design prior mean")

shiftManual <- pbf01seq(k1 = 1/3, k0 = 3, se = usd/sqrt(c(25, 50)),
                        n = c(25, 50), null = 0, pm = shiftPm - theta0,
                        psd = shiftPsd, dpm = shiftDpm - theta0,
                        dpsd = shiftDpsd, type = "directional")
shiftDirectional <- pbf01seq(k1 = 1/3, k0 = 3, se = usd/sqrt(c(25, 50)),
                             n = c(25, 50), null = theta0, pm = shiftPm,
                             psd = shiftPsd, dpm = shiftDpm,
                             dpsd = shiftDpsd, type = "directional")
expect_equal(shiftDirectional$cumpH1, shiftManual$cumpH1,
             tolerance = 1e-10,
             info = "directional sequential z design should match equivalent shifted zero-null design")
expect_equal(shiftDirectional$cumpH0, shiftManual$cumpH0,
             tolerance = 1e-10,
             info = "directional sequential z H0 probability should use the supplied null split point")

shiftMoment <- pbf01seq(k1 = 1/3, k0 = 3, se = usd/sqrt(shiftN), n = shiftN,
                        null = theta0, psd = 0.5, dpm = shiftDpm,
                        dpsd = shiftDpsd, type = "moment")
expect_equal(shiftMoment$cumpH1,
             pnmbf01(k = 1/3, n = shiftN, usd = usd, null = theta0,
                     psd = 0.5, dpm = shiftDpm, dpsd = shiftDpsd),
             tolerance = 1e-10,
             info = "moment sequential z H1 probability should match fixed pnmbf01 with nonzero null")
expect_equal(shiftMoment$cumpH0,
             pnmbf01(k = 3, n = shiftN, usd = usd, null = theta0,
                     psd = 0.5, dpm = shiftDpm, dpsd = shiftDpsd,
                     lower.tail = FALSE),
             tolerance = 1e-10,
             info = "moment sequential z H0 probability should match fixed pnmbf01 with nonzero null")

shiftMomentSearch <- nbf01seq(k1 = 1/3, k0 = 3, power = 0.1, usd = usd,
                              null = theta0, psd = 0.5, dpm = shiftDpm,
                              dpsd = shiftDpsd, type = "moment",
                              looks = 1, nrange = c(10, 500))
fixedMomentSearch <- nnmbf01(k = 1/3, power = 0.1, usd = usd,
                             null = theta0, psd = 0.5, dpm = shiftDpm,
                             dpsd = shiftDpsd, nrange = c(10, 500))
expect_equal(shiftMomentSearch, fixedMomentSearch,
             info = "one-look sequential moment z sample-size search should match fixed nnmbf01 with nonzero null")

shiftSearch <- nbf01seq(k1 = 1/3, k0 = 3, power = 0.4, usd = usd,
                        null = theta0, pm = shiftPm, psd = shiftPsd,
                        dpm = shiftDpm, dpsd = shiftDpsd, looks = 1,
                        nrange = c(10, 120))
fixedShiftSearch <- nbf01(k = 1/3, power = 0.4, usd = usd, null = theta0,
                          pm = shiftPm, psd = shiftPsd, dpm = shiftDpm,
                          dpsd = shiftDpsd, nrange = c(10, 120),
                          analytical = FALSE)
expect_equal(shiftSearch, fixedShiftSearch,
             info = "one-look sequential z sample-size search should match fixed nbf01 with nonzero null")

shiftPower <- powerbf01seq(n = shiftN, k1 = 1/3, k0 = 3, null = theta0,
                           pm = shiftPm, psd = shiftPsd, dpm = shiftDpm,
                           dpsd = shiftDpsd)
expect_equal(shiftPower$cumpH1, shiftSeq$cumpH1,
             tolerance = 1e-10,
             info = "powerbf01seq fixed-n mode should pass nonzero null to pbf01seq")
expect_equal(shiftPower$null, theta0,
             info = "powerbf01seq design should retain user-facing null value")

shiftMomentPower <- powerbf01seq(n = shiftN, k1 = 1/3, k0 = 3,
                                 null = theta0, psd = 0.5,
                                 dpm = shiftDpm, dpsd = shiftDpsd,
                                 bftype = "moment")
fixedMomentPower <- powernmbf01(n = shiftN, k = 1/3, null = theta0,
                                psd = 0.5, dpm = shiftDpm,
                                dpsd = shiftDpsd)
expect_equal(utils::tail(shiftMomentPower$cumpH1, 1),
             fixedMomentPower$power,
             tolerance = 1e-10,
             info = "powerbf01seq moment fixed-n mode should match fixed powernmbf01 with nonzero null")

fixedIncrease <- powerbf01seq(n = increaseSearch$n, k1 = 1/10, k0 = 10,
                              pm = 0, psd = 1/sqrt(2), dpm = 0, dpsd = 0,
                              type = "two.sample", bftype = "normal",
                              target = "H0", minN = 20, by = 20,
                              nrange = c(20, 1200), strict = FALSE)
expect_equal(fixedIncrease$n, increaseSearch$result$n,
             info = "powerbf01seq fixed-n mode should preserve searched increment schedule")
expect_equal(utils::tail(fixedIncrease$cumpH0, 1), increaseSearch$actualPower,
             tolerance = 1e-10,
             info = "powerbf01seq fixed-n mode should round-trip searched increment power")

fixedTicks <- 0L
powerbf01seq(n = 60, k1 = k1, k0 = k0, pm = pm, psd = psd,
             dpm = dpm, dpsd = dpsd, looks = 3,
             progress = function() fixedTicks <<- fixedTicks + 1L)
expect_equal(fixedTicks, 0L,
             info = "powerbf01seq fixed-n mode should not call search progress callback")

expect_error(
    nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd, pm = pm,
             psd = psd, dpm = dpm, dpsd = dpsd, looks = 3,
             nrange = c(2, 200), progress = 1),
    "progress",
    info = "sequential z progress callback should be NULL or a function"
)

expect_error(
    nbf01seq(k1 = k1, k0 = k0, power = pow, usd = usd, pm = pm,
             psd = psd, dpm = dpm, dpsd = dpsd, looks = 3,
             nrange = c(2, 200), search = "bad"),
    "should be one of",
    info = "sequential z search strategy should be validated"
)

expect_error(
    powerbf01seq(n = 60, k1 = k1, k0 = k0, pm = pm, psd = psd,
                 dpm = dpm, dpsd = dpsd, looks = 3, progress = 1),
    "progress",
    info = "powerbf01seq should validate progress in fixed-n mode"
)

expect_error(
    powerbf01seq(n = 60, k1 = k1, k0 = k0, pm = pm, psd = psd,
                 dpm = dpm, dpsd = dpsd, looks = 3, search = "bad"),
    "should be one of",
    info = "powerbf01seq should validate search strategy in fixed-n mode"
)
