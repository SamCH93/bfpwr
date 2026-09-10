library(tinytest)
library(bfpwr)

fixed_args <- list(k = 1/6, dpm = 0.5, dpsd = 0,
                   alternative = "greater", type = "two.sample",
                   drange = c(0, 2))
search_args <- c(fixed_args, list(nrange = c(5, 300)))
target_power <- 0.8
ratio <- 2

manual_ratio_root <- suppressWarnings(stats::uniroot(
    f = function(n) {
        do.call(ptbf01, c(fixed_args, list(
            n = n, n1 = n, n2 = ceiling(n*ratio)
        ))) - target_power
    },
    interval = search_args$nrange
)$root)

ratio_root <- suppressWarnings(do.call(ntbf01, c(search_args, list(
    power = target_power, ratio = ratio, integer = FALSE
))))

expect_equal(
    ratio_root,
    manual_ratio_root,
    tolerance = 1e-8,
    info = "ntbf01 should search group-1 n with group-2 n determined by ratio"
)

default_ratio <- suppressWarnings(do.call(ntbf01, c(search_args, list(
    power = target_power, integer = FALSE
))))
explicit_ratio_one <- suppressWarnings(do.call(ntbf01, c(search_args, list(
    power = target_power, ratio = 1, integer = FALSE
))))
expect_equal(
    explicit_ratio_one,
    default_ratio,
    tolerance = 1e-10,
    info = "ntbf01 ratio = 1 should preserve the default equal-allocation search"
)

low_nrange_root <- suppressWarnings(do.call(ntbf01, c(fixed_args, list(
    power = target_power, ratio = 1, integer = FALSE, nrange = c(2, 300)
))))
expect_equal(
    low_nrange_root,
    default_ratio,
    tolerance = 1e-8,
    info = "ntbf01 should skip low-n NaN endpoints for valid numeric drange searches"
)

power_search <- suppressWarnings(do.call(powertbf01, c(search_args, list(
    power = target_power, ratio = ratio
))))
expect_equal(
    power_search$n,
    ratio_root,
    tolerance = 1e-8,
    info = "powertbf01 should forward ratio and drange to ntbf01 in sample-size mode"
)
expect_equal(power_search$ratio, ratio,
             info = "powertbf01 should store the fixed-design allocation ratio")
expect_equal(power_search$drange, fixed_args$drange,
             info = "powertbf01 should store the fixed-design drange")

fixed_n <- 40
power_at_n <- suppressWarnings(do.call(powertbf01, c(search_args, list(
    n = fixed_n, ratio = ratio
))))
manual_power_at_n <- suppressWarnings(do.call(ptbf01, c(fixed_args, list(
    n = fixed_n, n1 = fixed_n, n2 = ceiling(fixed_n*ratio)
))))
expect_equal(
    power_at_n$power,
    manual_power_at_n,
    tolerance = 1e-10,
    info = "powertbf01 should forward ratio and drange to ptbf01 in fixed-n mode"
)
expect_error(
    powertbf01(n = 2, k = 1/6, dpm = 0.5, dpsd = 0,
               alternative = "greater", ratio = 0.1),
    "group-2 sample size",
    info = "powertbf01 should reject fixed-n allocations with n2 <= 1 clearly"
)

expect_error(
    powertbf01(n = 20, power = 0.8, k = 1/6, dpm = 0.5, dpsd = 0,
               alternative = "greater"),
    "exactly one of 'n' and 'power'",
    info = "powertbf01 should require exactly one of n and power"
)

small_ratio_args <- list(k = 1/6, dpm = 0.5, dpsd = 0,
                         alternative = "greater", type = "two.sample",
                         nrange = c(2, 200))
small_ratio_root <- suppressWarnings(do.call(ntbf01, c(small_ratio_args, list(
    power = 0.3, ratio = 0.2, integer = FALSE
))))
expect_true(
    is.finite(small_ratio_root) && ceiling(small_ratio_root*0.2) > 1,
    info = "ntbf01 should lift the search lower bound when small ratios would make n2 invalid"
)
