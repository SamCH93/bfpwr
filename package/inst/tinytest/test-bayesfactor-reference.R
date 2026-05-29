library(tinytest)
library(bfpwr)

source("helper-extended-tests.R", local = TRUE)
if (!bfpwr_run_extended_tests()) {
    exit_file(bfpwr_extended_skip_message(
        "BayesFactor reference checks are extended"
    ))
}

if (!requireNamespace("BayesFactor", quietly = TRUE)) {
    exit_file("BayesFactor is not installed")
}

## BayesFactor::ttest.tstat() returns log(BF10).  tbf01(..., log = TRUE)
## returns log(BF01), so the signs should be opposite for the default JZS
## two-sided t-test.
bf_cases <- data.frame(
    t = c(0.69, 3.20, 2.24, -0.90, 7.792904),
    n1 = c(100, 100, 80, 53, 500),
    n2 = c(0, 0, 0, 57, 500)
)

bf_reference <- vapply(seq_len(nrow(bf_cases)), function(i) {
    BayesFactor::ttest.tstat(t = bf_cases$t[i],
                             n1 = bf_cases$n1[i],
                             n2 = bf_cases$n2[i],
                             rscale = "medium")$bf
}, numeric(1))

bfpwr_reference <- vapply(seq_len(nrow(bf_cases)), function(i) {
    if (bf_cases$n2[i] == 0) {
        tbf01(t = bf_cases$t[i], n = bf_cases$n1[i],
              pscale = 1/sqrt(2), pdf = 1, type = "one.sample",
              alternative = "two.sided", log = TRUE)
    } else {
        tbf01(t = bf_cases$t[i], n1 = bf_cases$n1[i], n2 = bf_cases$n2[i],
              pscale = 1/sqrt(2), pdf = 1, type = "two.sample",
              alternative = "two.sided", log = TRUE)
    }
}, numeric(1))

expect_equal(
    bfpwr_reference,
    -bf_reference,
    tolerance = 1e-7,
    info = "two-sided JZS t BF agrees with BayesFactor on the log scale"
)
