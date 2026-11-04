library(tinytest)
library(bfpwr)

## R's existing matching of method/ngrid also accepts unambiguous aliases and
## positional controls. Metadata retains its original exact-name recording.
args <- list(k1 = 0.1, k0 = 3, se = c(0.3, 0.2), n = c(20, 40),
             null = 0, pm = 0, psd = 1, dpm = 0.2, dpsd = 0.1,
             type = "normal", strict = TRUE, alternative = "two.sided")
for (controls in list(list(meth = "pmvnorm"), list(ngr = 31),
                      list("pmvnorm"), list("lpmvnorm", 31))) {
    design <- do.call(pbf01seq, c(args, controls))
    referenceControls <- if (is.null(names(controls))) {
        if (length(controls) == 1) {
            list(method = controls[[1]])
        } else list(method = controls[[1]], ngrid = controls[[2]])
    } else if ("meth" %in% names(controls)) {
        list(method = controls$meth)
    } else list(ngrid = controls$ngr)
    reference <- do.call(pbf01seq, c(args, referenceControls))
    expect_equal(design$cumpH1, reference$cumpH1, tolerance = 1e-13)
    expect_equal(design$cumpH0, reference$cumpH0, tolerance = 1e-13)
    expect_equal(design$cumpInc, reference$cumpInc, tolerance = 1e-13)
    expect_equal(design$integration,
                 bfpwr:::.bfseq_integration_settings(controls))
}
