library(tinytest)
library(bfpwr)

## Compare the full-design evaluator with the retained single-stage path used
## by adaptive searches. These references independently prepare every prefix.
expect_stage_equivalence <- function(design) {
    isZ <- design$test == "z"
    oneCritical <- if (isZ) {
        bfpwr:::.bfseq_z_one_critical(design$type, design$psd, design$alternative)
    } else design$alternative != "two.sided"
    se <- if (isZ) design$se else if (design$type == "two.sample") {
        sqrt(1/design$n1 + 1/design$n2)
    } else 1/sqrt(design$n1)
    bounds <- lapply(seq_along(se), function(i) {
        list(se = se[i],
             zk0 = if (oneCritical) design$zk0[i] else design$zk0[, i],
             zk1 = if (oneCritical) design$zk1[i] else design$zk1[, i])
    })
    dots <- design$integration[!names(design$integration) %in% "bf.control"]
    stages <- lapply(seq_along(se), function(i) {
        bfpwr:::.bfseq_stage_probabilities_from_bounds(
            bounds = bounds[seq_len(i)], oneCritical = oneCritical,
            strict = design$strict,
            direction = bfpwr:::.bfseq_region_direction(design$alternative),
            null = if (isZ) design$null else 0, dpm = design$dpm,
            dpsd = design$dpsd, dots = dots)
    })
    pH1 <- vapply(stages, `[[`, numeric(1), "pH1")
    pH0 <- vapply(stages, `[[`, numeric(1), "pH0")
    expect_equal(design$cumpH1, cumsum(pH1), tolerance = 1e-13)
    expect_equal(design$cumpH0, cumsum(pH0), tolerance = 1e-13)
    expect_equal(design$cumpInc, 1 - cumsum(pH1) - cumsum(pH0), tolerance = 1e-13)
    for (group in if (isZ) "" else c("1", "2")) {
        n <- design[[paste0("n", group)]]
        if (!is.null(n)) {
            moments <- bfpwr:::.bfseq_sample_size_moments(pH1, pH0, n)
            expect_equal(design[[paste0("EN", group)]], moments$EN,
                         tolerance = 1e-13)
            expect_equal(design[[paste0("VarN", group)]], moments$VarN,
                         tolerance = 1e-13)
        }
    }
}

## The information ratios and design-prior covariance keep their original
## units, including a shifted null and nonuniform standard errors.
for (se in list(0.3, c(0.9, 0.6, 0.2, 0.07))) {
    for (dpsd in c(0, 0.13, 2)) {
        inf <- 1/se^2
        sigma <- matrix(NA_real_, length(se), length(se))
        for (i in seq_along(se)) {
            for (j in seq_along(se)) {
                sigma[i, j] <- sqrt(min(inf[i], inf[j])/max(inf[i], inf[j])) +
                    dpsd^2*sqrt(inf[i])*sqrt(inf[j])
            }
        }
        pars <- bfpwr:::predpars(se, null = 0.4, dpm = -0.1, dpsd = dpsd)
        expect_equal(pars$mean, (-0.1 - 0.4)/se, tolerance = 1e-15)
        expect_equal(pars$sigma, sigma, tolerance = 1e-15)
    }
}

local({
    oldKind <- RNGkind()
    hadSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (hadSeed) oldSeed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(oldKind))
        if (hadSeed) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    })
    zbase <- list(k1 = 0.1, k0 = 3, se = 1/sqrt(c(20, 50, 90)),
                  n = c(20, 50, 90), null = 0.2, pm = 0.25, psd = 0.5,
                  dpm = 0.35, dpsd = 0.1, strict = TRUE, ngrid = 37)
    variants <- list(
        list(type = "normal", alternative = "two.sided"),
        list(type = "normal", alternative = "greater"),
        list(type = "normal", alternative = "less"),
        list(type = "normal", psd = 0, pm = 0.4, alternative = "greater"),
        list(type = "directional"),
        list(type = "moment", alternative = "two.sided"),
        list(type = "moment", alternative = "greater"),
        list(type = "moment", alternative = "less"))
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG", "Wichmann-Hill")) {
        RNGkind(kind)
        set.seed(3901)
        seed <- get(".Random.seed", envir = .GlobalEnv)
        kinds <- RNGkind()
        for (variant in variants) {
            design <- do.call(pbf01seq, modifyList(zbase, variant))
            expect_stage_equivalence(design)
        }
        tbase <- list(k1 = 0.1, k0 = 2, n1 = c(20, 35, 60),
                      n2 = c(25, 40, 70), dpm = 0.3, dpsd = 0.1,
                      ngrid = 37, bf.control = list(rel.tol = 1e-7))
        for (alternative in c("two.sided", "greater", "less")) {
            expect_stage_equivalence(do.call(ptbf01seq,
                c(tbase, list(alternative = alternative))))
        }
        expect_stage_equivalence(ptbf01seq(
            k1 = 0.1, k0 = 3, n = c(20, 35, 60), type = "paired",
            alternative = "greater", dpm = 0.3, dpsd = 0, ngrid = 37))
        ## pmvnorm keeps its algorithm and numerical controls, and restores
        ## RNG state under the caller's generator as before.
        pmv <- list(method = "pmvnorm",
                    algorithm = mvtnorm::GenzBretz(maxpts = 1234,
                        abseps = 1e-7, releps = 1e-7))
        expect_stage_equivalence(do.call(pbf01seq, c(zbase, pmv)))
        expect_stage_equivalence(do.call(ptbf01seq,
            c(tbase, list(alternative = "greater"), pmv)))
        expect_identical(RNGkind(), kinds)
        expect_identical(get(".Random.seed", envir = .GlobalEnv), seed)
        expect_error(do.call(pbf01seq, c(zbase, list(unknown.control = 1))))
        expect_identical(RNGkind(), kinds)
        expect_identical(get(".Random.seed", envir = .GlobalEnv), seed)
    }
    kinds <- RNGkind()
    rm(".Random.seed", envir = .GlobalEnv)
    do.call(pbf01seq, zbase)
    ## withr initializes an absent seed in the existing single-stage path too.
    expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    expect_identical(RNGkind(), kinds)
    rm(".Random.seed", envir = .GlobalEnv)
    do.call(ptbf01seq, tbase)
    expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    expect_identical(RNGkind(), kinds)
})

## The cached path retains integration-control validation, empty-region
## handling and stage-probability validation.
for (controls in list(list(ngrid = 0), list(ngrid = 1.5),
                      list(method = "unknown"), list(ngrid = NA_real_))) {
    expect_error(do.call(pbf01seq, c(list(k1 = 0.1, se = c(0.3, 0.2),
        pm = 0, psd = 1), controls)))
    expect_error(do.call(ptbf01seq, c(list(k1 = 0.1, n = c(20, 40)), controls)))
}
local({
    pars <- bfpwr:::predpars(c(0.3, 0.2), dpm = 0, dpsd = 0.1)
    evaluator <- bfpwr:::.bfseq_region_evaluator(pars, ngrid = 17)
    stage <- function(regions, evalRegions = evaluator) {
        bfpwr:::.bfseq_stage_stop_probabilities(
            regions, se = c(0.3, 0.2), dpm = 0, dpsd = 0.1,
            dots = list(ngrid = 17), evalRegions = evalRegions)
    }
    regions <- list(H1 = list(matrix(c(-Inf, Inf, 0, Inf), nrow = 2)),
                    H0 = list(matrix(NaN, nrow = 2, ncol = 2)))
    expect_equal(stage(regions)$pH0, 0)
    regions$H1[[1]][1, 1] <- NA_real_
    expect_error(stage(regions), "NA bounds")
    expect_error(stage(regions, function(...) c(H1 = Inf, H0 = 0)),
                 "non-finite sequential stage probability")
})
