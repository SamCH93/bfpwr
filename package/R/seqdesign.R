## Shared builders for sequential BF design objects
## -----------------------------------------------------------------------------

## Collect per-look standard errors and BF critical values in the shape needed
## by the region generators.
.bfseq_boundary_data <- function(bounds, oneCritical) {
    se <- vapply(bounds, `[[`, numeric(1), "se")
    if (oneCritical) {
        zk0 <- vapply(bounds, function(x) x$zk0[[1]], numeric(1))
        zk1 <- vapply(bounds, function(x) x$zk1[[1]], numeric(1))
    } else {
        zk0 <- do.call(cbind, lapply(bounds, `[[`, "zk0"))
        zk1 <- do.call(cbind, lapply(bounds, `[[`, "zk1"))
    }
    list(se = se, zk0 = zk0, zk1 = zk1)
}

## Translate BF stopping boundaries into H1/H0 integration regions for one
## stage, dispatching to the one- or two-critical-value region geometry.
.bfseq_stage_regions <- function(boundaries, oneCritical, strict,
                                 direction = NULL) {
    if (oneCritical) {
        return(.bfseq_genregions1_stage(zcrit0 = boundaries$zk0,
                                        zcrit1 = boundaries$zk1,
                                        direction = direction))
    }
    .bfseq_genregions2_stage(zcrit0 = boundaries$zk0,
                             zcrit1 = boundaries$zk1,
                             strict = strict)
}

## Integrate the predictive distribution over the H1 and H0 stopping regions
## for a single stage.
.bfseq_stage_stop_probabilities <- function(regions, se, null = 0, dpm, dpsd,
                                            dots) {
    pars <- predpars(se = se, null = null, dpm = dpm, dpsd = dpsd)
    pH1 <- do.call(.bfseq_intstage,
                   c(list(stageregions = regions$H1,
                          mean = pars$mean,
                          sigma = pars$sigma),
                     dots))
    pH0 <- do.call(.bfseq_intstage,
                   c(list(stageregions = regions$H0,
                          mean = pars$mean,
                          sigma = pars$sigma),
                     dots))
    if (!is.numeric(pH1) || length(pH1) != 1 || !is.finite(pH1) ||
        !is.numeric(pH0) || length(pH0) != 1 || !is.finite(pH0)) {
        .bfseq_candidate_invalid(
            "non-finite sequential stage probability",
            reason = "stage_probability",
            terminal = TRUE
        )
    }
    list(pH1 = pH1, pH0 = pH0)
}

## Full stage calculation from raw boundary objects: reshape boundaries,
## construct stopping regions, then integrate them.
.bfseq_stage_probabilities_from_bounds <- function(bounds, oneCritical, strict,
                                                   direction, null = 0, dpm,
                                                   dpsd, dots) {
    boundaries <- .bfseq_boundary_data(bounds = bounds,
                                       oneCritical = oneCritical)
    regions <- .bfseq_stage_regions(boundaries = boundaries,
                                    oneCritical = oneCritical,
                                    strict = strict,
                                    direction = direction)
    .bfseq_stage_stop_probabilities(regions = regions, se = boundaries$se,
                                    null = null, dpm = dpm, dpsd = dpsd,
                                    dots = dots)
}

## First two moments of the stopping sample size under the stage-wise stopping
## probabilities, with non-stoppers assigned the maximum planned sample size.
.bfseq_sample_size_moments <- function(pH1, pH0, n) {
    stopProb <- pH1 + pH0
    EN <- sum(stopProb*n) + (1 - sum(stopProb))*max(n)
    EN2 <- sum(stopProb*n^2) + (1 - sum(stopProb))*max(n^2)
    list(EN = EN, VarN = EN2 - EN^2)
}

## Build a sequential z-test design object from a look schedule. Optional
## boundary and stage callbacks let sample-size searches reuse cached work.
.bfseq_build_z_design <- function(k1, k0, se, n = NULL, null = 0, pm, psd,
                                  dpm, dpsd, type, strict, dots,
                                  getBoundary = NULL, evalStage = NULL) {
    oneCritical <- (type == "normal" && psd == 0) || type == "directional"
    integration <- .bfseq_integration_settings(dots)

    if (is.null(getBoundary)) {
        getBoundary <- function(i) {
            list(
                n = if (is.null(n)) NA_real_ else n[[i]],
                se = se[[i]],
                zk0 = zcrit(k = k0, se = se[[i]], null = null, mu = pm,
                            tau = psd, type = type),
                zk1 = zcrit(k = k1, se = se[[i]], null = null, mu = pm,
                            tau = psd, type = type)
            )
        }
    }

    bounds <- lapply(seq_along(se), getBoundary)
    boundaries <- .bfseq_boundary_data(bounds = bounds,
                                       oneCritical = oneCritical)
    stages <- lapply(seq_along(se), function(i) {
        if (!is.null(evalStage)) {
            return(evalStage(i = i, bounds = bounds[seq_len(i)]))
        }
        .bfseq_stage_probabilities_from_bounds(
            bounds = bounds[seq_len(i)], oneCritical = oneCritical,
            strict = strict, direction = NULL, null = null, dpm = dpm,
            dpsd = dpsd, dots = dots
        )
    })
    pH1 <- vapply(stages, `[[`, numeric(1), "pH1")
    pH0 <- vapply(stages, `[[`, numeric(1), "pH0")
    cumpH1 <- cumsum(pH1)
    cumpH0 <- cumsum(pH0)
    cumpInc <- 1 - cumpH1 - cumpH0

    if (is.null(n)) {
        EN <- NA_real_
        VarN <- NA_real_
    } else {
        moments <- .bfseq_sample_size_moments(pH1 = pH1, pH0 = pH0, n = n)
        EN <- moments$EN
        VarN <- moments$VarN
    }

    structure(list(
        k1 = k1, k0 = k0, se = boundaries$se, n = n, null = null, pm = pm,
        psd = psd, dpm = dpm, dpsd = dpsd, type = type,
        strict = strict, integration = integration, test = "z",
        zk1 = boundaries$zk1,
        zk0 = boundaries$zk0, EN = EN, VarN = VarN,
        cumpH1 = cumpH1, cumpH0 = cumpH0, cumpInc = cumpInc
    ), class = "bfseqdesign")
}

## Build a sequential t-test design object from precomputed t boundaries and
## derive cumulative stopping probabilities and expected sample sizes.
.bfseq_build_t_design <- function(k1, k0, bounds, dpm, dpsd, plocation,
                                  pscale, pdf, alternative, type, trange,
                                  strict, tail.eps, tail.nquad, dots,
                                  evalStage = NULL) {
    oneCritical <- alternative != "two.sided"
    integration <- .bfseq_integration_settings(dots)
    regionDirection <- if (alternative == "greater") {
        "positive"
    } else if (alternative == "less") {
        "negative"
    } else {
        NULL
    }

    boundaries <- .bfseq_boundary_data(bounds = bounds,
                                       oneCritical = oneCritical)
    n1 <- vapply(bounds, `[[`, numeric(1), "n1")
    n2 <- vapply(bounds, `[[`, numeric(1), "n2")
    stages <- lapply(seq_along(bounds), function(i) {
        if (!is.null(evalStage)) {
            return(evalStage(i = i, bounds = bounds[seq_len(i)]))
        }
        .bfseq_stage_probabilities_from_bounds(
            bounds = bounds[seq_len(i)], oneCritical = oneCritical,
            strict = strict, direction = regionDirection, dpm = dpm,
            dpsd = dpsd, dots = dots
        )
    })
    pH1 <- vapply(stages, `[[`, numeric(1), "pH1")
    pH0 <- vapply(stages, `[[`, numeric(1), "pH0")
    cumpH1 <- cumsum(pH1)
    cumpH0 <- cumsum(pH0)
    cumpInc <- 1 - cumpH1 - cumpH0

    moments1 <- .bfseq_sample_size_moments(pH1 = pH1, pH0 = pH0, n = n1)
    moments2 <- .bfseq_sample_size_moments(pH1 = pH1, pH0 = pH0, n = n2)
    structure(list(
        k1 = k1, k0 = k0, n1 = n1, n2 = n2, dpm = dpm,
        dpsd = dpsd, plocation = plocation, pscale = pscale,
        pdf = pdf, alternative = alternative, type = type,
        trange = trange, strict = strict, integration = integration, test = "t",
        tail.eps = tail.eps, tail.nquad = tail.nquad,
        zk1 = boundaries$zk1, zk0 = boundaries$zk0,
        EN1 = moments1$EN, EN2 = moments2$EN,
        VarN1 = moments1$VarN, VarN2 = moments2$VarN,
        cumpH1 = cumpH1, cumpH0 = cumpH0, cumpInc = cumpInc
    ), class = "bfseqdesign")
}
