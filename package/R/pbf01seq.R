#' @title Sequential Bayes Factor Design
#'
#' @description Computes cumulative probabilities of observing \eqn{z}-test
#'     Bayes factors that provide evidence for the null hypothesis
#'     \eqn{H_0}{H0}, the alternative hypothesis \eqn{H_1}{H1}, or remain
#'     inconclusive in a sequential design. Optionally, also computes the
#'     expected sample size.
#'
#' @param k1 Bayes factor threshold in favor of \eqn{H_1}{H1} (i.e.,
#'     \eqn{\text{BF}_{01} \leq \code{k1} < 1}{BF01 < \code{k1} < 1} implies
#'     evidence for \eqn{H_1})
#' @param k0 Bayes factor threshold in favor of \eqn{H_0}{H0} (i.e.,
#'     \eqn{\text{BF}_{01} \geq \code{k0} > 1}{BF01 > \code{k0} > 1} implies
#'     evidence for \eqn{H_1})
#' @param se Numeric vector of standard errors for each sequential stage
#' @param n Optional numeric vector of sample sizes corresponding to \code{se}.
#'     If supplied, the expected sample size is computed
#' @param pm Analysis prior mean. Not taken into account for \code{type =
#'     "moment"}
#' @param psd Analysis prior standard deviation (\code{type = "moment"} and
#'     \code{type = "directional"}) or scale (\code{type = "moment"})
#' @param dpm Mean of the normal design prior
#' @param dpsd Standard deviation of the normal design prior. Set \code{dpsd =
#'     0} to obtain a point prior at \code{dpm}
#'
#' @param type Character string. One of \itemize{
#' \item \code{"normal"}
#'     (default): point null vs. normal alternative (set \code{psd = 0} to
#'     obtain a point alternative) \item \code{"directional"}: directional null
#'     vs. directional alternative with a marginal normal prior \item
#'     \code{"moment"}: point null vs. normal moment alternative which is
#'     centered around 0
#' }
#'
#' @param strict Logical. If \code{TRUE} and there are more than two critical
#'     values per stage, integrate over all possible region combinations (slow
#'     but exact). If \code{FALSE}, only integrates over the main regions where
#'     the sign of the z-statistics does not change across stages (faster,
#'     recommended when many interim analyses, e.g., more than 10, are
#'     performed). Defaults to \code{TRUE}
#' @param ... Additional arguments passed to \code{mvtnorm::lpmvnorm}
#'
#' @return An object of class \code{"bfseqdesign"}, which is a list containing
#'     the input arguments, the critical z-values, the expected sample size, the
#'     cumulative probabilities of stopping for \eqn{H_1}{H1} and \eqn{H_0}{H0}
#'     by each stage, and the cumulative probabilities of remaining inconclusive
#'     by each stage.
#'
#' @details The function constructs per-stage integration regions for cumulative
#'     z-statistics based on the Bayes factor thresholds \code{k1} and
#'     \code{k0}, then computes the probability of these regions under a
#'     predictive distribution defined by \code{se} and the normal design prior
#'     with \code{dpm} and \code{dpsd}. Integration is performed via
#'     \code{mvtnorm::lpmvnorm}.
#'
#' @examples
#' n <- seq(50, 200, 50) # sample size per stage
#' se <- sqrt(2/n) # standard errors per stage
#' res <- pbf01seq(k1 = 1/10, k0 = 3, se = se, n = n, pm = 0, psd = 1,
#'                 dpm = 0.5, dpsd = 0.05, type = "normal")
#' res # print summary
#' plot(res) # plot summary
#' res$cumpH1 # cumulative probability to stop for H1 by each stage
#' res$cumpH0 # cumulative probability to stop for H0 by each stage
#' res$EN # expected sample size
#'
#' @author Samuel Pawel
#'
#' @export
pbf01seq <- function(k1, k0 = 1/k1, se, n = NULL, pm = NULL, psd, dpm = pm,
                     dpsd = psd, type = c("normal", "directional", "moment"),
                     strict = TRUE, ...) {

    ## input checks
    stopifnot(
        length(k1) == 1,
        is.numeric(k1),
        is.finite(k1),
        k1 <= 1,

        length(k0) == 1,
        is.numeric(k0),
        is.finite(k0),
        k0 >= 1,

        length(se) >= 1,
        is.numeric(se),
        all(is.finite(se)),
        all(se > 0)
    )
    if (!is.null(n)) {
        stopifnot(
            is.numeric(n),
            length(n) == length(se),
            all(is.finite(n)),
            all(n >= 1)
        )
    }
    type <- match.arg(type)
    if (type != "moment") {
        stopifnot(
            length(pm) == 1,
            is.numeric(pm),
            is.finite(pm)
        )
    }
    stopifnot(

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        psd >= 0,

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        dpsd >= 0,

        length(strict) == 1,
        is.logical(strict),
        !is.na(strict)

    )
    if (type != "normal") {
        stopifnot(psd > 0)
    }

    ## get marginal mean and covariance matrix
    pars <- predpars(se = se, dpm = dpm, dpsd = dpsd)
    mean <- pars$mean
    sigma <- pars$sigma


    ## get integration regions based on BFs with one critical value
    if ((type == "normal" & psd == 0) | type == "directional") {
        ## get region where evidence for H1 in each stage
        zk0 <- zcrit(k = k0, se = se, mu = pm, tau = psd, type = type)
        zk1 <- zcrit(k = k1, se = se, mu = pm, tau = psd, type = type)
        intregions <- genregions1(zcrit0 = zk0, zcrit1 = zk1)
    } else {
        ## get integration regions based on BFs with two critical values
        zk0 <- sapply(X = se, FUN = function(sei) {
            zcrit(k = k0, se = sei, mu = pm, tau = psd, type = type)
        })
        zk1 <- sapply(X = se, FUN = function(sei) {
            zcrit(k = k1, se = sei, mu = pm, tau = psd, type = type)
        })
        intregions <- genregions2(zcrit0 = zk0, zcrit1 = zk1, strict = strict)
    }

    ## compute stage-wise stopping probabilities
    pH1 <- intstages(intregions = intregions$H1, mean = mean, sigma = sigma,
                     ...)
    pH0 <- intstages(intregions = intregions$H0, mean = mean, sigma = sigma,
                     ...)

    ## compute cumulate stopping probabilities
    cumpH1 <- cumsum(pH1)
    cumpH0 <- cumsum(pH0)
    cumpInc <- 1 - cumpH1 - cumpH0 # inconclusive evidence

    ## compute expected sample size and variance of sample size
    if (!is.null(n)) {
    EN <- sum((pH1 + pH0)*n) + # stopping evidence for H0/H1 in stage n
        (1 - sum(pH1 + pH0))*max(n) # no evidence until last stage
    EN2 <- sum((pH1 + pH0)*n^2) +
        (1 - sum(pH1 + pH0))*max(n^2)
    VarN <- EN2 - EN^2
    } else {
        EN <- NA
        VarN <- NA
    }

    ## put everything together
    out <- structure(list("k1" = k1, "k0" = k0, "se" = se, "n" = n, "pm" = pm,
                          "psd" = psd, "dpm" = dpm, "dpsd" = dpsd,
                          "type" = type, "strict" = strict, "test" = "z",
                          "zk1" = zk1, "zk0" = zk0, "EN" = EN, "VarN" = VarN,
                          "cumpH1" = cumpH1, "cumpH0" = cumpH0,
                          "cumpInc" = cumpInc),
                     class = "bfseqdesign")
    return(out)
}

## ## compare to simulation-based probabilities
## set.seed(142)
## n <- seq(10, 50, 5) # sample size per stage
## se <- sqrt(2/n) # standard errors per stage
## dpm <- 0.3
## dpsd <- 0.1
## pm <- 0.3
## psd <- 0.5
## type <- "normal"
## k0 <- 6
## k1 <- 1/30
## nsim <- 10000
## results <- replicate(n = nsim, expr = {
##     smd <- rnorm(n = 1, mean = dpm, sd = dpsd)
##     y1 <- rnorm(n = max(n), mean = 0, sd = 1)
##     y2 <- rnorm(n = max(n), mean = smd, sd = 1)
##     smd <- sapply(seq_along(n), FUN = function(i) {
##         (mean(y2[1:n[i]]) - mean(y1[1:n[i]]))
##     })
##     se <- sapply(seq_along(n), FUN = function(i) {
##         sqrt(2/n[i])
##     })
##     bf <- sapply(seq_along(n), FUN = function(i) {
##         if (type == "normal") {
##             bf01(estimate = smd[i], se = se[i], null = 0, pm = pm, psd = psd)
##         } else if (type == "directional") {
##             dirbf01(estimate = smd[i], se = se[i], null = 0, pm = pm, psd = psd)
##         } else {
##             nmbf01(estimate = smd[i], se = se[i], null = 0, psd = psd)
##         }
##     })
##     result <- "inconclusive"
##     nfinal <- 0
##     for (i in seq_along(bf)) {
##         nfinal <- n[i]
##         if (bf[i] >= k0) {
##             result <- "H0"
##             break
##         }
##         if (bf[i] <= k1) {
##         ## if ((bf[i] < k1) & t[i] > 0) {
##             result <- "H1"
##             break
##         }
##     }
##     list("result" = result, "nfinal" = nfinal)
## }, simplify = FALSE)
## decisions <- sapply(results, function(x) x$result)
## nfinal <- sapply(results, function(x) x$nfinal)

## pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm, psd = psd, dpm = dpm,
##          dpsd = dpsd, type = type, strict = TRUE)
## mean(decisions == "H1")
## mean(decisions == "H0")
## mean(decisions == "inconclusive")
## mean(nfinal)
## var(nfinal)

## ## checks: sequential with one stage should give the same as the fixed N functions
## ## TODO implement as real tests
## k1 <- 1/10
## k0 <- 3
## pm <- 0.2
## psd <- 1
## dpm <- -0.2
## dpsd <- 0.05
## n <- 50
## usd <- sqrt(2)
## se <- usd/sqrt(50)

## ## normal alternative
## pbf01(k = k1, n = n, usd = usd, pm = pm, psd = psd, dpm = dpm, dpsd = dpsd)
## pbf01(k = k0, n = n, usd = usd, pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, lower.tail = FALSE)
## pbf01seq(k1 = k1, k0 = k0, se = se, pm = pm, psd = psd, dpm = dpm, dpsd = dpsd, type = "normal")

## ## point alternative
## pbf01(k = k1, n = n, usd = usd, pm = pm, psd = 0, dpm = dpm, dpsd = dpsd)
## pbf01(k = k0, n = n, usd = usd, pm = pm, psd = 0, dpm = dpm, dpsd = dpsd, lower.tail = FALSE)
## pbf01seq(k1 = k1, k0 = k0, se = se, pm = pm, psd = 0, dpm = dpm, dpsd = dpsd, type = "normal")

## ## normal moment alternative
## pnmbf01(k = k1, n = n, usd = usd, psd = psd, dpm = dpm, dpsd = dpsd)
## pnmbf01(k = k0, n = n, usd = usd, psd = psd, dpm = dpm, dpsd = dpsd, lower.tail = FALSE)
## pbf01seq(k1 = k1, k0 = k0, se = se, psd = psd, dpm = dpm, dpsd = dpsd, type = "moment")


#' Print method for class \code{"bfseqdesign"}
#' @method print bfseqdesign
#'
#' @param x Object of class \code{"bfseqdesign"}
#' @param digits Number of digits for formatting of numbers
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"bfseqdesign"} object
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbf01seq}
#'
#' @examples
#' n <- seq(50, 300, 50) # sample size per stage
#' se <- sqrt(2/n) # standard errors per stage
#' res <- pbf01seq(k1 = 1/10, k0 = 5, se = se, n = n, pm = 0, psd = 1,
#'                 dpm = 0.5, dpsd = 0.1, type = "normal")
#' res
#'
#' @export
print.bfseqdesign <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

    cat("\nSequential Bayes Factor Design\n")
    cat("--------------------------------\n")

    ## Hypotheses
    if (x$test == "t") {
        null <- " ="
        if (x$alternative == "two.sided") {
            alt <- "!="
        } else if (x$alternative == "greater") {
            alt <- " >"
        } else {
            alt <- " <"
        }
        if (x$type == "one.sample") {
            parlong <- "SM (stand. mean)"
            par <- "SM"
        } else {
            parlong <- "SMD (stand. mean diff.)"
            par <- "SMD"
        }

    } else {
        par <- parlong <- "parameter"
        if (x$type %in% c("normal", "moment")) {
            null <- " ="
            alt <- "!="
        }
        if (x$type == "directional") {
            null <- " <"
            alt <- " >"
        }
    }
    cat(paste0("H0:               ", parlong,  " ", null, " 0\n"))
    cat(paste0("H1:               ", parlong, " ", alt, " 0\n"))

    ## Analysis prior
    if (x$test == "t") {
        aprior <- paste0(par,
                         "|H1 ~ t(location = ", round(x$plocation, digits = digits),
                         ", scale = ",round(x$pscale, digits = digits),
                         ", df = ", round(x$pdf, digits = digits), ")")
        if (x$alternative == "greater") {
            aprior <- paste0(aprior, "_+")
        }
        if (x$alternative == "less") {
            aprior <- paste0(aprior, "_-")
        }
    } else if (x$type %in% c("normal", "directional")) {
        if (x$psd == 0) {
            aprior <- paste0("parameter = ", round(x$pm, digits = digits))
        } else {
            if (x$type == "normal") {
                parameter <- "parameter|H1" # prior is conditional on H1
            } else {
                parameter <- "parameter" # prior is marginal
            }
            aprior <- paste0(parameter,
                             " ~ N(mean = ", round(x$pm, digits = digits),
                             ", sd = ", round(x$psd, digits = digits), ")")
        }
    } else {
        aprior <- paste0("parameter|H1 ~ NM(location = 0, scale = ",
                         round(x$psd, digits = digits), ")")
    }
    cat(paste0("Analysis prior:   ", aprior, "\n"))

    ## Design prior
    if (x$dpsd == 0) {
            dprior <- paste0(par, " = ", round(x$dpm, digits = digits))
        } else {
            dprior <- paste0(par, " ~ N(mean = ", round(x$dpm, digits = digits),
                             ", sd = ", round(x$dpsd, digits = digits), ")")
        }
    cat(paste0("Design prior:     ", dprior, "\n"))

    ## Bayes factor thresholds
    k1char <- paste0("1/", round(1/x$k1, digits = digits))
    k0char <- as.character(round(x$k0, digits = digits))

    cat(sprintf("BF thresholds:    H1 if BF01 <= %s, H0 if BF01 >= %s\n",
                k1char, k0char))

    ## Stages and sample sizes
    if (x$test == "t") {
        m <- length(x$n1)
        cat(sprintf("Number of looks:  %d\n", m))
        if (x$type != "two.sample") {
            cat(sprintf("Sample sizes:     %s\n",
                        paste(round(x$n1, digits = digits), collapse = ", ")))
        } else {
            cat(sprintf("Sample sizes 1:   %s\n",
                        paste(round(x$n1, digits = digits), collapse = ", ")))
            cat(sprintf("Sample sizes 2:   %s\n",
                        paste(round(x$n2, digits = digits), collapse = ", ")))
        }
    } else {
        m <- length(x$se)
        cat(sprintf("Number of looks:  %d\n", m))
        if (!is.null(x$n)) {
            cat(sprintf("Sample sizes:     %s\n", paste(round(x$n, digits = digits),
                                                        collapse = ", ")))
        }
    }

    ## Stagewise results
    cat("\n\nStagewise cumulative probabilities:\n")
    tab <- data.frame(
        Stage = seq_len(m),
        `Pr(H1 stop)` = round(x$cumpH1, digits = digits),
        `Pr(H0 stop)` = round(x$cumpH0, digits = digits),
        `Pr(inconclusive)` = round(x$cumpInc, digits = digits),
        check.names = FALSE
    )
    print(tab, row.names = FALSE)

    ## Expected sample size
    if (x$test == "t") {
        if (x$type != "two.sample") {
            cat("\nExpected sample size: ",
                round(x$EN1, digits = digits), "\n", sep = "")
            cat("Standard deviation of sample size: ",
                round(sqrt(x$VarN1), digits = digits), "\n", sep = "")
        } else {
            cat("\nExpected sample size 1: ",
                round(x$EN1, digits = digits), "\n", sep = "")
            cat("Expected sample size 2: ",
                round(x$EN2, digits = digits), "\n", sep = "")
            cat("Standard deviation of sample size 1: ",
                round(sqrt(x$VarN1), digits = digits), "\n", sep = "")
            cat("Standard deviation of sample size 2: ",
                round(sqrt(x$VarN2), digits = digits), "\n", sep = "")
        }
    } else {
        if (!is.null(x$n) && !is.na(x$EN)) {
            cat("\nExpected sample size: ",
                round(x$EN, digits = digits), "\n", sep = "")
        }
        if (!is.null(x$n) && !is.na(x$VarN)) {
            cat("Standard deviation of sample size: ",
                round(sqrt(x$VarN), digits = digits), "\n", sep = "")
        }
    }
    ## Note
    cat("\nNOTE:  BF01 < 1 indicates evidence for H1 over H0\n\n")

    invisible(x)
}


#' Plot method for class \code{"bfseqdesign"}
#' @method plot bfseqdesign
#'
#' @param x Object of class \code{"power.bftest"}
#' @param plot Logical indicating whether data should be plotted. If
#'     \code{FALSE} only the data used for plotting are returned.
#' @param nullplot Logcal indicating whether a second plot with the stopping
#'     probabilities computed under the null hypothesis should also be produced.
#'     Defaults to \code{TRUE}
#' @param zplot Logcal indicating whether a plot of the critical z-values should
#'     be produced
#' @param digits Number of digits for formatting of numbers
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Plots stopping curves (if specified) and invisibly returns a list of
#'     data frames containing the data underlying the stopping curves
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbf01seq}
#'
#'
#' @examples
#' n <- seq(25, 150, 25) # sample size per stage
#' se <- sqrt(2/n) # standard errors per stage
#' res <- pbf01seq(k1 = 1/10, k0 = 3, se = se, n = n, pm = 0, psd = 1,
#'                 dpm = 0.5, dpsd = 0.1, type = "moment")
#' plot(res, nullplot = FALSE) # only under design prior
#' plot(res, nullplot = TRUE) # also plot under null hypothesis
#' plot(res, zplot = TRUE) # show critical z-values
#'
#' ## case which causes trouble
#' pm <- 0
#' psd <- 5
#' dpm <- 24
#' dpsd <- 0
#' seseq <- c(12.030479, 8.506833, 6.945800, 6.015239, 5.380194, 4.911422, 4.547094,
#'            4.253417, 4.010160, 3.804371)
#' nseq <- seq(2, 20, 2)
#' pbf01seq(k1 = 1/10, k0 = 3, se = seseq, n = nseq, pm = pm, psd = psd,
#'          dpm = dpm, dpsd = dpsd, type = "directional")
#'
#' @export
plot.bfseqdesign <- function(x, plot = TRUE, nullplot = TRUE, zplot = FALSE,
                             digits = max(3L, getOption("digits") - 3L), ...) {
    ## input checks
    stopifnot(
        length(plot) == 1,
        is.logical(plot),
        !is.na(plot),

        length(nullplot) == 1,
        is.logical(nullplot),
        !is.na(nullplot),

        length(nullplot) == 1,
        is.logical(nullplot),
        !is.na(nullplot)
    )

    ## data frame to return
    if (x$test == "t") {
        m <- length(x$n1)
        stages <- seq_len(m)
        xvar <- x$n1
        x$n <- x$n1
        if (x$type == "two.sample") {
            xlab <- "Sample size (group 1)"
        } else {
            xlab <- "Sample size"
        }
    } else {
        m <- length(x$se)
        stages <- seq_len(m)
        if (!is.null(x$n)) {
            xvar <- x$n
            xlab <- "Sample size"
        } else {
            x$n <- rep(NA, m)
            xvar <- stages
            xlab <- "Stage"
        }
    }
    if (zplot == FALSE) {
        plotDF <- data.frame(stage = stages, n = x$n, pH0 = x$cumpH0,
                             pH1 = x$cumpH1, pInc = x$cumpInc)
        if (nullplot == TRUE) {
            if (x$test == "t") {
                x0 <- ptbf01seq(k1 = x$k1, k0 = x$k0, n1 = x$n1, n2 = x$n2,
                                plocation = x$plocation, pscale = x$pscale,
                                pdf = x$pdf, dpm = 0, dpsd = 0, type = x$type,
                                alternative = x$alternative, trange = x$trange,
                                strict = x$strict)
            } else {
                x0 <- pbf01seq(k1 = x$k1, k0 = x$k0, se = x$se, pm = x$pm,
                               psd = x$psd, dpm = 0, dpsd = 0, type = x$type,
                               strict = x$strict)
            }
            plotDF0 <- data.frame(stage = stages, n = x$n, pH0 = x0$cumpH0,
                                  pH1 = x0$cumpH1, pInc = x0$cumpInc)
        }
    } else {
        zvals <- t(rbind(x$zk1, x$zk0))
        if (ncol(zvals) == 2) {
            colnames(zvals) <- c("zH1", "zH0")
        } else {
            colnames(zvals) <- c("zH1.1", "zH1.2", "zH0.1", "zH0.2")
        }
        plotDF <- data.frame(stage = stages, n = x$n, zvals)
    }
    if (plot == TRUE) {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))

        k1char <- paste0("1/", round(1/x$k1, digits = digits))
        k0char <- as.character(round(x$k0, digits = digits))

        if (zplot == TRUE) {
            graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(1, 10))
            graphics::par(mar = c(0, 0, 0, 0))
            graphics::plot.new()
            graphics::legend("center",
                             legend = c(bquote("Stop for" ~ italic(H)[1] ~
                                                   "(BF"["01"] <= .(k1char)*")  "),
                                        bquote("Stop for" ~ italic(H)[0] ~
                                                   "(BF"["01"] >= .(k0char)*")  ")),
                             lty = 1, pch = 20, lwd = 1.5, col = c(4, 2),
                             horiz = TRUE, xpd = TRUE, bty = "n", text.width = NA)

            graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(xvar, zvals[,1], type = "n", xlab = xlab,
                 ylab = bquote("Critical" ~ italic(z) * "-value"),
                 ylim = c(min(c(zvals, 0), na.rm = TRUE), max(c(zvals, 0), na.rm = TRUE)),
                 las = 1,
                 panel.first = graphics::grid(lty = 3, col = "#0000001A"))
            graphics::matlines(xvar, zvals, type = "b", pch = 20, lwd = 1.5,
                               lty = 1, cex = 1.5,
                               col = c(rep(4, ncol(zvals)/2), rep(2, ncol(zvals)/2)))
        } else {


            if (nullplot == TRUE) {
                graphics::layout(matrix(c(1, 2, 3), ncol = 1), heights = c(1, 10, 10))
            } else {
                graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(1, 10))
            }

            if (x$test != "t") {
                parameter <- "Assuming parameter"
            } else {
                if (x$type == "one.sample") {
                    parameter <- "Assuming standardized mean"
                } else {
                    parameter <- "Assuming standardized mean difference"
                }
            }

            graphics::par(mar = c(0, 0, 0, 0))
            graphics::plot.new()
            graphics::legend("center",
                             legend = c(bquote("Stop for" ~ italic(H)[1] ~
                                                   "(BF"["01"] <= .(k1char)*")  "),
                                        bquote("Inconclusive  "),
                                        bquote("Stop for" ~ italic(H)[0] ~
                                                   "(BF"["01"] >= .(k0char)*")  ")),
                             pch = 20, lwd = 1.5, lty = 1, col = c(4, 1, 2),
                             horiz = TRUE, xpd = TRUE, bty = "n",
                             text.width = NA)

            if (x$dpsd == 0) {
                title <- paste0(parameter, " = ", round(x$dpm, digits = digits))
            } else {
                title <- paste0(parameter,
                                " ~ N(mean = ", round(x$dpm, digits = digits),
                                ", sd = ", round(x$dpsd, digits = digits), ")")
            }
            graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
            plot(xvar, x$cumpH1, type = "n", xlab = xlab, ylab = "Probability",
                 ylim = c(0, 100), yaxt = "n",
                 panel.first = graphics::grid(lty = 3, col = "#0000001A"), main = title)
            graphics::matlines(xvar, cbind(x$cumpH1, x$cumpInc, x$cumpH0)*100,
                               type = "b", pch = 20, col = c(4, 1, 2), lwd = 1.5,
                               lty = 1, cex = 1.5)
            graphics::axis(side = 2, at = seq(0, 100, 20),
                           labels = paste0(seq(0, 100, 20), "%"), las = 1)

            if (nullplot == TRUE) {
                plot(xvar, x0$cumpH1, type = "n", xlab = xlab, ylab = "Probability",
                     ylim = c(0, 100), yaxt = "n",
                     panel.first = graphics::grid(lty = 3, col = "#0000001A"),
                     main = paste0(parameter, " = 0"))
                graphics::matlines(xvar, cbind(x0$cumpH1, x0$cumpInc, x0$cumpH0)*100,
                                   type = "b", pch = 20, col = c(4, 1, 2), lwd = 1.5,
                                   lty = 1, cex = 1.5)
                graphics::axis(side = 2, at = seq(0, 100, 20),
                               labels = paste0(seq(0, 100, 20), "%"), las = 1)
            }
        }

    }
    ## invisibly return data frames
    if (zplot == TRUE) {
        ret <- list("zDF" = plotDF)
    } else if (nullplot == TRUE) {
        ret <- list("pDF1" = plotDF, "pDF2" = plotDF0)
    } else {
        ret <- list("pDF1" = plotDF)
    }
    invisible(ret)
}

#' ## problematic example (-> now fixed with specified integration grid)
#' k1 <- 1/10
#' k0 <- 3
#' pm <- 0
#' psd <- 5
#' dpm <- 23.5
#' dpsd <- 0
#' nseq <- seq(2, 20, 2)
#' type <- "directional"
#' se <- c(12.03, 8.51, 6.94, 6.01, 5.38, 4.91, 4.54, 4.25, 4.01, 3.80)
#' pbf01seq(k1 = k1, k0 = k0, se = se, n = nseq, pm = pm, psd = psd, dpm = dpm,
#' dpsd = dpsd, type = type)
