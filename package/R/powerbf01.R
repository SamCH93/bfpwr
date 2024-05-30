#' @title Power and sample size calculations for Bayes factor under normality
#'
#' @description Compute probability that Bayes factor under normality is smaller
#'     than a specified threshold (the power), or determine sample size to
#'     obtain a target power
#'
#' @note An error message will be displayed in case that the specified target
#'     power is not achievable under the specified analysis and design priors.
#'
#' @param k Bayes factor threshold. Defaults to \code{1/10}, Jeffreys' threshold
#'     for 'strong evidence' against the null hypothesis
#' @param n Sample size. Has to be \code{NULL} if \code{power} is specified.
#'     Defaults to \code{NULL}
#' @param power Target power. Has to be \code{NULL} if \code{n} is specified.
#'     Defaults to \code{NULL}
#' @param sd Standard deviation of one observation (for \code{type =
#'     "two.sample"} or \code{type = "one.sample"}) or of one difference within
#'     a pair of observations (\code{type = "paired"}). Is assumed to be known
#' @param null Mean difference under the point null hypothesis. Defaults to
#'     \code{0}
#' @param pm Mean of the normal prior assigned to the mean difference under the
#'     alternative in the analysis
#' @param psd Standard deviation of the normal prior assigned to the mean
#'     difference under the alternative in the analysis. Set to \code{0} to
#'     obtain a point prior at the prior mean
#' @param dpm Mean of the normal design prior assigned to the mean difference.
#'     Defaults to the same value as the analysis prior \code{pm}
#' @param dpsd Standard deviation of the normal design prior assigned to the
#'     mean difference. Defaults to the same value as the analysis prior
#'     \code{psd}
#' @param nrange Sample size search range (only taken into account when \code{n}
#'     is \code{NULL}). Defaults to \code{c(1, 10^5)}
#' @param type The type of test. One of \code{"two.sample"},
#'     \code{"one.sample"}, \code{"paired"}. Defaults to \code{"two.sample"}
#'
#' @return Object of class \code{"power.bftest"}, a list of the arguments
#'     (including the computed one) augmented with \code{method} and \code{note}
#'     elements
#'
#' @author Samuel Pawel
#'
#' @seealso \link{plot.power.bftest}, \link{nbf01}, \link{pbf01}
#'
#' @examples
#' ## determine power
#' powerbf01(n = 100, pm = 0, psd = 1, dpm = 0.5, dpsd = 0)
#'
#' ## determine sample size
#' powerbf01(power = 0.99, pm = 0, psd = 1, dpm = 0.5, dpsd = 0)
#'
#' @export

powerbf01 <- function(n = NULL, power = NULL, k = 1/10, sd = 1, null = 0, pm, psd,
                      dpm = pm, dpsd = psd, nrange = c(1, 10^5),
                      type = c("two.sample", "one.sample", "paired")) {
    ## input checks
    if (is.null(n) && is.null(power)) {
        stop("exactly one of 'n' and 'power' must be NULL")
    }
    if (is.null(n)) {
        stopifnot(
            length(power) == 1,
            is.numeric(power),
            is.finite(power),
            0 < power, power < 1
        )
    } else {
        stopifnot(
            length(n) == 1,
            is.numeric(n),
            is.finite(n),
            0 < n
        )
    }
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(sd) == 1,
        is.numeric(sd),
        is.finite(sd),
        0 < sd,

        length(null) == 1,
        is.numeric(null),
        is.finite(null),

        length(pm) == 1,
        is.numeric(pm),
        is.finite(pm),

        length(psd) == 1,
        is.numeric(psd),
        is.finite(psd),
        0 <= psd,

        length(dpm) == 1,
        is.numeric(dpm),
        is.finite(dpm),

        length(dpsd) == 1,
        is.numeric(dpsd),
        is.finite(dpsd),
        0 <= dpsd,

        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1]
    )
    type <- match.arg(type)

    ## determine unit variance
    if (type == "two.sample") {
        uv <- 2*sd^2
    } else {
        uv <- sd^2
    }

    ## determine sample size
    if (is.null(n)) {
        n <- nbf01(k = k, power = power, sd = sqrt(uv), null = null, pm = pm,
                   psd = psd, dpm = dpm, dpsd = dpsd, nrange = nrange,
                   integer = FALSE, analytical = TRUE)
    } else {
        ## determine power
        power <- pbf01(k = k, n = n, sd = sqrt(uv), null = null, pm = pm,
                       psd = psd, dpm = dpm, dpsd = dpsd, lower.tail = TRUE)
    }

    ## return object
    structure(list(n = n, power = power, sd = sd, null = null, pm = pm,
                   psd = psd, dpm = dpm, dpsd = dpsd, k = k, nrange = nrange,
                   type = type),
              class = "power.bftest")

}

#' Print method for class \code{"power.bftest"}
#' @method print power.bftest
#'
#' @param x Object of class \code{"power.bftest"}
#' @param digits Number of digits for formatting of numbers
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"power.bftest"} object
#'
#' @note Function adapted from \code{stats:::print.power.htest} written by Peter
#'     Dalgaard
#'
#' @author Samuel Pawel
#'
#' @seealso \link{powerbf01}
#'
#' @examples
#' powerbf01(power = 0.95, pm = 0, psd = 1, dpm = 0.5, dpsd = 0)
#' powerbf01(power = 0.95, pm = 0, psd = 1, dpm = 0.5, dpsd = 0, type = "one.sample")
#' powerbf01(power = 0.95, pm = 0, psd = 1, dpm = 0.5, dpsd = 0, type = "paired")
#' powerbf01(power = 0.95, pm = 1, psd = 0, dpm = 0.8, dpsd = 0, type = "paired")
#'
#' @export
print.power.bftest <- function(x, digits = getOption("digits"), ...) {

    method <- c("Bayes factor power calculation")
    note <- "BF oriented in favor of H0 (BF < 1 indicates evidence for H1 over H0)"
    if (x$type == "paired") {
        note <- paste(note,
                      "n is number of *pairs*",
                      "sd is standard deviation of *differences* between pairs",
                      sep = "\n      ")
        method <- paste("Paired", method)
    } else if (x$type == "one.sample") {
        note <- paste(note,
                      "n is number of *samples*",
                      "sd is standard deviation of one sample",
                      sep = "\n      ")
        method <- paste("One-sample", method)
    } else {
        note <- paste(note,
                      "n is number of *samples per group*",
                      "sd is standard deviation of one sample (assumed equal in both groups)",
                      sep = "\n      ")
        method <- paste("Two-sample", method)
    }
    cat("\n    ", method, "\n\n")

     ## formatting BF threshold k
    if (x$k < 1) {
        x$k <- paste0("1/", format(x = 1/x$k, digits = digits))
    }
    printx <- x[c("n", "power", "sd", "null", "pm", "psd", "dpm", "dpsd", "k")]
    names(printx) <- c("n", "power", "sd", "null",
                       "analysis prior mean", "analysis prior sd", "design prior mean",
                       "design prior sd", "BF threshold k")

    cat(paste(paste("     ", format(names(printx), width = 10L, justify = "right")),
              format(printx, digits = digits), sep = " = "), sep = "\n")
    cat("\n", "NOTE: ", note, "\n\n", sep = "")
    invisible(x)
}


#' Plot method for class \code{"power.bftest"}
#' @method plot power.bftest
#'
#' @param x Object of class \code{"power.bftest"}
#' @param nlim Range of samples sizes over which the power should be computed.
#'     Defaults to \code{c(1, 500)}
#' @param plot Logical indicating whether data should be plotted. If
#'     \code{FALSE} only the data used for plotting are returned.
#' @param nullplot Logcal indicating whether a second plot with the power in
#'     favor of the null (using a Bayes factor threshold of 1/k) should be
#'     created. Defaults to \code{TRUE}
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Plots power curves (if specified) and invisibly returns a list of
#'     data frames containing the data underlying the power curves
#'
#' @author Samuel Pawel
#'
#' @seealso \link{powerbf01}
#'
#' @examples
#' ssd1 <- powerbf01(k = 1/6, power = 0.95, pm = 0, psd = 1/sqrt(2), dpm = 0.5, dpsd = 0)
#' plot(ssd1, nlim = c(1, 8000))
#'
#' power1 <- powerbf01(k = 1/2, n = 120, pm = 0, psd = 1/sqrt(2), dpm = 0.5, dpsd = 0)
#' plot(power1, nlim = c(1, 1000))
#'
#' @export
plot.power.bftest <- function(x, nlim = c(1, 500), plot = TRUE, nullplot = TRUE,
                              ...) {
    ## input checks
    stopifnot(
        length(nlim) == 2,
        all(is.numeric(nlim)),
        all(is.finite(nlim)),
        nlim[2] > nlim[1],
        nlim[1] >= 1,

        length(plot) == 1,
        is.logical(plot),
        !is.na(plot),

        length(nullplot) == 1,
        is.logical(nullplot),
        !is.na(nullplot)
    )


    ## determine unit standard deviation
    if (x$type == "two.sample") {
        usd <- sqrt(2)*x$sd
    } else {
        usd <- x$sd
    }
    ## compute power curves
    nseq <- seq(from = nlim[1], to = nlim[2], by = 1)
    pow <- pbf01(k = x$k, n = nseq, sd = usd, null = x$null, pm = x$pm,
                 psd = x$psd, dpm = x$dpm, dpsd = x$dpsd)
    powNull <- pbf01(k = x$k, n = nseq, sd = usd, null = x$null, pm = x$pm,
                     psd = x$psd, dpm = x$null, dpsd = x$null)
    powDF <- rbind(data.frame(n = nseq, power = pow, prior = "Design prior"),
                   data.frame(n = nseq, power = powNull, prior = "Null"))
    if (nullplot) {
        ## compute power in favor of H0
        powH0 <- pbf01(k = 1/x$k, n = nseq, sd = usd, null = x$null, pm = x$pm,
                       psd = x$psd, dpm = x$dpm, dpsd = x$dpsd,
                       lower.tail = FALSE)
        powNullH0 <- pbf01(k = 1/x$k, n = nseq, sd = usd, null = x$null,
                           pm = x$pm, psd = x$psd, dpm = x$null, dpsd = x$null,
                           lower.tail = FALSE)
        powDFH0 <- rbind(data.frame(n = nseq, power = pow, prior = "Design prior"),
                         data.frame(n = nseq, power = powNull, prior = "Null"))
        ## compute n to achieve target power in favor of the null
        nH0 <- nbf01(k = 1/x$k, power = x$power, sd = usd, null = x$null,
                     pm = x$pm, psd = x$psd, dpm = x$null, dpsd = x$null,
                     lower.tail = FALSE)

    }

    ## plot power curves if specified
    if (plot == TRUE) {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))
        graphics::par(mar = c(4, 5, 2.5, 3))
        if (nullplot) {
            graphics::par(mfrow = c(2, 1))
        }
        if (x$type == "one.sample") {
            xlab <-  bquote("Sample size" ~ italic(n))
        } else if (x$type == "two.sample") {
            xlab <-  bquote("Sample size per group" ~ italic(n))
        } else {
            xlab <-  bquote("Number of pairs" ~ italic(n))
        }
        if (x$k < 1) {
            kformat <- paste0("1/", format(x = 1/x$k, digits = getOption("digits")))
            if (nullplot) {
                kformatH0 <- paste0(format(x = 1/x$k, digits = getOption("digits")))
            }
        } else {
            kformat <- format(x = x$k, digits = getOption("digits"))
            if (nullplot) {
                kformatH0 <- paste0("1/", format(x = 1/x$k, digits = getOption("digits")))
            }
        }
        plot(nseq, pow*100, xlab = xlab,
             ylab = bquote("Pr(BF"["01"] < .(kformat) * " )"), type = "s",
             ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = 4,
             panel.first = graphics::grid(lty = 3, col = "#0000001A"))
        graphics::lines(nseq, powNull*100, type = "l", col = 2, lwd = 1.5)
        graphics::axis(side = 2, at = seq(0, 100, 20),
                       labels = paste0(seq(0, 100, 20), "%"), las = 1)
        graphics::axis(side = 4, at = x$power*100,
                       labels = paste0(round(x$power*100, 1), "%"), las = 1,
                       col = "#00000033", cex.axis = 0.8)
        graphics::abline(h = x$power*100, col = "#00000033")
        if (is.finite(x$n)) {
            graphics::axis(side = 3, at = ceiling(x$n), col = "#00000033", cex.axis = 0.8)
            graphics::abline(v = ceiling(x$n), col = "#00000033")
        }
        graphics::legend("right", title = "Data distribution",
                         legend = c("Design prior", "Null hypothesis"), lty = 1,
                         lwd = 1.5, col = c(4, 2), bg = "white", cex = 0.8)
        if (nullplot) {
            plot(nseq, powH0*100, xlab = xlab,
                 ylab = bquote("Pr(BF"["01"] > .(kformatH0) * " )"), type = "s",
                 ylim = c(0, 100), lwd = 1.5, yaxt = "n", col = 4,
                 panel.first = graphics::grid(lty = 3, col = "#0000001A"))
            graphics::lines(nseq, powNullH0*100, type = "l", col = 2, lwd = 1.5)
            graphics::axis(side = 2, at = seq(0, 100, 20),
                           labels = paste0(seq(0, 100, 20), "%"), las = 1)
            graphics::axis(side = 4, at = x$power*100,
                           labels = paste0(round(x$power*100, 1), "%"), las = 1,
                           col = "#00000033", cex.axis = 0.8)
            graphics::abline(h = x$power*100, col = "#00000033")
            if (is.finite(nH0)) {
                graphics::axis(side = 3, at = ceiling(nH0), col = "#00000033", cex.axis = 0.8)
                graphics::abline(v = ceiling(nH0), col = "#00000033")
            }
        }
    }

    ## invisibly return power curves
    if (nullplot) {
        ret <- list("powDFH1" = powDF, "powDFH0" = powDFH0)
    } else {
        ret <- list("powDFH1" = powDF)
    }
    invisible(ret)
}
