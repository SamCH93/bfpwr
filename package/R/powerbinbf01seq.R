#' @title Design characteristics of two-stage binomial Bayes factor designs
#'
#' @description Compute design characteristics of a two-stage binomial Bayes
#'     factor design (\link{binbf01}) that stops for futility at the interim
#'     stage if the Bayes factor BF01 is greater than \code{kf}, and declares
#'     efficacy at the final stage if BF01 is equal to or below \code{k}.
#'
#' @details For some users, the low-level function \link{pbinbf01seq} that is
#'     used internally may also be useful.
#'
#' @inheritParams pbinbf01seq
#'
#' @return An object of class \code{"power.bftestseq"}, i.e., a list of the
#'     arguments, the probabilities of stopping for futility at the interim
#'     stage (\code{pintfut}, \code{pintfut0}) and of efficacy at the final
#'     stage (\code{pfineff}, \code{pfineff0}), and the expected sample sizes
#'     (\code{nexp}, \code{nexp0}), each under the design prior and under
#'     p = p0
#'
#' @author Samuel Pawel
#'
#' @seealso \link{binbf01}, \link{pbinbf01seq}
#'
#' @examples
#' ## point prior under the alternative
#' powerbinbf01seq(n1 = 20, n2 = 40, k = 1/10, kf = 3, p0 = 0.2, dp = 0.4)
#'
#' ## uniform prior [0.2, 1] under the alternative
#' powerbinbf01seq(n1 = 20, n2 = 40, k = 1/10, kf = 3, p0 = 0.2, dl = 0.2, du = 1)
#'
#' @export
powerbinbf01seq <- function(n1, n2, k = 1/10, kf = 1/k, p0 = 0.5,
                            type = c("direction"), a = 1, b = 1, dp = NA,
                            da = a, db = b, dl = 0, du = 1) {
    type <- match.arg(arg = type)

    ## probability to stop at final for efficacy (this checks all arguments)
    pfineff <- pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0,
                           type = type, a = a, b = b, dp = dp, da = da, db = db,
                           dl = dl, du = du)

    ## n1 and n2 have to be integers
    n1 <- ceiling(n1)
    n2 <- ceiling(n2)

    pfineff0 <- pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0,
                            type = type, a = a, b = b, dp = p0)

    ## probability to stop at interim for futility
    pintfut <- pbinbf01(n = n1, k = kf, p0 = p0, type = type, a = a, b = b,
                        dp = dp, da = da, db = db, dl = dl, du = du,
                        lower.tail = FALSE)
    pintfut0 <- pbinbf01(n = n1, k = kf, p0 = p0, type = type, a = a, b = b,
                         dp = p0, lower.tail = FALSE)

    ## expected sample size
    nexp <- pintfut*n1 + (1 - pintfut)*n2

    ## expected sample size under H0: p = p0
    nexp0 <- pintfut0*n1 + (1 - pintfut0)*n2

    ## return object
    structure(list(n1 = n1, n2 = n2, pintfut = pintfut, pintfut0 = pintfut0,
                   pfineff = pfineff, pfineff0 = pfineff0, nexp = nexp,
                   nexp0 = nexp0, p0 = p0, type = type, a = a, b = b, dp = dp,
                   da = da, db = db, dl = dl, du = du, k = k, kf = kf,
                   test = "binomial"),
              class = "power.bftestseq")
}

#' Print method for class \code{"power.bftestseq"}
#' @method print power.bftestseq
#'
#' @param x Object of class \code{"power.bftestseq"}
#' @param digits Number of digits for formatting of numbers
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"power.bftestseq"} object
#'
#' @author Samuel Pawel
#'
#' @export
print.power.bftestseq <- function(x, digits = getOption("digits"), ...) {
    method <- paste("One-sample directional-null binomial two-stage Bayes",
                    "factor power calculation")
    note <- paste0(
        "BF oriented in favor of H0 (BF01 < 1 indicates evidence for H1 over H0)",
        "\n      n1 and n2 are number of *observations*"
    )
    cat("\n    ", method, "\n\n")

    ## formatting BF thresholds kf and k (the object itself is not modified)
    formatThreshold <- function(threshold) {
        if (threshold < 1) {
            paste0("1/", format(x = 1/threshold, digits = digits))
        } else {
            threshold
        }
    }
    kf <- formatThreshold(x$kf)
    k <- formatThreshold(x$k)

    if (is.na(x$dp)) {
        dpvars <- c("da", "db", "dl", "du")
        dpnames <- c("design prior successes", "design prior failures",
                     "design prior lower limit", "design prior upper limit")
        condition <- "design prior"
    } else {
        dpvars <- "dp"
        dpnames <- "assumed probability"
        condition <- paste0("p = ", x$dp)
    }
    printx <- x[c("n1", "n2", "p0", "a", "b", dpvars, "kf", "k", "pintfut",
                  "pintfut0", "pfineff", "pfineff0", "nexp", "nexp0")]
    printx$kf <- kf
    printx$k <- k
    conditionH0 <- paste0("p = ", x$p0)
    names(printx) <- c("n1", "n2",
                       "null value",
                       "analysis prior successes",
                       "analysis prior failures", dpnames,
                       "BF futility threshold kf",
                       "BF efficacy threshold k",
                       paste0("Pr(BF01 > ", kf, " at interim | ", condition, ")"),
                       paste0("Pr(BF01 > ", kf, " at interim | ", conditionH0, ")"),
                       paste0("Pr(BF01 <= ", k, " at final | ", condition, ")"),
                       paste0("Pr(BF01 <= ", k, " at final | ", conditionH0, ")"),
                       paste0("E[n | ", condition, "]"),
                       paste0("E[n | ", conditionH0, "]"))

    cat(paste(paste("     ", format(names(printx), width = 10L,
                                    justify = "right")),
              format(printx, digits = digits), sep = " = "), sep = "\n")
    cat("\n", "NOTE: ", note, "\n\n", sep = "")
    invisible(x)
}
