#' @title Power calculations for binomial Bayes factor with two-stage designs
#'
#' @description Compute design characteristics for two-stage binomial Bayes
#'     factor design which stops for futility if the binomial Bayes factor
#'     (\link{binbf01}) is equal or below a threshold at final stage, taking
#'     into account that the study can stop if it is greater than a specified
#'     threshold at interim.
#'
#' @details For some users, the low-level function \link{pbinbf01seq} that is
#'     used internally may also be useful.
#'
#' @inheritParams pbinbf01seq
#'
#' @return An object of type \code{"power.bftestseq"}, i.e., a list of the
#'     arguments, computed probabilitites and expected samples sizes, and
#'     \code{method} and \code{note} elements
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


    ## probability to stop at interim for futility
    pintfut <- pbinbf01(n = n1, k = kf, p0 = p0, type = type, a = a, b = b,
                        dp = dp, da = da, db = db, dl = dl, du = du,
                        lower.tail = FALSE)
    pintfut0 <- pbinbf01(n = n1, k = kf, p0 = p0, type = type, a = a, b = b,
                         dp = p0, lower.tail = FALSE)

    ## probability to stop at final for efficacy
    pfineff <- pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0,
                           type = type, a = a, b = b, dp = dp, da = da, db = db,
                           dl = dl, du = du)
    pfineff0 <- pbinbf01seq(n1 = n1, n2 = n2, k = k, kf = kf, p0 = p0,
                            type = type, a = a, b = b, dp = p0)

    ## expected sample size
    nexp <- pintfut*n1 + (1 - pintfut)*n2

    ## expected sample size under H0: p = p0
    nexp0 <- pintfut0*n1 + (1 - pintfut0)*n2

    ## return object
    structure(list(n1 = n1, n2 = n2, pintfut = pintfut, pintfut0 = pintfut0,
                   pfineff = pfineff, pfineff0 = pfineff0, nexp = nexp,
                   nexp0 = nexp0, p0 = p0, type = type, a = a, b = b, dp = dp,
                   da = da, db = db, dl = dl, du = du, k = k, kf = kf,
                   type = type, test = "binomial"),
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

    if (x$test %in% c("binomial")) {
        methodstring <- "binomial"
    } else {
        methodstring <- NA
    }
    method <- paste(methodstring, c("two-stage Bayes factor power calculation"))
    note <- "BF oriented in favor of H0 (BF01 < 1 indicates evidence for H1 over H0)"
    if (x$test == "binomial") {
        note <- paste(note,
                      "n1 and n2 are number of *observations*",
                      "Plot power to check that it does not decrease below target for larger n!",
                      sep = "\n      ")
        if (x$type == "point") {
            testtype <- "point-null"
        } else {
            testtype <- "directional-null"
        }
        method <- paste("One-sample", testtype, method)
    } else {
        note <- NA
        method <- NA
    }
    cat("\n    ", method, "\n\n")

     ## formatting BF threshold k
    if (x$k < 1) {
        x$k <- paste0("1/", format(x = 1/x$k, digits = digits))
    }
    if (x$kf < 1) {
        x$kf <- paste0("1/", format(x = 1/x$kf, digits = digits))
    }
    if (x$test == "binomial") {
        if (is.na(x$dp)) {
            dpvars <- c("da", "db", "dl", "du")
            dpnames <- c("design prior successes", "design prior failures",
                         "design prior lower limit", "design prior upper limit")
            condition <- "design prior"
        } else {
            dpvars <- "dp"
            dpnames <- "assumed propability"
            condition <- paste0("p = ", x$dp)
        }
        printx <- x[c("n1", "n2", "p0", "a", "b", dpvars, "kf", "k", "pintfut",
                      "pintfut0", "pfineff", "pfineff0", "nexp", "nexp0")]
        names(printx) <- c("n1", "n2",
                           "null value",
                           "analysis prior successes",
                           "analysis prior failures", dpnames,
                           "BF futility threshold k",
                           "BF efficacy threshold k",
                           ## paste0("Pr(futility at interim |", condition, ")"),
                           ## paste0("Pr(futility at interim | p = ", x$p0, ")"),
                           ## paste0("Pr(efficacy at final |", condition, ")"),
                           ## paste0("Pr(efficacy at final | p = ", x$p0, ")"),
                           paste0("Pr(BF01 > ", x$kf, " at interim | ", condition, ")"),
                           paste0("Pr(BF01 > ", x$kf, " at interim | p = ", x$p0, ")"),
                           paste0("Pr(BF01 < ", x$k, " at final | ", condition, ")"),
                           paste0("Pr(BF01 < ", x$k, " at final | p = ", x$p0, ")"),
                           paste0("E[n | ", condition, "]"),
                           paste0("E[n | p = ", x$p0, "]"))
    }

    cat(paste(paste("     ", format(names(printx), width = 10L, justify = "right")),
              format(printx, digits = digits), sep = " = "), sep = "\n")
    cat("\n", "NOTE: ", note, "\n\n", sep = "")
    invisible(x)
}
