nbinbf01. <- function(k, power, p0 = 0.5, type = c("point", "direction"), a = 1,
                      b = 1, dp = NA, da = a, db = b, dl = 0, du = 1,
                      lower.tail = TRUE, nrange = c(1, 10^4), ...) {
    ## input checks
    stopifnot(
        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

        length(nrange) == 2,
        all(is.numeric(nrange)),
        all(is.finite(nrange)),
        nrange[2] > nrange[1],
        nrange[1] > 0
    )


    ## define function for numerical root-finding
    rootFun <- function(n) {
        suppressWarnings({
            pbinbf01(k = k, n = n, p0 = p0, type = type, a = a, b = b, dp = dp,
                     da = da, db = db, dl = dl, du = du,
                     lower.tail = lower.tail) - power
        })
    }

    n <- searchNoscil(rootFun = rootFun, nrange = nrange, nextend = 10,
                      maxcycles = 5, ... = ...)
    return(ceiling(n))
}


#' @title Sample size determination for binomial Bayes factor
#'
#' @description This function computes the required sample size to obtain a
#'     binomial Bayes factor (\link{binbf01}) more extreme than a threshold
#'     \code{k} with a specified target power.
#'
#' @inheritParams pbinbf01
#' @param power Target power
#' @param nrange Sample size search range over which numerical search is
#'     performed. Defaults to \code{c(1, 10^4)}
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbinbf01}, \link{binbf01}
#'
#' @examples
#' ## sample size parameters
#' pow <- 0.9
#' p0 <- 3/4
#' a <- 1
#' b <- 1
#' k <- 1/10
#'
#' \dontrun{
#' ## sample sizes for directional testing
#' (nH1 <- nbinbf01(k = k, power = pow, p0 = p0, type = "direction", a = a,
#'                  b = b, da = a, db = b, dl = p0, du = 1))
#' (nH0 <- nbinbf01(k = 1/k, power = pow, p0 = p0, type = "direction", a = a,
#'                  b = b, da = a, db = b, dl = 0, du = p0, lower.tail = FALSE))
#' nseq <- seq(1, 1.1*max(c(nH1, nH0)), length.out = 100)
#' powH1 <- pbinbf01(k = k, n = nseq, p0 = p0, type = "direction", a = a,
#'                   b = b, da = a, db = b, dl = p0, du = 1)
#' powH0 <- pbinbf01(k = 1/k, n = nseq, p0 = p0, type = "direction", a = a,
#'                   b = b, da = a, db = b, dl = 0, du = p0, lower.tail = FALSE)
#' matplot(nseq, cbind(powH1, powH0), type = "s", xlab = "n", ylab = "Power", lty = 1,
#'         ylim = c(0, 1), col = c(2, 4), las = 1)
#' abline(h = pow, lty = 2)
#' abline(v = c(nH1, nH0), col = c(2, 4), lty = 2)
#' legend("topleft", legend = c("H1", "H0"), lty = 1, col = c(2, 4))
#'
#' ## sample sizes for point null testing
#' (nH1 <- nbinbf01(k = k, power = pow, p0 = p0, type = "point", a = a,
#'                  b = b, da = a, db = b))
#' (nH0 <- nbinbf01(k = 1/k, power = pow, p0 = p0, type = "point", a = a,
#'                  b = b, dp = p0, lower.tail = FALSE, nrange = c(1, 10^5)))
#' nseq <- seq(1, max(c(nH1, nH0)), length.out = 100)
#' powH1 <- pbinbf01(k = k, n = nseq, p0 = p0, type = "point", a = a,
#'                   b = b, da = a, db = b, dl = 0, du = 1)
#' powH0 <- pbinbf01(k = 1/k, n = nseq, p0 = p0, type = "point", a = a,
#'                   b = b, dp = p0, lower.tail = FALSE)
#' matplot(nseq, cbind(powH1, powH0), type = "s", xlab = "n", ylab = "Power", lty = 1,
#'         ylim = c(0, 1), col = c(2, 4), las = 1)
#' abline(h = pow, lty = 2)
#' abline(v = c(nH1, nH0), col = c(2, 4), lty = 2)
#' legend("topleft", legend = c("H1", "H0"), lty = 1, col = c(2, 4))
#' }
#' @export
nbinbf01 <- Vectorize(FUN = nbinbf01.,
                      vectorize.args = c("k", "power", "p0", "type", "a", "b",
                                         "dp", "da", "db", "dl", "du",
                                         "lower.tail"))
