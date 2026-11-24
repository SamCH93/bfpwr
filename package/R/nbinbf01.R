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


    type <- match.arg(type)

    ## define function for numerical root-finding
    rootFun <- function(n) {
        suppressWarnings({
            pbinbf01(k = k, n = n, p0 = p0, type = type, a = a, b = b, dp = dp,
                     da = da, db = db, dl = dl, du = du,
                     lower.tail = lower.tail) - power
        })
    }

    ## Binomial power can oscillate in n (saw-tooth) because the data are
    ## discrete, a known feature of discrete-data models that is qualitatively
    ## different from a non-monotone power curve. As is standard for
    ## discrete-data sample sizes, searchNoscil() deliberately checks that
    ## power does not drop below the target over the next 10 sample sizes
    ## after a crossing instead of returning the first n that reaches the
    ## target (see Details in ?nbinbf01 for the one-step advance and the
    ## restart limit). Keep this rule; do not replace it with the plain
    ## searchN() crossing used by the continuous-data functions.
    n <- searchNoscil(rootFun = rootFun, nrange = nrange, nextend = 10,
                      maxcycles = 5, ...)
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
#' @param ... Other arguments passed to \code{stats::uniroot}. The sample-size
#'     root tolerance defaults to \code{tol = 1e-8} and can be overridden.
#'
#' @details Because the number of successes is discrete, binomial power is not
#'     a smooth function of the sample size but can oscillate in a saw-tooth
#'     (zig-zag) pattern: power can reach the target at some \code{n} and fall
#'     below it again at a slightly larger \code{n} (Chernick and Liu, 2002).
#'     Such oscillation is a known feature of the binomial and other
#'     discrete-data models, where it also affects, for example, the coverage
#'     of confidence intervals and the type I error rate of exact tests. It is
#'     qualitatively different from a non-monotone power curve for continuous
#'     data, such as the probability of misleading evidence (see
#'     \link{nbf01}). As is common in sample size calculations for discrete
#'     data (e.g., Julious, 2023, Chapter 17), the function therefore aims
#'     for a sample size after which power does not drop below the target
#'     again, by checking a window of sample sizes after each crossing, rather
#'     than simply returning the first sample size at which power reaches the
#'     target. The two can coincide when power does not dip.
#'
#'     Specifically, the search locates a crossing of the target power with
#'     \code{stats::uniroot} and checks power at the crossing and at the next
#'     10 sample sizes. If all of them meet the target, the crossing is
#'     returned. If only the crossing itself falls below the target, the next
#'     sample size is returned (so 9 further sample sizes have been checked).
#'     Otherwise, the search is restarted from the last checked sample size
#'     below the target, with at most 5 restarts; if these are used up, a
#'     warning is issued and the last crossing is returned, at which power may
#'     be below the target. The returned sample size can thus be larger than
#'     the first sample size at which power reaches the target, and power can
#'     still fall below the target beyond the checked sample sizes;
#'     \link{pbinbf01} shows power over a range of sample sizes. \link{nbf01},
#'     \link{nnmbf01}, and \link{ntbf01} determine the sample size without
#'     this check.
#'
#' @return The sample size selected by the search, rounded up to an integer,
#'     or \code{NaN} with a warning if the numerical search fails
#'
#' @references Chernick, M. R., Liu, C. Y. (2002). The saw-toothed behavior of
#'     power versus sample size and software solutions: Single binomial
#'     proportion using exact methods. The American Statistician,
#'     56(2):149-155. \doi{10.1198/000313002317572835}
#'
#' Julious, S. A. (2023). Sample Sizes for Clinical Trials (2nd ed.). Chapman
#'     and Hall/CRC. \doi{10.1201/9780429503658}
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbinbf01}, \link{binbf01}
#'
#' @examples
#' ## binomial power oscillates with the sample size: power reaches 0.8 at
#' ## n = 74, 76, and 77, but nbinbf01() returns n = 79, after which power
#' ## stays at or above 0.8 for the next 10 sample sizes
#' (n <- nbinbf01(k = 1/6, power = 0.8, p0 = 0.5, type = "point", dp = 0.7))
#' nseq <- 72:(n + 10)
#' round(setNames(pbinbf01(k = 1/6, n = nseq, p0 = 0.5, type = "point",
#'                         dp = 0.7), nseq), 3)
#'
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
