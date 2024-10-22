nbf01. <- function(k, power, usd, null = 0, pm, psd, dpm = pm, dpsd = psd,
                   nrange = c(1, 10^5), lower.tail = TRUE, integer = TRUE,
                   analytical = TRUE, ...) {
    ## input checks
    stopifnot(
        length(k) == 1,
        is.numeric(k),
        is.finite(k),
        0 < k,

        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

        length(usd) == 1,
        is.numeric(usd),
        is.finite(usd),
        0 < usd,

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
        nrange[2] > nrange[1],

        length(lower.tail) == 1,
        is.logical(lower.tail),
        !is.na(lower.tail),

        length(integer) == 1,
        is.logical(integer),
        !is.na(integer),

        length(analytical) == 1,
        is.logical(analytical),
        !is.na(analytical)
    )

    ## use analytical solution if specified and available
    if (analytical == TRUE) {
        available <- TRUE # is analytical solution available?
        if (psd == 0)  {
            ## check whether power < limiting power
            zlim <- (null + pm - 2*dpm)*0.5/dpsd
            if (!is.nan(zlim) && pm > null) {
                powlim <- 1 - stats::pnorm(q = zlim)
            } else if (!is.nan(zlim) && pm < null) {
                powlim <- stats::pnorm(q = zlim)
            } else {
                powlim <- 0.5
            }
            if (lower.tail == FALSE) powlim <- 1 - powlim
            if (power > powlim) {
                warnmessage <- paste0("specified power (", round(power, 2),
                                      ") higher than limiting power (", round(powlim, 2),
                                      ") for specified parameters")
                warning(warnmessage)
                n <- NaN
            } else {
                zb <- stats::qnorm(p = power)
                a <- ((null + pm)*0.5 - dpm)^2 - zb^2*dpsd^2
                b <- usd^2*((null + pm - 2*dpm)*log(k)/(null - pm) - zb^2)
                c <- (usd^2*log(k)/(null - pm))^2
                if (power >= 0.5) n <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
                else n <- (-b - sqrt(b^2 - 4*a*c))/(2*a)
            }

        } else {
            available <- FALSE # analytical solution not available (currently...!)
        }
        ## else if (null == pm && pm == dpm && psd == dpsd) {
        ##     ## this is only a (quite accurate) approximation, so let's rather
        ##     ## use the more exact numerical procedure (also safes the lamW dependency!)
        ##     n <- usd^2/psd^2*k^2*exp(-lamW::lambertWm1(-k^2*stats::qnorm(p = 0.5*power)^2))
        ## }

    }

    ## otherwise use numerical solution
    if (analytical == FALSE || available == FALSE) {
        ## define function for numerical root-finding
        rootFun <- function(n) {
            pbf01(k = k, n = n, usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
                  dpsd = dpsd, lower.tail = lower.tail) - power
        }

        n <- searchN(rootFun = rootFun, nrange = nrange, ... = ...)
    }
    if (integer) return(ceiling(n))
    else return(n)
}


#' @title Sample size determination for z-test Bayes factor
#'
#' @description This function computes the required sample size to obtain a
#'     Bayes factor (\link{bf01}) more extreme than a threshold \code{k} with a
#'     specified target power.
#'
#' @inheritParams pbf01
#' @param power Target power
#' @param nrange Sample size search range over which numerical search is
#'     performed. Defaults to \code{c(1, 10^5)}
#' @param integer Logical indicating whether only integer valued sample sizes
#'     should be returned. If \code{TRUE} the required sample size is rounded to
#'     the next larger integer. Defaults to \code{TRUE}
#' @param analytical Logical indicating whether analytical (if available) or
#'     numerical method should be used. Defaults to \code{TRUE}
#' @param ... Other arguments passed to \code{stats::uniroot}
#'
#' @inherit pbf01 details
#'
#' @note A warning message will be displayed in case that the specified target
#'     power is not achievable under the specified analysis and design priors.
#'
#' @return The required sample size to achieve the specified power
#'
#' @author Samuel Pawel
#'
#' @seealso \link{pbf01}, \link{powerbf01}, \link{bf01}
#'
#' @examples
#' ## point alternative (analytical and numerical solution available)
#' nbf01(k = 1/10, power = 0.9, usd = 1, null = 0, pm = 0.5, psd = 0,
#'       analytical = c(TRUE, FALSE), integer = FALSE)
#'
#' ## standardized mean difference (usd = sqrt(2), effective sample size = per group size)
#' nbf01(k = 1/10, power = 0.9, usd = sqrt(2), null = 0, pm = 0, psd = 1)
#' ## this is the sample size per group (assuming equally sized groups)
#'
#' ## z-transformed correlation (usd = 1, effective sample size = n - 3)
#' nbf01(k = 1/10, power = 0.9, usd = 1, null = 0, pm = 0.2, psd = 0.5)
#' ## have to add 3 to obtain the actual sample size
#'
#' ## log hazard/odds ratio (usd = 2, effective sample size = total number of events)
#' nbf01(k = 1/10, power = 0.9, usd = 2, null = 0, pm = 0, psd = sqrt(0.5))
#' ## have to convert the number of events to a sample size
#'
#' @export
nbf01 <- Vectorize(FUN = nbf01.,
                   vectorize.args = c("k", "power", "usd", "null", "pm", "psd",
                                      "dpm", "dpsd", "integer", "analytical"))



## ## should give the same as closed-form solutions
## ## - for point alternatives and design priors (should be exact)
## usd <- 2
## pm <- 1
## psd <- 0
## dpm <- 2
## dpsd <- 0
## null <- 0.5
## beta <- 0.2
## k <- 1/10
## zb <- stats::qnorm(p = 1 - beta)
## (zb + sqrt(zb^2 - log(k^2)*(null + pm - 2*dpm)/(null - pm)))^2/(null + pm - 2*dpm)^2*usd^2
## nbf01(k = k, power = 1 - beta, usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
##       dpsd = dpsd, integer = FALSE)
## ## - for local alternatives and design priors (this is an approximation)
## psd <- dpsd <- 1
## pm <- dpm <- null
## usd^2/psd^2*k^2*exp(-lamW::lambertWm1(x = -stats::qnorm(p = (1 + beta)/2)^2*k^2))
## nbf01(k = k, power = 1 - beta, usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
##       dpsd = dpsd, integer = FALSE)

## ## - for point alternatives and normal design priors (should be exact)
## usd <- 2
## pm <- 1
## psd <- 0
## dpm <- 0.8
## dpsd <- 0
## null <- -0.5
## beta <- 0.2
## k <- 1/10
## zb <- stats::qnorm(p = 1 - beta)
## (n1 <- nbf01(k = k, power = 1 - beta, usd = usd, null = null, pm = pm, psd = psd,
##              dpm = dpm, dpsd = dpsd, integer = FALSE))
## pbf01(k = k, n = n1, usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
##       dpsd = dpsd)

## A <- zb^2*dpsd^2 - ((null + pm)/2 - dpm)^2
## B <- zb^2*usd^2 - (null + pm - 2*dpm)*usd^2*log(k)/(null - pm)
## C <- -(usd^2*log(k)/(null - pm))^2
## (n2 <- (-B + c(-1, 1)*sqrt(B^2 - 4*A*C))/(2*A))
## pbf01(k = k, n = n2[1], usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
##       dpsd = dpsd)

## ## - for local alternatives and design priors (this is an approximation)
## psd <- dpsd <- 1
## pm <- dpm <- null
## (n2 <- usd^2/psd^2*k^2*exp(-lamW::lambertWm1(x = -stats::qnorm(p = (1 + beta)/2)^2*k^2)))
## n2b <- nbf01(k = k, power = 1 - beta, usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
##       dpsd = dpsd, integer = FALSE, analytical = FALSE)
## pbf01(k = k, n = n2, usd = usd, null = null, pm = pm, psd = psd, dpm = dpm,
##       dpsd = dpsd)

## ## produce some tables
## ## - point alternatives with SMD = 1
## kseq <- rev(c(1/1000, 1/300, 1/100, 1/30, 1/10, 1/3, 1/2))
## powseq <- seq(0.5, 0.95, 0.05)
## tab1 <- sapply(X = kseq, FUN = function(k) {
##     beta <- 1 - powseq
##     zb <- stats::qnorm(p = 1 - beta)
##     (zb + sqrt(zb^2 - log(k^2)))^2
## })
## colnames(tab1) <- BayesRep::formatBF(kseq)
## rownames(tab1) <- powseq*100
## xtab1 <- xtable::xtable(ceiling(tab1), digits = 0)
## xtable::print.xtable(xtab1, booktabs = TRUE)
## ## - local normal unit-information alternative
## tab2 <- sapply(X = kseq, FUN = function(k) {
##     beta <- 1 - powseq
##     k^2*exp(-lamW::lambertWm1(x = -k^2*qnorm(p = (1 + beta)/2)^2))
## })
## colnames(tab2) <- BayesRep::formatBF(kseq)
## rownames(tab2) <- powseq*100
## xtab2 <- xtable::xtable(ceiling(tab2), digits = 0)
## xtable::print.xtable(xtab2, booktabs = TRUE)
