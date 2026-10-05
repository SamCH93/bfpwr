## Reproducible example for the PR review: log.p alone does not prevent
## cancellation when both log probabilities are approximately -5e11.
library(bfpwr)

logbf_simple <- function(z) {
    m <- -1
    r <- 1e-6
    priorz <- m/r
    postz <- (priorz + r*z)/sqrt(1 + r^2)
    log1p(r^2)/2 + (m*(m - 2*z) - r^2*z^2)/(2*(1 + r^2)) +
        pnorm(priorz, log.p = TRUE) - pnorm(postz, log.p = TRUE)
}
critical <- uniroot(logbf_simple, c(-1, 1), tol = 1e-10)$root
c(simplified = pnorm(critical, lower.tail = FALSE),
  safeguarded = pbf01(k = 1, n = 1, usd = 1, pm = -1, psd = 1e-6,
      dpm = 0, dpsd = 0, alternative = "greater"))
## simplified safeguarded
##  0.8413447   0.5000000
