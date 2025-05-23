% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbf01.R
\name{nbf01}
\alias{nbf01}
\title{Sample size determination for z-test Bayes factor}
\usage{
nbf01(
  k,
  power,
  usd,
  null = 0,
  pm,
  psd,
  dpm = pm,
  dpsd = psd,
  nrange = c(1, 10^5),
  lower.tail = TRUE,
  integer = TRUE,
  analytical = TRUE,
  ...
)
}
\arguments{
\item{k}{Bayes factor threshold}

\item{power}{Target power}

\item{usd}{Unit standard deviation, the (approximate) standard error of the
parameter estimate based on \eqn{\code{n}=1}{n=1}, see details}

\item{null}{Parameter value under the point null hypothesis. Defaults to
\code{0}}

\item{pm}{Mean of the normal prior assigned to the parameter under the
alternative in the analysis}

\item{psd}{Standard deviation of the normal prior assigned to the parameter
under the alternative in the analysis. Set to \code{0} to obtain a point
prior at the prior mean}

\item{dpm}{Mean of the normal design prior assigned to the parameter.
Defaults to the same value as the analysis prior \code{pm}}

\item{dpsd}{Standard deviation of the normal design prior assigned to the
parameter. Defaults to the same value as the analysis prior \code{psd}}

\item{nrange}{Sample size search range over which numerical search is
performed. Defaults to \code{c(1, 10^5)}}

\item{lower.tail}{Logical indicating whether Pr(\eqn{\mathrm{BF}_{01}}{BF01}
\eqn{\leq}{<=} \code{k}) (\code{TRUE}) or Pr(\eqn{\mathrm{BF}_{01}}{BF01}
\eqn{>} \code{k}) (\code{FALSE}) should be computed. Defaults to
\code{TRUE}}

\item{integer}{Logical indicating whether only integer valued sample sizes
should be returned. If \code{TRUE} the required sample size is rounded to
the next larger integer. Defaults to \code{TRUE}}

\item{analytical}{Logical indicating whether analytical (if available) or
numerical method should be used. Defaults to \code{TRUE}}

\item{...}{Other arguments passed to \code{stats::uniroot}}
}
\value{
The required sample size to achieve the specified power
}
\description{
This function computes the required sample size to obtain a
    Bayes factor (\link{bf01}) more extreme than a threshold \code{k} with a
    specified target power.
}
\details{
It is assumed that the standard error of the future parameter
    estimate is of the form \eqn{\code{se} =\code{usd}/\sqrt{\code{n}}}{se =
    usd/sqrt(n)}. For example, for normally distributed data with known
    standard deviation \code{sd} and two equally sized groups of size
    \code{n}, the standard error of an estimated standardized mean difference
    is \eqn{\code{se} = \code{sd}\sqrt{2/n}}{se = sd*sqrt(2/n)}, so the
    corresponding unit standard deviation is \eqn{\code{usd} =
    \code{sd}\sqrt{2}}{usd = sd*sqrt(2)}. See the vignette for more
    information.
}
\note{
A warning message will be displayed in case that the specified target
    power is not achievable under the specified analysis and design priors.
}
\examples{
## point alternative (analytical and numerical solution available)
nbf01(k = 1/10, power = 0.9, usd = 1, null = 0, pm = 0.5, psd = 0,
      analytical = c(TRUE, FALSE), integer = FALSE)

## standardized mean difference (usd = sqrt(2), effective sample size = per group size)
nbf01(k = 1/10, power = 0.9, usd = sqrt(2), null = 0, pm = 0, psd = 1)
## this is the sample size per group (assuming equally sized groups)

## z-transformed correlation (usd = 1, effective sample size = n - 3)
nbf01(k = 1/10, power = 0.9, usd = 1, null = 0, pm = 0.2, psd = 0.5)
## have to add 3 to obtain the actual sample size

## log hazard/odds ratio (usd = 2, effective sample size = total number of events)
nbf01(k = 1/10, power = 0.9, usd = 2, null = 0, pm = 0, psd = sqrt(0.5))
## have to convert the number of events to a sample size

}
\seealso{
\link{pbf01}, \link{powerbf01}, \link{bf01}
}
\author{
Samuel Pawel
}
