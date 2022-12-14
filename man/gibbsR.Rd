% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp22012R.R
\name{gibbsR}
\alias{gibbsR}
\title{A Gibbs sampler to generate a bivariate normal chain}
\usage{
gibbsR(N, rho, mu1, mu2, sigma1, sigma2)
}
\arguments{
\item{N}{the number of samples}

\item{rho}{the correlation of two dimensions}

\item{mu1}{the mean of the first dimension}

\item{mu2}{the mean of the second dimension}

\item{sigma1}{the standard deviation of the first dimension}

\item{sigma2}{the standard deviation of the second dimension}
}
\value{
a random sample of size \code{N}
}
\description{
A Gibbs sampler to generate a bivariate normal chain
}
\examples{
\dontrun{
rnR <- gibbsR(N=2000, mu1=0, mu2=0, sigma1=1, sigma2=1, rho=0.9)
par(mfrow=c(1,3));
plot(rnR[,1],type='l')
plot(rnR[,2],type='l')
plot(rnR)
}
}
