% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitPeak.R
\name{fitPeak}
\alias{fitPeak}
\title{fitPeak}
\usage{
fitPeak(x, y, startPars = c(mu = 70.065, I = 5000, sigma = 3e-04),
  maxAttempts = 1000, sigma = 0)
}
\arguments{
\item{x}{vector of x-values;}

\item{y}{vector of y-values;}

\item{startPars}{names vector of initial parameters
\code{mu}, \code{I} and \code{sigma};}

\item{maxAttempts}{maximum number of attempts, see below;}

\item{sigma}{if \code{sigma>0} provided, two-parameter fitting is performed,
fixing \code{sigma} at provided value;}
}
\value{
Returns an object of class \code{nls} if the fitting was successful.
}
\description{
Fits gaussian peak for the data
\code{ y = I/sigma * exp(-1/2*(x-mu)^2 / sigma^2) }
using nls function with initial parameters specified by named vector
\code{startPars}. If \code{sigma > 0} orivuded, then only two parameters are
fitted using nls - intensity \code{I} and position \code{mu}.
}
\examples{
fitPeak(x=x, y=y, startPars=c('mu' = 70.065, 'I' = 5e3, 'sigma' = 3e-4))
fitPeak(x=x, y=y, startPars=c('mu' = 70.065, 'I' = 5e3, 'sigma' = 3e-4), sigma=3e-4)
}
