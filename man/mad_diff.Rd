% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mad-diff.R
\name{mad_diff}
\alias{mad_diff}
\title{Median Absolute Deviation Estimator Using \eqn{(p+1)}-th Differences}
\usage{
mad_diff(yy, degree)
}
\arguments{
\item{yy}{A numeric vector containing the data.}

\item{degree}{The degree of polynomial parametrization of the mean of the data.}
}
\description{
Estimates the standard deviation using the median of the \eqn{(p+1)}-th difference of the data. Suitable for estimating the scale of the noise when the data sequence is Gaussian with a piecewise polynomial mean of degree \eqn{p}.
}
\examples{
degree <- 3 

nn <- 500 

yy <- rnorm(nn) + (1:nn)^{degree}

mad_diff(yy, degree)
}
