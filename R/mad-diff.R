#' Median Absolute Deviation Estimator Using \eqn{(p+1)}-th Differences
#'@description Estimates the standard deviation using the median of the \eqn{(p+1)}-th difference of the data. Suitable for estimating the scale of the noise when the data sequence is Gaussian with a piecewise polynomial mean of degree \eqn{p}. 
#'@param yy A numeric vector containing the data. 
#'@param degree The degree of polynomial parametrization of the mean of the data. 
#'@examples
#' degree <- 3 
#' 
#' nn <- 500 
#' 
#' yy <- rnorm(nn) + (1:nn)^{degree}
#' 
#' mad_diff(yy, degree)
#'@export 
mad_diff <- function(yy, degree)
{
  yy_diff <- diff(yy, differences = (degree+1))
  
  cp <- sum(choose(degree+1, 0:(degree+1))**2)
  
  return(mad(yy_diff/ sqrt(cp)))
}
