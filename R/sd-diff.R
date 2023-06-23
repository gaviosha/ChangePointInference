#'Standard Deviation Estimator Based On \eqn{(p+1)}-th Differences
#'@description Estimates the standard deviation of the data using the mean on the squared \eqn{(p+1)}-th difference of the data sequence. Suitable for estimating the scale of the noise when the data sequence consists of a piecewise polynomial function contaminated with independently distributed noise having bounded fourth moment. 
#'@param yy A numeric vector containing the data. 
#'@param degree The degree of polynomial parametrization of the mean of the data. 
#'@examples
#' degree <- 3 
#' 
#' nn <- 500 
#' 
#' yy <- rt(nn, df = 5) * sqrt(3/5) + (1:nn)^{degree}
#' 
#' sd_diff(yy, degree)
#'@export
sd_diff <- function(yy, degree)
{
  nn <- length(yy)
  
  cp <- sum(choose(degree+1, 0:(degree+1))**2)
  
  yy_diff <- diff(yy, differences = (degree+1))
  
  var_est <- sum((yy_diff**2)/cp) / (nn - (degree+1))
  
  return(sqrt(var_est))
}