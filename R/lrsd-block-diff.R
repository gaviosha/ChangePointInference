#' (Long Run) Standard Deviation Estimator Based on \eqn{(p+1)}-th Differences of Local Sums of the Data
#'@description Estimates the (long run) standard deviation based on \eqn{(p+1)}-th differences of non-overlapping local sums of the data. 
#'@param yy A numeric vector containing the data. 
#'@param ww The scale at which local sums of the data will be calculated. 
#'@param degree The degree of polynomial parametrization of the mean of the data. 
#'@examples 
#' degree <- 3 
#' 
#' nn <- 500 
#' 
#' yy <- arima.sim(model = list(ar = 0.5), n = nn) + (1:nn)^{degree}
#' 
#' sd_diff(yy, degree)
#'@export
lrsd_block_diff <- function(yy, ww = length(yy)^{(1/3)}, degree = 0)
{
  nn <- length(yy)
  
  mm <- floor(nn/ww)
  
  Cp <- sum(choose(degree + 1, 0:(degree + 1))^2)
  
  blocks <- split(yy, cut(seq_along(yy), mm, labels = FALSE)) 
  
  block_partial_sums <- sapply(blocks, function(ii) sum(ii) / sqrt(ww))
  
  block_diffs <- diff(block_partial_sums, differences = degree + 1)
  
  sqrt(mean(block_diffs^2 / Cp))
}
