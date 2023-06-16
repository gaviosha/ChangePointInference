#'@title Change point inference via differencing
#'@description
#' TO-DO!
#'
#'@param xx signal to be segmented
#'@param degree polynomial degree of the underlying signal
#'@param alpha desired coverage 
#'@param gaussian_noise whether the noise is Gaussian
#'@param independent_noise whether the noise is serially dependent
#'@param tau noise level if know, if null the noise level will be estimated via MAD if the noise in independent and via TACV estimator if the noise is dependent
#'@param aa controls grid density (a in the paper)
#'@param min_scale minimum width of intervals tested for a change
#'@param HH supply numeric constant H if pre-compute
#'@export  

diffInf <- function(xx, degree, alpha = 0.1, gaussian_noise = TRUE, independent_noise = TRUE, tau = NULL, aa = sqrt(2), min_scale = floor(sqrt(length(xx))), HH = NULL)
{
  nn <- length(xx)
  
  xx_cumsum <- c(0,cumsum(xx))
  
  scaling <- sum(choose(degree+1,0:(degree+1))**2) / (degree+2)
  
  
  noise_level_know <- !is.null(tau)
  
  if (!noise_level_know)
  {
    if (gaussian_noise && independent_noise) tau <- gen_diff_mad(xx, degree)
    
    if (!gaussian_noise && independent_noise) tau <- gen_diff_sd(xx, degree)
    
    if (!independent_noise) tau <- gen_diff_catoni_tavc(xx, (degree+2)*min_scale, degree)
  }
  
  
  if (gaussian_noise && independent_noise) min_scale <- log(nn)
  
  thresh <- tau * get_thresh(nn, min_scale, alpha, degree, aa, HH, (gaussian_noise && independent_noise))
  
  ints_df <- data.frame(matrix(nrow = 0, ncol = 3))
  
  names(ints_df) <- c("start", "end", "value")
  
  ints_out <- diff_bin_seg(xx_cumsum, 1, nn, degree, aa, min_scale, thresh, scaling)

  ints_out <- rbind(ints_df,ints_out)
  
  
  out <- list(intervals = ints_out, thresh = thresh, data = xx, degree = degree)
  
  class(out) <- c("cptInference", class(out))
  
  return(out)
}
