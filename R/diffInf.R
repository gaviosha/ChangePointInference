
diffInf <- function(xx, degree, alpha, noise_type = c("gaussian","non_gaussian_dependent")[1], dependent_noise = FALSE, tau = NULL, aa = sqrt(2), min_scale = log(length(xx)), HH = NULL, WW = NULL, tacv_max_scale = NULL)
{
  #' Change point inference via differencing
  #'
  #' Main routine
  #'
  #'@param xx signal to be segmented
  #'@param degree polynomial degree of the underlying signal
  #'@param alpha desired coverage 
  #'@param noise_type one of "gaussian" or "non_gaussian_dependent"
  #'@param dependent_noise whether the noise is serially dependent
  #'@param tau noise level if know, if null the noise level will be estimated via MAD if the noise in independent and via TACV estimator if the noise is dependent
  #'@param aa controls grid density (a in the paper)
  #'@param min_scale minimum width of intervals tested for a change
  #'@param HH supply numeric constant H if pre-computed
  #'@param WW minimum segment length in the case of non-gaussian or dependent noises (W in the paper)
  #'@param tacv_max_scale coarsest scale on which TAVC will be calculated, all coarser scales will use TAVC calculated on this scale
  #'
  #'@export  
  
  nn <- length(xx)
  
  xx_cumsum <- c(0,cumsum(xx))
  
  scaling <- sum(choose(degree+1,0:(degree+1))**2) / (degree+2)
  
  thresh <- get_thresh(nn, min_scale, alpha, degree, aa, HH, noise_type)
  
  
  noise_level_know <- !is.null(tau)
  
  if (dependent_noise && is.null(tacv_max_scale)) tacv_max_scale <- floor(2.5*sqrt(nn))
  
  if (!noise_level_know) tau <- ifelse(dependent_noise, generalised_tavc_est(xx_cumsum, tacv_max_scale, degree, scaling, tacv_max_scale), stats::mad(diff(xx/sqrt(2))))
  
  
  ints_df <- data.frame(matrix(nrow = 0, ncol = 3))
  
  names(ints_df) <- c("start", "end", "value")
  
  ints_out <- diff_bin_seg(xx_cumsum, 1, nn, degree, aa, min_scale, thresh, scaling, tau, tacv_max_scale, noise_level_know, noise_type, dependent_noise)
  
  ints_out <- rbind(ints_df,ints_out)
  
  return(
    list(
    intervals = ints_out, 
    thresh = thresh, 
    data = xx
    ))
}
