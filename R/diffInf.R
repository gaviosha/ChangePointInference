
diffInf <- function(xx, degree, alpha, sigma = NULL, aa = 2, LL = 0, HH = NULL)
{
  #' Change point inference via differencing
  #'
  #' Main routine
  #'
  #'@param xx signal to be segmented
  #'@param degree polynomial degree of the underlying signal
  #'@param alpha desired coverage 
  #'@param sigma noise level if know, if null sigma will be estimated using 
  #'@param aa controls grid density (a in the paper)
  #'@param LL control smallest scale (L in the paper)
  #'@param HH supply numeric constant H if pre-computed
  #'
  #'@export  
  
  nn <- length(xx)
  
  xx_cumsum <- c(0,cumsum(xx))
  
  sigma <- ifelse(is.null(sigma), stats::mad(diff(xx/sqrt(2))), sigma)
  
  thresh <- sigma * get_thresh(nn, alpha, degree, aa, LL, HH)
  
  min_width <- ceiling(log(nn)/aa**LL)
  
  scaling <- sum(choose(degree+1,0:(degree+1))**2) / (degree+2)
  
  ints_df <- data.frame(matrix(nrow = 0, ncol = 3))
  
  names(ints_df) <- c("start", "end", "value")
  
  ints_out <- diff_bin_seg(xx_cumsum, 1, nn, degree, aa, min_width, thresh, scaling)
  
  ints_out <- rbind(ints_df,ints_out)
  
  return(
    list(
    intervals = ints_out, 
    thresh = thresh, 
    data = xx
    ))
}
