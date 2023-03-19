loc_diff <- function(xx_cumsum, ll, ww, degree)
{
  #' Local differences
  #'
  #'@param xx_cumsum cumulative sum of data
  #'@param ll start of the segment
  #'@param ww number of data points in the segment
  #'@param degree polynomial degree of the underlying signal
  #'
  #'@export
  
  local_sums <- sapply(0:(degree+1), function(jj) xx_cumsum[ll+((jj+1)*ww/(degree+2))] - xx_cumsum[ll+(jj*ww/(degree+2))])
  
  local_diff <- diff(local_sums, differences = degree + 1)
  
  return(abs(local_diff))
}


diff_bin_seg <- function(xx_cumsum, ss, ee, degree, aa, min_width, thresh, scaling)
{
  #' diffInf binary segmentation
  #'
  #' Recursively search for narrowest interval containing an change point over a particular logarithmic grid using binary segmentation algorithm; not not be called directly by the user. 
  #'
  #'@param xx_cumsum cumulative sum of data
  #'@param ss start of search area 
  #'@param ee end of search area
  #'@param degree polynomial degree of the underlying signal
  #'@param aa controls grid density (a in the paper)
  #'@param min_width minimum width of intervals tested for a change
  #'@param thresh threshold to control FWE computed by `diffInf`
  #'@param scaling numeric calling constant in local tests pre-computed by `diffInf`
  #'
  #'@export  
  
  nn_loc <- (ee - ss + 1)
  
  if (nn_loc <= min(min_width, degree + 2)) return(NULL)
  
  end_search <- FALSE
  
  test_passed <- FALSE
  
  for (jj in 1:floor(logb(nn_loc,aa)))
  {
    ww <- (degree+2) * floor((aa**jj)/(degree+2))
    
    if ((ww > min_width) & (ww < nn_loc)) 
    {
      for (ll in ss:(ee-ww))
      {
        
        local_test <- loc_diff(xx_cumsum, ll, ww, degree) / sqrt(ww*scaling)
        
        if (local_test > thresh)
        {
          end_search <- TRUE
          
          test_passed <- TRUE
          
          break
        }
      }
    }
    
    if (end_search) break
  }
  
  if (!test_passed) return(NULL)
  
  left <- diff_bin_seg(xx_cumsum, ss, max(ss,ll-1), degree, aa, min_width, thresh, scaling)
  
  right <- diff_bin_seg(xx_cumsum, min(ee,ll+ww), ee, degree, aa, min_width, thresh, scaling)
  
  return(rbind(left, c(ll, ll+ww-1, local_test), right))
}