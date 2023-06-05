#'@title Robust estimation of TAVC with unknown piecewise polynomial signal
#'@description A short description...
#'@param xx_cumus, cumulative sum of the data 
#'@param ww, scale at which TVAC will be calculated  
#'@param degree polynomial degree of the underlying signal
#'@param scaling numeric calling constant in local tests pre-computed by `diffInf`
#'@param tacv_max_scale, maximum scale at which to calculate TVAC
#'@param b_max number of starting points to use for the global TAVC estimator; helps with robustness.
#'@references McGonigle, Euan T., and Haeran Cho. "Robust multiscale estimation of time-average variance for time series segmentation." Computational Statistics & Data Analysis 179 (2023): 107648. 
#'@references https://github.com/EuanMcGonigle/TAVC.seg
#'@export
gen_diff_catoni_tavc <- function(xx, ww, degree, tacv_max_scale = NULL, b_max = NULL)
{
  
  xx_cumsum <- c(0,cumsum(xx))
  
  scaling <- sum(choose(degree+1,0:(degree+1))**2) / (degree+2)
  
  nn <- length(xx_cumsum) - 1
  
  if (is.null(tacv_max_scale)) tacv_max_scale <- floor(2.5*sqrt(nn))
  
  if (ww > tacv_max_scale) ww <- tacv_max_scale
  
  if (!(ww %% (degree+2) == 0)) ww <- (degree+2) * floor(ww/(degree+2))
  
  if (is.null(b_max)) b_max <- ww
  
  
  tavc_est_vec <- numeric(b_max)
  
  for (ii in seq_along(tavc_est_vec))
  {
    zz <- sapply(seq(ii,nn,ww), function(ll) loc_diff(xx_cumsum, ll, ww, degree) / sqrt(ww*scaling)) ** 2
    
    zz <- na.omit(zz)
    
    qq_est <- sqrt(ww/(nn)) / (mean(zz, trim = 0.25))

    tavc_uniroot <- try(uniroot(catoni_sum, xx = zz, qq = qq_est, interval = c(-10,max(zz))), silent = TRUE)

    if(class(tavc_uniroot) == "try-error")
    {
      tavc_est_vec[ii] <- NA
    } else {
      tavc_est_vec[ii] = tavc_uniroot$root
    }
  }
  
  return(sqrt(median(tavc_est_vec, na.rm = TRUE)))
}


catoni_influence_fun <- function(xx)
{
  #' Influence function from Catoni (2012)
  #'
  #'@param xx data vector
  #'
  #'@references Catoni, Olivier. "Challenging the empirical mean and empirical variance: a deviation study." Annales de l'IHP Probabilités et statistiques. Vol. 48. No. 4. 2012.
  #'@references McGonigle, Euan T., and Haeran Cho. "Robust multiscale estimation of time-average variance for time series segmentation." Computational Statistics & Data Analysis 179 (2023): 107648. 
  #'@references https://github.com/EuanMcGonigle/TAVC.seg
  
  if (xx >= 1) return(log(2))
  
  if ((xx >= 0) && (xx < 1)) return(-log(1-xx+xx^2/2))
  
  if ((xx <= 0) && (xx > -1)) return(log(1+xx+xx^2/2))
  
  return(-log(2))
}

catoni_sum <- function(tau,xx,qq)
{
  #' Function to be optimized to obtain an estimate for \tau
  #'
  #'@param tau, parameter
  #'@param xx, data vector
  #'@param qq, scaling
  #'
  #'@references Catoni, Olivier. "Challenging the empirical mean and empirical variance: a deviation study." Annales de l'IHP Probabilités et statistiques. Vol. 48. No. 4. 2012.
  #'@references McGonigle, Euan T., and Haeran Cho. "Robust multiscale estimation of time-average variance for time series segmentation." Computational Statistics & Data Analysis 179 (2023): 107648. 
  #'@references https://github.com/EuanMcGonigle/TAVC.seg
  
  out <-sapply(xx, function(ii) (1/qq) * catoni_influence_fun(qq*(ii-tau)))
  
  return(mean(out))
}

#'Median Absolute Deviation estimator based on differenced data
#'@description TO DO! 
#'@param xx data
#'@param degree degree of the underlying signal
#'@export 
gen_diff_mad <- function(xx, degree)
{
  xx_diff <- diff(xx, differences = (degree+1))
  
  cp <- sum(choose(degree+1, 0:(degree+1))**2)
  
  return(mad(xx_diff/ sqrt(cp)))
}


#'Standard deviation estimator based on differenced data
#'@description TO DO!
#'@param xx data 
#'@param degree degree of the underlying 
#'@export
gen_diff_sd <- function(xx, degree)
{
  nn <- length(xx)
  
  cp <- sum(choose(degree+1, 0:(degree+1))**2)
  
  xx_diff <- diff(xx, differences = (degree+1))
  
  var_est <- sum((xx_diff**2)/cp) / (nn - (degree+1))
  
  return(sqrt(var_est))
}

