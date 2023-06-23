#' #'@title Robust estimation of TAVC with unknown piecewise polynomial signal
#' #'@description A short description...
#' #'@param yy_cumus, cumulative sum of the data 
#' #'@param ww, scale at which TVAC will be calculated  
#' #'@param degree polynomial degree of the underlying signal
#' #'@param scaling numeric calling constant in local tests pre-computed by `diffInf`
#' #'@param tacv_max_scale, maximum scale at which to calculate TVAC
#' #'@param b_max number of starting points to use for the global TAVC estimator; helps with robustness.
#' #'@references McGonigle, Euan T., and Haeran Cho. "Robust multiscale estimation of time-average variance for time series segmentation." Computational Statistics & Data Analysis 179 (2023): 107648. 
#' #'@references https://github.com/EuanMcGonigle/TAVC.seg
#' #'@export
#' gen_diff_catoni_tavc <- function(yy, ww, degree, tacv_max_scale = NULL, b_max = NULL)
#' {
#'   
#'   yy_cumsum <- c(0,cumsum(yy))
#'   
#'   scaling <- sum(choose(degree+1,0:(degree+1))**2) / (degree+2)
#'   
#'   nn <- length(yy_cumsum) - 1
#'   
#'   if (is.null(tacv_max_scale)) tacv_max_scale <- floor(2.5*sqrt(nn))
#'   
#'   if (ww > tacv_max_scale) ww <- tacv_max_scale
#'   
#'   if (!(ww %% (degree+2) == 0)) ww <- (degree+2) * floor(ww/(degree+2))
#'   
#'   if (is.null(b_max)) b_max <- ww
#'   
#'   
#'   tavc_est_vec <- numeric(b_max)
#'   
#'   for (ii in seq_along(tavc_est_vec))
#'   {
#'     zz <- sapply(seq(ii,nn,ww), function(ll) loc_diff(yy_cumsum, ll, ww, degree) / sqrt(ww*scaling)) ** 2
#'     
#'     zz <- na.omit(zz)
#'     
#'     qq_est <- sqrt(ww/(nn)) / (mean(zz, trim = 0.25))
#' 
#'     tavc_uniroot <- try(uniroot(catoni_sum, yy = zz, qq = qq_est, interval = c(-10,max(zz))), silent = TRUE)
#' 
#'     if(class(tavc_uniroot) == "try-error")
#'     {
#'       tavc_est_vec[ii] <- NA
#'     } else {
#'       tavc_est_vec[ii] = tavc_uniroot$root
#'     }
#'   }
#'   
#'   return(sqrt(median(tavc_est_vec, na.rm = TRUE)))
#' }
#' 
#' 
#' catoni_influence_fun <- function(yy)
#' {
#'   #' Influence function from Catoni (2012)
#'   #'
#'   #'@param yy data vector
#'   #'
#'   #'@references Catoni, Olivier. "Challenging the empirical mean and empirical variance: a deviation study." Annales de l'IHP Probabilités et statistiques. Vol. 48. No. 4. 2012.
#'   #'@references McGonigle, Euan T., and Haeran Cho. "Robust multiscale estimation of time-average variance for time series segmentation." Computational Statistics & Data Analysis 179 (2023): 107648. 
#'   #'@references https://github.com/EuanMcGonigle/TAVC.seg
#'   
#'   if (yy >= 1) return(log(2))
#'   
#'   if ((yy >= 0) && (yy < 1)) return(-log(1-yy+yy^2/2))
#'   
#'   if ((yy <= 0) && (yy > -1)) return(log(1+yy+yy^2/2))
#'   
#'   return(-log(2))
#' }
#' 
#' catoni_sum <- function(tau,yy,qq)
#' {
#'   #' Function to be optimized to obtain an estimate for \tau
#'   #'
#'   #'@param tau, parameter
#'   #'@param yy, data vector
#'   #'@param qq, scaling
#'   #'
#'   #'@references Catoni, Olivier. "Challenging the empirical mean and empirical variance: a deviation study." Annales de l'IHP Probabilités et statistiques. Vol. 48. No. 4. 2012.
#'   #'@references McGonigle, Euan T., and Haeran Cho. "Robust multiscale estimation of time-average variance for time series segmentation." Computational Statistics & Data Analysis 179 (2023): 107648. 
#'   #'@references https://github.com/EuanMcGonigle/TAVC.seg
#'   
#'   out <-sapply(yy, function(ii) (1/qq) * catoni_influence_fun(qq*(ii-tau)))
#'   
#'   return(mean(out))

