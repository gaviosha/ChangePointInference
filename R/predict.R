#' @title Predict a 'cptInference' Object
#' @description Fits a piecewise polynomial signal to the to the data, using change point locations recovered by \code{\link{cpt}}
#' @param obj An object of class 'cptInference'. 
#' @param cpt_loc Determines how the change points will be estimated. If \code{cpt_loc = "RSS"} the change points are taken to be the split points within each interval which provide the lowest RSS when a piecewise polynomial function is fitted on the same interval. If \code{cpt_loc = "midpoint"} the change points are taken to be the midpoint of each interval. 
#' @method predict cptInference
#' @examples
#' #' Piecewise linear mean with i.i.d. Gaussian noise
#' 
#' set.seed(42)
#' 
#' waves_signal <- c((1:150) * (2**-3), (150:1) * (2**-3), (1:150) * (2**-3), (150:1) * (2**-3))
#' 
#' yy <- waves_signal + rnorm(length(waves_signal), sd = 5)
#' 
#' 
#' #' Recover intervals of significance
#' 
#' diffInf_obj <- diffInf(yy, degree = 1)
#' 
#' diffInf_obj
#' 
#' #' plot the intervals of significance
#' 
#' diffInf_obj |> plot(type = "l", col = "grey")
#' 
#' waves_signal |> lines(lty = 2, lwd = 2)
#' 
#' 
#' #' recover fitted signal
#' 
#' lines(predict(diffInf_obj, "RSS"), col = "red", lty = 2, lwd = 2)
#' 
#' lines(predict(diffInf_obj, "midpoint"), col = "blue", lty = 3, lwd = 2)
#' 
#' @export 
predict.cptInference <- function(obj, cpt_loc = "RSS")
{
  data <- obj$data
  
  nn <- length(data)
  
  fit <- c()
  
  degree <- obj$degree
  
  intervals <- obj$intervals
  
  
  if ((nrow(intervals) == 0) && (degree == 0)) return(list(fit = rep(mean(data),nn), est_cpt_loc = NULL))
  
  if ((nrow(intervals) == 0) && (degree > 0)) return(list(fit = predict(lm(data ~ poly(1:nn, degree = degree))), est_cpt_loc = NULL))
  
  est_cpt_locs <- cpt(obj, cpt_loc)
  
  NN <- length(est_cpt_locs)
  
  est_cpt_locs <- unique(c(0, est_cpt_locs, nn))
  
  for (ii in 1:(NN+1))
  {
    ss <- est_cpt_locs[ii] + 1
    ee <- est_cpt_locs[ii+1]
    
    if (degree == 0) fit <- c(fit, rep(mean(data[ss:ee]), ee - ss + 1))
    if (degree > 0) fit <- c(fit, predict(lm(data[ss:ee] ~ poly(ss:ee, degree = degree))))
  }
  
  return(fit)
}







