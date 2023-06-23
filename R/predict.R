#' @title Predict a 'cptInference' object
#' @description 123
#' @method predict cptInference
#' @export 
#' @param obj An object of class 'cptInference', returned by \code{diffInf}.

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







