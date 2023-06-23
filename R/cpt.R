#' @title Estimate Change Point Locations
#' @description Estimates the most likely change point location within each of the intervals obtained by calling \code{\link{diffInf}}.
#' @param obj An object of class \code{cptInference} returned by \code{\link{diffInf}}.
#' @param cpt_loc Determines how the change point locations will be estimated. 
#' Should be one of "midpoint" or "RSS". 
#' If \code{cpt_loc == "midpoint"} the change point locations are simply taken to be the midpoints of each interval. 
#' If \code{cpt_loc == "RSS"} the change point locations are chosen to be the split points within each interval associated with the piecewsie polynomial fit providing the lowest sum of squared residuals.
#' @returns A vector of estimated change point locations. 
#' @examples 
#' # Piecewise linear mean with i.i.d. Gaussian noise
#' set.seed(42)
#'  
#' waves_signal <- c((1:150) * (2**-3), (150:1) * (2**-3), (1:150) * (2**-3), (150:1) * (2**-3))
#'  
#' yy <- waves_signal + rnorm(length(waves_signal), sd = 5)
#'  
#'  
#' #' Recover intervals of significance
#'  
#'  diffInf_obj <- diffInf(yy, degree = 1)
#'  
#'  diffInf_obj
#'  
#' #' plot the intervals of significance
#'  
#'  diffInf_obj |> plot(type = "l", col = "grey")
#'  
#'  waves_signal |> lines(lty = 2, lwd = 2)
#'  
#'  
#' #' recover likely change point locations
#'  
#' cpt(diffInf_obj, "RSS")
#'  
#' cpt(diffInf_obj, "midpoint")
#'  
#' abline(v = cpt(diffInf_obj, "RSS"), col = "red", lty = 2)
#'  
#' abline(v = cpt(diffInf_obj, "midpoint"), col = "blue", lty = 3)
#' @export 

cpt <- function(obj, cpt_loc = "RSS")
{
  intervals <- obj$intervals
  
  if (nrow(intervals) == 0) return(NULL)
  
  data <- obj$data
  
  nn <- length(data)
  
  degree <- obj$degree
  
  if (cpt_loc == "midpoint") est_cpt_locs <- apply(intervals, 1, function(ii) round((ii[1]+ii[2])/2))
  
  if (cpt_loc == "RSS") est_cpt_locs <- apply(intervals, 1, function(ii) min_rss_cpt(data, ii[1],ii[2], degree))
  
  return(est_cpt_locs)
}