#' @title Change point locations
#' @description Estimate change point locations 
#' @export 
#' @param obj An object of class 'cptInference', returned by \code{diffInf}.

cpt <- function(obj, cpt_loc = c("midpoint","RSS")[1])
{
  data <- obj$data
  
  nn <- length(data)
  
  degree <- obj$degree
  
  intervals <- obj$intervals
  
  if (cpt_loc == "midpoint") est_cpt_locs <- apply(intervals, 1, function(ii) round((ii[1]+ii[2])/2))
  
  if (cpt_loc == "RSS") est_cpt_locs <- apply(intervals, 1, function(ii) min_rss_cpt(data,ii[1],ii[2], degree))
  
  return(est_cpt_locs)
}