#' @title Plot a 'cptInference' object
#' @description Plots the data and as well as intervals of significance returned by the algorithm
#' @method plot cptInference
#' @export 
#' @param obj An object of class 'cptInference', returned by \code{diffInf}.

plot.cptInference <- function(obj, cpt_loc = c(FALSE,"midpoint")[1], ...)
{
  data <- obj$data
  
  regions <- obj$intervals
  
  plot(data, ...)
  
  for (ii in 1:nrow(regions))
  {
    rect(
      xleft = regions[ii,1],
      ybottom = min(data),
      xright =  regions[ii,2],
      ytop = max(data),
      col = rgb(1,0,0,.2),
      border = NA
    )
    
    if (is.character(cpt_loc))
    {
      if (cpt_loc == "midpoint")
      {
        abline(
          v = (regions[ii,1] + regions[ii,2])/2,
          lwd = 1.5,
          lty = 3
        )
      }
      
      if (cpt_loc == "RSS")
      {
        
      }
    }
  }
}
