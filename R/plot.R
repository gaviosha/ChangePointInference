#' @title Plot a 'cptInference' object
#' @description Plots the data and as well as intervals of significance returned by the algorithm
#' @method plot cptInference
#' @export 
#' @param obj An object of class 'cptInference', returned by \code{diffInf}.

plot.cptInference <- function(obj, cpt_loc_est = FALSE, ...)
{
  data <- obj$data
  
  regions <- obj$intervals
  
  plot(data, ...)
  
  if (nrow(regions) > 0)
  {
    for (ii in 1:nrow(regions))
    {
      rect(
        xleft = regions[ii,1],
        ybottom = -2*min(data)*sign(min(data)),
        xright =  regions[ii,2],
        ytop = 2*max(data)*sign(max(data)),
        col = rgb(1,0,0,.2),
        border = NA
      )
      
      if (is.character(cpt_loc_est))
      {
        if (!(cpt_loc_est %in% c("midpoint","RSS"))) stop("cpt_loc_est should be one of 'midpoint' or 'RSS'")
        
        cpt_locs <- cpt(obj, cpt_loc_est)
        
        abline(
          v = cpt_locs,
          lwd = 1.5,
          lty = 3,
          col = "red"
        )
      }
    }
  }
}
