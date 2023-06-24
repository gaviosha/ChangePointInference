#' @title Plot a 'cptInference' Object
#' @description Plots the data and as well as intervals of significance returned by the \code{\link{diffInf}}.
#' @param obj An object of class 'cptInference'. 
#' @param cpt_loc_est One of \code{"RSS"}, \code{"midpoint"}, or \code{FALSE}. If \code{"RSS"} or \code{"midpoint"} the mosyt likely change pint location within each interval will be plotted. If \code{FALSE} only the intervals are plotted. 
#' @param ... Additional graphical parameters passed to \code{\link{plot}}. 
#' @method plot cptInference
#' @examples
#' # Piecewise linear mean with i.i.d. Gaussian noise
#' 
#' set.seed(42)
#' 
#' waves_signal <- c((1:150) * (2**-3), (150:1) * (2**-3), (1:150) * (2**-3), (150:1) * (2**-3))
#' 
#' yy <- waves_signal + rnorm(length(waves_signal), sd = 5)
#' 
#' 
#' # Recover intervals of significance
#' 
#' diffInf_obj <- diffInf(yy, degree = 1)
#' 
#' diffInf_obj
#' 
#' # plot the intervals of significance
#' 
#' diffInf_obj |> plot("RSS", type = "l", col = "grey")
#' 
#' # plot intervals and minimum RSS split points
#' 
#' diffInf_obj |> plot("RSS", type = "l", col = "grey")
#' 
#' # plot the intervals and their midpoints
#' 
#' diffInf_obj |> plot("midpoint", type = "l", col = "grey")
#' @export 
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
