#'@title Inference for Change Points in Piecewise Polynomials via Differencing
#'@description Identifies sub-intervals of the data sequence which each must contain a change point, in the sense that the mean function cannot be described as a polynomial of degree \eqn{p}, uniformly with probability asymptotically larger than \eqn{1 - \alpha + o(1)} where \eqn{\alpha \in (0,1)} can be set by adjusting the parameter \code{alpha}.
#' The object returned by this function can be passed to:
#' \itemize{
#'  \item \code{\link{plot}}: for potting the data along with the intervals of significance returned. 
#'  \item \code{\link{cpt}}: for estimating the most probable change point location within each interval recovered. 
#'  \item \code{\link{predict}}: for estimating the unobserved piecewise polynomial mean function, using the most likely change point locations recovered by \code{cpt}.
#' }
#'@details The data supplied in \code{yy} are assumed to follow the 'signal + noise' model: \deqn{Y_t = f_\circ (t/n) + \zeta_t \hspace{2em} t = 1, \dots, n.} Where \eqn{n} is the length of the data sequence, \eqn{f_\circ (t/n)} is a piecewise polynomial function of know degree \eqn{p}, and the noise terms are one of: 
#' \enumerate{
#' \item Independently distributed and Gaussian, with common variance. 
#' \item Independently distributed with common variance, but not necessarily Gaussian.
#' \item Weakly stationary with strictly positive long run variance.
#' }
#'@param yy A numeric vector containing the data to be inspected for change points. 
#'@param degree The degree of polynomial parametrization of the mean of the data. 
#'@param alpha Desired maximum probability of obtaining an interval that does not contain a change point. 
#'@param gaussian_noise Set to \code{TRUE} if the contaminating noise noise is assumed to be independently distributed and Gaussian, with common variance, else set to \code{FALSE}. 
#'@param independent_noise Set to \code{TRUE} if the contaminating noise is assumed to be independently distributed, else set to \code{FALSE}. 
#'@param tau Noise level in as measured by the (long run) standard deviation of the noise, if know. If \code{tau=NULL} the noise level will be estimated via: 
#'\itemize{
#'  \item \link{mad_diff} if \code{gaussian_noise = TRUE} and \code{independent_noise = TRUE}. 
#'  \item \link{sd_diff} if \code{gaussian_noise = FALSE} and \code{independent_noise = TRUE}. 
#'  \item \link{lrsd_block_diff} if \code{independent_noise = FALSE}. 
#'}
#'@param aa Decay parameter controlling the density of the grid on which local tests for the presence of a change point will be performed. Tests will be performed on all contiguous sub-intervals of \eqn{\left \{ 1, \dots, n \right \}} whose length is larger than \code{min_scale} and can be expressed as an integer power of \code{aa}. 
#'@param min_scale Minimum scale at which local tests for the presence of a change point will be performed. 
#'@param HH Numeric constant appearing in the extreme value limit of the supremum over all local tests, under the null of no change points. If pre-computed the value can be passed to the function. Otherwise, the value is computed automatically using. 
#'@return An object of class "list" and "not", which contains the following fields:
#'\itemize{
#'  \item \code{intervals}: intervals of significance returned by the procedure. 
#'  \item \code{thresh}: the threshold used to for each local test. 
#'  \item \code{data}: the original input data \code{yy}. 
#'  \item \code{degree}: degree of polynomial parametrization of the mean of the data. 
#'}
#'@examples
#'## Piecewise constant mean with i.i.d. Gaussian noise
#'
#'set.seed(42)
#'
#'blocks_signal <- c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40))
#'
#'yy <- blocks_signal + rnorm(length(blocks_signal), sd = 5)
#'
#'
#'## Recover intervals of significance
#'
#'diffInf_obj <- diffInf(yy, degree = 0)
#'
#'diffInf_obj
#'
#'
#'## some examples of how the `diffInf` object can be used
#'
#'# plot the intervals of significance
#'
#'diffInf_obj |> plot(type = "l", col = "grey")
#'
#'blocks_signal |> lines(lty = 2, lwd = 2)
#'
#'
#'# recover likely change point locations
#'
#'cpt(diffInf_obj)
#'
#'abline(v = cpt(diffInf_obj), col = "red", lty = 2)
#'
#'
#'# estimate the best piecewise polynomial fit based on the change point locations
#'
#'predict(diffInf_obj)
#'
#'lines(predict(diffInf_obj), col = "red")
#'
#'@export  

diffInf <- function(yy, degree, alpha = 0.1, gaussian_noise = TRUE, independent_noise = TRUE, tau = NULL, aa = sqrt(2), min_scale = floor(sqrt(length(yy))/2), HH = NULL)
{
  nn <- length(yy)
  
  yy_cumsum <- c(0,cumsum(yy))
  
  scaling <- sum(choose(degree+1,0:(degree+1))**2) / (degree+2)
  
  
  noise_level_know <- !is.null(tau)
  
  if (!noise_level_know)
  {
    if (gaussian_noise && independent_noise) tau <- mad_diff(yy, degree)
    
    if (!gaussian_noise && independent_noise) tau <- sd_diff(yy, degree)
    
    if (!independent_noise) tau <- lrsd_block_diff(yy, min_scale, degree)
  }
  
  
  if (gaussian_noise && independent_noise) min_scale <- log(nn)
  
  thresh <- tau * get_thresh(nn, min_scale, alpha, degree, aa, HH, (gaussian_noise && independent_noise))
  
  ints_df <- data.frame(matrix(nrow = 0, ncol = 3))
  
  names(ints_df) <- c("start", "end", "value")
  
  ints_out <- diff_bin_seg(yy_cumsum, 1, nn, degree, aa, min_scale, thresh, scaling)

  ints_out <- rbind(ints_df,ints_out)
  
  
  out <- list(intervals = ints_out, thresh = thresh, data = yy, degree = degree)
  
  class(out) <- c("cptInference", class(out))
  
  return(out)
}
