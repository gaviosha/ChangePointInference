
p_inf_sqrd <- function(xx)
{
  #' P-inf
  #' 
  #' Numerically approximates the square of function p_\inf on page 32 of Kabluchko (2018)
  #'
  #'@param xx a positive number
  #'
  #'@references Kabluchko, Zakhar. "Extreme-value analysis of standardized Gaussian increments." arXiv preprint arXiv:0706.1849 (2018).
  #'
  #'@export
  
  exp(-2*sum((1/(1:1000)) * pnorm(-sqrt((xx*(1:1000)/4)))))
}


H_num_est <- function(degree, aa, LL)
{
  #' H-num-est
  #'
  #' Numeric estimate of the constant H_2 in the paper
  #'
  #'@param degree polynomial degree of the underlying signal
  #'@param aa controls grid density (a in the paper)
  #'@param LL control smallest scale (L in the paper)
  #'
  #'@export
  
  c1 <- sum(choose(degree+1,1:(degree+1))*choose(degree+1,0:degree))
  
  c2 <- sum(choose(degree+1,0:(degree+1))**2)
  
  C_xi_p <- (degree+2)*(1+(c1/c2))
  
  
  up_seq <- sapply(-LL:100, function(jj) p_inf_sqrd((2*1*C_xi_p)/(aa**jj)))
  
  return(sum(up_seq))
}

get_thresh <- function(nn, alpha, degree, aa, LL, HH)
{
  #' Get threshold
  #'
  #' Calculate threshold for FWE controll 
  #'
  #'@param nn length of the data sequence
  #'@param degree polynomial degree of the underlying signal
  #'@param aa controls grid density (a in the paper)
  #'@param LL control smallest scale (L in the paper)
  #'@param HH numeric constant 
  #'
  #'@export 
  
  HH <- ifelse(is.null(HH), H_num_est(degree, aa, LL), HH)
  
  a_n <- sqrt(2*log(nn)) + (-0.5*log(log(nn)) - log(2*sqrt(pi)) + log(HH)) / sqrt(2*log(nn))
  
  b_n <- 1/sqrt(2*log(nn))
  
  tau <- log(1/log(1/sqrt(1-alpha)))
  
  return(a_n + tau * b_n)
}
