
p_inf_sqrd <- function(yy)
{
  #' P-inf
  #' 
  #' Numerically approximates the square of function p_\inf on page 32 of Kabluchko (2018)
  #'
  #'@param yy a positive number
  #'
  #'@references Kabluchko, Zakhar. "Extreme-value analysis of standardized Gaussian increments." arXiv preprint arXiv:0706.1849 (2018).
  #'@noRd
  
  exp(-2*sum((1/(1:1000)) * pnorm(-sqrt((yy*(1:1000)/4)))))
}


H_num_est_gauss <- function(degree, aa)
{
  #' Estimate the numeric constant H under Gaussian noise
  #'
  #' Numeric estimate of the constant H_2 in the paper
  #'
  #'@param degree polynomial degree of the underlying signal
  #'@param aa controls grid density (a in the paper)
  #'@noRd
  
  c1 <- sum(choose(degree+1,1:(degree+1))*choose(degree+1,0:degree))
  
  c2 <- sum(choose(degree+1,0:(degree+1))**2)
  
  C_xi_p <- (degree+2)*(1+(c1/c2))
  
  
  up_seq <- sapply(0:100, function(jj) p_inf_sqrd((2*1*C_xi_p)/(aa**jj)))
  
  return(sum(up_seq))
}


H_nongauss <- function(degree, aa)
{
  #' Get the numeric constant H under non-Gaussian and / or dependent noise
  #' 
  #' Calculates the constant H_4 in the paper
  #'
  #'@param degree polynomial degree of the underlying signal
  #'@param aa controls grid density (a in the paper)
  #'@noRd
  
  c1 <- sum(choose(degree+1,1:(degree+1))*choose(degree+1,0:degree))
  
  c2 <- sum(choose(degree+1,0:(degree+1))**2)
  
  C_xi_p <- (degree+2)*(1+(c1/c2))
  
  return(C_xi_p / (1-(1/aa)))
}

get_thresh <- function(nn, WW, alpha, degree, aa, HH, gauss_indep_noise)
{
  #' Get threshold
  #'
  #' Calculate threshold for FWE control 
  #'
  #'@param nn length of the data sequence
  #'@param WW minimum segment length
  #'@param alpha desired coverage 
  #'@param degree polynomial degree of the underlying signal
  #'@param aa controls grid density (a in the paper)
  #'@param HH numeric constant 
  #'@param gauss_indep_noise true if noise is Gaussian and independent
  #'@noRd
  
  if (gauss_indep_noise)
  {
    HH <- ifelse(is.null(HH), H_num_est_gauss(degree, aa), HH)
    
    a_n <- sqrt(2*log(nn)) + (-0.5*log(log(nn)) - log(2*sqrt(pi)) + log(HH)) / sqrt(2*log(nn))
    
    b_n <- 1/sqrt(2*log(nn))
    
  } else {
    
    HH <- ifelse(is.null(HH), H_nongauss(degree, aa), HH)
    
    a_n <- sqrt(2*log(nn/WW)) + (0.5*log(log(nn/WW)) - log(sqrt(pi)) + log(HH)) / sqrt(2*log(nn/WW))
    
    b_n <- 1/sqrt(2*log(nn/WW))
  }
    
  tau <- log(1/log(1/sqrt(1-alpha)))
  
  return(a_n + tau * b_n)
}
