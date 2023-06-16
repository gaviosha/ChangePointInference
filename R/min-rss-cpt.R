
min_rss_cpt <- function(data, ss, ee, degree)
{
  if (ee - ss < (degree+1)) return(round(mean(c(ss,ee))))
  
  cpt_locs <- (ss + degree + 1):(ee - degree - 1)
  
  rss <- numeric(length(cpt_locs))
  
  for (ii in seq_along(cpt_locs))
  {
    if (degree == 0) rss[ii] <- sum((data[ss:cpt_locs[ii]] - mean(data[ss:cpt_locs[ii]]))**2) + sum((data[(cpt_locs[ii]+1):ee] - mean(data[(cpt_locs[ii]+1):ee]))**2)
    
    if (degree > 0) rss[ii] <- deviance(lm(data[ss:cpt_locs[ii]] ~ poly(ss:cpt_locs[ii], degree))) + deviance(lm(data[(cpt_locs[ii]+1):ee] ~ poly((cpt_locs[ii]+1):ee, degree)))
  }
  
  return(cpt_locs[which.min(rss)])
}