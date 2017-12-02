WBIC <- function(lerca, dta) {
  
  dta <- as.data.frame(dta)
  N <- nrow(dta)
  num_conf <- dim(lerca$alphas)[5]
  covs_cols <- which((names(dta) %in% paste0('C', 1 : num_conf)))
  post_samples <- dim(lerca$cutoffs)[2]
  chains <- dim(lerca$cutoffs)[1]

  des_mat <- cbind(Int = 1, X = dta$X, dta[, cov_cols])
  
  wbic <- 0
  
  for (ii in 1 : N) {
    for (ss in 1 : post_samples) {
      for (cc in 1 : chains) {
        
        curr_exper <- sum(lerca$cutoffs[cc, ss, ] < dta$X[ii]) + 1
        curr_coefsX <- lerca$coefs[1, cc, ss, curr_exper, - 2]
        curr_coefsY <- lerca$coefs[2, cc, ss, curr_exper, ]
        curr_sdX <- sqrt(lerca$variances[1, cc, ss, curr_exper])
        curr_sdY <- sqrt(lerca$variances[2, cc, ss, curr_exper])
        
        wbic <- wbic +
          dnorm(x = dta$Y[ii], mean = sum(des_mat[ii, ] * curr_coefsY),
                sd = curr_sdY, log = TRUE) + 
          dnorm(x = dta$X[ii], mean = sum(des_mat[ii, - 2] * curr_coefsX),
                sd = curr_sdX, log = TRUE)
      }
    }
  }
  
  wbic <- - wbic / (post_samples * chains)
  
  return(wbic)
}