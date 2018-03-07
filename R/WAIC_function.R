#' WAIC of LERCA fit using the full conditionals.
#'
#' @param lerca The LERCA output including cutoffs, coefficients and variances
#' of the exposure and outcome models.
#' @param dta The full data set including the confounders as C1, C2, ..., the
#' outcome as Y and the exposure as X.
#' 
#' @export
WAIC <- function(lerca, dta) {
  
  dta <- as.data.frame(dta)
  N <- nrow(dta)
  num_conf <- dim(lerca$alphas)[5]
  covs_cols <- which((names(dta) %in% paste0('C', 1 : num_conf)))
  post_samples <- dim(lerca$cutoffs)[2]
  chains <- dim(lerca$cutoffs)[1]
  
  des_mat <- cbind(Int = 1, X = dta$X, dta[, cov_cols])

  lppd <- 0
  Elog <- 0
  logE <- 0
  pwaic1 <- 0
  pwaic2 <- 0
  
  for (ii in 1 : N) {
    
    like <- matrix(NA, nrow = post_samples, ncol = chains)
    
    for (ss in 1 : post_samples) {
      for (cc in 1 : chains) {
        
        curr_exper <- sum(lerca$cutoffs[cc, ss, ] < dta$X[ii]) + 1
        curr_coefsX <- lerca$coefs[1, cc, ss, curr_exper, - 2]
        curr_coefsY <- lerca$coefs[2, cc, ss, curr_exper, ]
        curr_sdX <- sqrt(lerca$variances[1, cc, ss, curr_exper])
        curr_sdY <- sqrt(lerca$variances[2, cc, ss, curr_exper])
        
        like[ss, cc] <- dnorm(x = dta$Y[ii], sd = curr_sdY,
                              mean = sum(des_mat[ii, ] * curr_coefsY)) *
          dnorm(x = dta$X[ii], mean = sum(des_mat[ii, - 2] * curr_coefsX),
                sd = curr_sdX)
      }
    }
    
    logE_ii <- log(mean(like))
    Elog_ii <- mean(log(like))
    
    lppd <- lppd + logE_ii
    pwaic1 <- pwaic1 + 2 * (logE_ii - Elog_ii)
    pwaic2 <- pwaic2 + var(as.numeric(log(like)))
  }
  r <- c(lppd = lppd, pwaic1 = pwaic1, pwaic2 = pwaic2,
         waic1 = - 2 * (lppd - pwaic1), waic2 = - 2 * (lppd - pwaic2))
  return(r)
}
