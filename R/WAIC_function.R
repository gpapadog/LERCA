#' WAIC of LERCA fit using the full conditionals.
#'
#' @param lerca The LERCA output including cutoffs, coefficients and variances
#' of the exposure and outcome models.
#' @param dta The full data set including the confounders as C1, C2, ..., the
#' outcome as Y and the exposure as X.
#' 
#' @return Numeric vector with entries: the estimated log predictive density,
#' the penalty term based on both criteria, and the WAIC based on both penalty
#' criteria.
#' 
#' @export
WAIC <- function(lerca, dta) {
  
  dta <- as.data.frame(dta)
  N <- nrow(dta)
  num_conf <- dim(lerca$alphas)[5]
  cov_cols <- which((names(dta) %in% paste0('C', 1 : num_conf)))
  post_samples <- dim(lerca$cutoffs)[2]
  chains <- dim(lerca$cutoffs)[1]
  K <- dim(lerca$cutoffs)[3]
  
  des_mat <- as.matrix(cbind(Int = 1, X = dta$X, dta[, cov_cols]))
  
  progress <- floor(seq(0, post_samples, length.out = 11)[- 1])
  likelihood <- array(NA, dim = c(chains, post_samples, N))
  
  for (ss in 1 : post_samples) {
    
    if (ss %in% progress) {
      print(paste0(which(progress == ss) * 10, '% completed.'))
    }
    
    for (cc in 1 : chains) {
      
      curr_cutoffs <- lerca$cutoffs[cc, ss, ]
      dta_exper <- sapply(dta$X, function(x) sum(curr_cutoffs < x) + 1)
      
      for (ee in 1 : (K + 1)) {
        
        wh_obs <- which(dta_exper == ee)
        curr_coefsX <- lerca$coefs[1, cc, ss, ee, - 2]
        curr_coefsY <- lerca$coefs[2, cc, ss, ee, ]
        curr_sdX <- sqrt(lerca$variances[1, cc, ss, ee])
        curr_sdY <- sqrt(lerca$variances[2, cc, ss, ee])
        curr_meanX <- des_mat[wh_obs, - 2] %*% curr_coefsX
        curr_meanY <- des_mat[wh_obs, ] %*% curr_coefsY
        
        likeY <- dnorm(x = dta$Y[wh_obs], sd = curr_sdY, mean = curr_meanY)
        likeX <- dnorm(x = dta$X[wh_obs], sd = curr_sdX, mean = curr_meanX)
        
        likelihood[cc, ss, wh_obs] <- likeY * likeX
      }
    }
  }
  
  logE <- log(apply(likelihood, 3, mean))
  Elog <- apply(log(likelihood), 3, mean)
  
  pwaic1 <- 2 * sum(logE - Elog)
  pwaic2 <- sum(apply(log(likelihood), 3, function(x) var(as.numeric(x))))
  lppd <- sum(logE)
  
  r <- c(lppd = lppd, pwaic1 = pwaic1, pwaic2 = pwaic2,
         waic1 = - 2 * (lppd - pwaic1), waic2 = - 2 * (lppd - pwaic2))
  return(r)
}
