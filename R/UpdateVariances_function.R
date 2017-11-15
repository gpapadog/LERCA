UpdateVariances <- function(dta, current_cutoffs, current_coefs,
                            alpha_priorX, beta_priorX, alpha_priorY,
                            beta_priorY) {
  
  minX <- min(dta$X) - 0.001
  maxX <- max(dta$X) + 0.001
  K <- length(current_cutoffs)
  cuts <- c(minX, current_cutoffs, maxX)
  
  r <- array(0, dim = c(2, K + 1))
  
  for (ee in 1 : (K + 1)) {
    
    D <- dta
    D <- subset(D, X >= cuts[ee] & X < cuts[ee + 1])
    N_ee <- nrow(D)  # number of observations in current experiment.

    # For the exposure model.
    current_coefsX <- current_coefs[1, ee, - 2]
    resid <- D$X - cbind(1, as.matrix(D[, cov_cols])) %*% current_coefsX
    
    alpha_new <- alpha_priorX + N_ee / 2
    beta_new <- beta_priorX + sum(resid ^ 2) / 2
    r[1, ee] <- invgamma::rinvgamma(1, shape = alpha_new, rate = beta_new)
    
    # For the outcome model.
    current_coefsY <- current_coefs[2, ee, ]
    resid <- D$Y - cbind(1, D$X, as.matrix(D[, cov_cols])) %*% current_coefsY
    
    alpha_new <- alpha_priorY + N_ee / 2
    beta_new <- beta_priorY + sum(resid ^ 2) / 2
    r[2, ee] <- invgamma::rinvgamma(1, shape = alpha_new, rate = beta_new)
  }
  return(r)
}
  
