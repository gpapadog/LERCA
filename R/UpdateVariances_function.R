UpdateVariances <- function(dta, current_cutoffs, current_coefs, cov_cols,
                            alpha_priorX, beta_priorX, alpha_priorY,
                            beta_priorY) {
  
  K <- length(current_cutoffs)
  exact_cuts <- c(min(dta$X), current_cutoffs, max(dta$X))
  
  r <- array(0, dim = c(2, K + 1))
  
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, E == ee)
    n_k <- nrow(D)  # number of observations in current experiment.

    # For the exposure model.
    current_coefsX <- current_coefs[1, ee, - 2]
    des_mat <- matrix(1, nrow = n_k, ncol = 1)
    if (!is.null(cov_cols)) {
      des_mat <- cbind(des_mat, as.matrix(D[, cov_cols]))
    }
    resid <- D$X - des_mat %*% current_coefsX

    alpha_new <- alpha_priorX + n_k / 2
    beta_new <- beta_priorX + sum(resid ^ 2) / 2
    r[1, ee] <- invgamma::rinvgamma(1, shape = alpha_new, rate = beta_new)
    
    # For the outcome model.
    current_coefsY <- current_coefs[2, ee, ]
    des_mat <- cbind(1, D$X - exact_cuts[ee])
    if  (!is.null(cov_cols)) {
      des_mat <- cbind(des_mat, as.matrix(D[, cov_cols]))
    }
    resid <- D$Y - des_mat %*% current_coefsY
    
    alpha_new <- alpha_priorY + n_k / 2
    beta_new <- beta_priorY + sum(resid ^ 2) / 2
    r[2, ee] <- invgamma::rinvgamma(1, shape = alpha_new, rate = beta_new)
  }
  return(r)
}
  
