UpdateCoefficients <- function(dta, cov_cols, current_cutoffs,
                               current_alphas, current_vars,
                               Sigma_priorX, mu_priorX,
                               Sigma_priorY, mu_priorY) {
  
  minX <- min(dta$X) - 0.001
  maxX <- max(dta$X) + 0.001
  num_conf <- length(cov_cols)
  K <- length(current_cutoffs)
  
  r <- array(0, dim = c(2, K + 1, num_conf + 2))
  r[1, , 2] <- NA
  
  cuts <- c(minX, current_cutoffs, maxX)
  
  for (ee in 1 : (K + 1)) {
    
    D <- dta
    D <- subset(D, X >= cuts[ee] & X < cuts[ee + 1])

    # For the exposure model.
    
    current_alphaX <- current_alphas[1, ee, ]
    curr_variance <- current_vars[1, ee]
    
    which_in <- which(current_alphaX == 1)
    cov_cols_in <- cov_cols[which_in]
    C <- cbind(1, as.matrix(D[, cov_cols_in]))
    
    prior_var <- Sigma_priorX[c(1, which_in + 1), c(1, which_in + 1)]
    prior_mean <- mu_priorX[c(1, which_in + 1)]

    prior_var_inv <- chol2inv(chol(prior_var))
    post_var <- t(C) %*% C / curr_variance + prior_var_inv
    post_var <- chol2inv(chol(post_var))
    
    post_mean <- t(C) %*% matrix(D$X, ncol = 1) / curr_variance
    post_mean <- post_var %*% (post_mean + prior_var_inv %*% prior_mean)
    
    gen_coef <- mvnfast::rmvn(1, mu = post_mean, sigma = post_var)
    r[1, ee, c(1, which_in + 2)] <- gen_coef
    
    
    # For the outcome model.
    current_alphaY <- current_alphas[2, ee, ]
    curr_variance <- current_vars[2, ee]
    
    which_in <- which(current_alphaY == 1)
    cov_cols_in <- cov_cols[which_in]
    
    prior_var <- Sigma_priorY[c(1, 2, which_in + 2), c(1, 2, which_in + 2)]
    prior_mean <- mu_priorY[c(1, 2, which_in + 2)]
    C <- cbind(1, X = D$X, as.matrix(D[, cov_cols_in]))
    
    prior_var_inv <- chol2inv(chol(prior_var))
    post_var <- t(C) %*% C / curr_variance + prior_var_inv
    post_var <- chol2inv(chol(post_var))
    
    post_mean <- t(C) %*% matrix(D$Y, ncol = 1) / curr_variance
    post_mean <- post_var %*% (post_mean + prior_var_inv %*% prior_mean)
    
    gen_coef <- mvnfast::rmvn(1, mu = post_mean, sigma = post_var)
    r[2, ee, c(1, 2, which_in + 2)] <- gen_coef
    
  }
  return(r)
}
