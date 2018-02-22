UpdateCovCoef <- function(dta, cov_cols, current_cutoffs, current_coefs,
                          current_alphas, current_vars, Sigma_priorX,
                          mu_priorX, Sigma_priorY, mu_priorY) {
  
  num_conf <- ifelse(is.null(cov_cols), 0, length(cov_cols))
  K <- length(current_cutoffs)
  exact_cuts <- c(min(dta$X), current_cutoffs, max(dta$X))
  
  # Coefficients of the covariates including intercept (of exposure only).
  r <- array(0, dim = c(2, K + 1, num_conf + 1))
  cov_names <- 'Int'
  if (num_conf > 0) {
    cov_names <- c(cov_names, paste0('C', 1 : num_conf))
  }
  dimnames(r) <- list(model = c('Exposure', 'Outcome'),
                      exper = 1 : (K + 1), covar = cov_names)
  r[2, , 1] <- NA  # Intercept of the outcome model is not udpated here.
  
  
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, E == ee)
    
    # For the exposure model.
    curr_variance <- current_vars[1, ee]
    C <- matrix(1, nrow = nrow(D), ncol = 1)
    which_in <- NULL
    
    if (num_conf > 0) {
      current_alphaX <- current_alphas[1, ee, ]
      which_in <- which(current_alphaX == 1)
      cov_cols_in <- cov_cols[which_in]
      C <- cbind(C, as.matrix(D[, cov_cols_in]))
    }
    
    prior_var <- Sigma_priorX[c(1, which_in + 1), c(1, which_in + 1)]
    prior_var_inv <- chol2inv(chol(prior_var))
    post_var <- prior_var_inv + t(C) %*% C / curr_variance
    post_var <- chol2inv(chol(post_var))
    
    prior_mean <- mu_priorX[c(1, which_in + 1)]
    post_mean <- prior_var_inv %*% prior_mean
    post_mean <- post_mean + t(C) %*% matrix(D$X, ncol = 1) / curr_variance
    post_mean <- post_var %*% post_mean
    
    gen_coef <- mvnfast::rmvn(1, mu = post_mean, sigma = post_var)
    r[1, ee, c(1, which_in + 1)] <- gen_coef
    
    
    # For the outcome model, only update coefficients of covariates.
    if (num_conf > 0) {
      
      curr_variance <- current_vars[2, ee]
      current_alphaY <- current_alphas[2, ee, ]
      which_in <- which(current_alphaY == 1)
      cov_cols_in <- cov_cols[which_in]
      C <- as.matrix(D[, cov_cols_in, drop = FALSE])
      
      # Prior of the coefficients for the covariates.
      prior_var <- Sigma_priorY[which_in + 2, which_in + 2]
      prior_var_inv <- chol2inv(chol(prior_var))
      post_var <- prior_var_inv + t(C) %*% C / curr_variance
      post_var <- chol2inv(chol(post_var))
      
      resid <- D$Y - cbind(1, D$X - exact_cuts[ee]) %*% current_coefs[2, ee, 1 : 2]
      prior_mean <- mu_priorY[which_in + 2]
      post_mean <- prior_var_inv %*% prior_mean
      post_mean <- post_mean + t(C) %*% matrix(resid, ncol = 1) / curr_variance
      post_mean <- post_var %*% post_mean
      
      gen_coef <- mvnfast::rmvn(1, mu = post_mean, sigma = post_var)
      r[2, ee, which_in + 1] <- gen_coef
    }
  }
  
  return(r)
}