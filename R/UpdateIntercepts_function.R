UpdateIntercepts <- function(dta, cov_cols, current_cutoffs, current_coefs,
                             current_vars, Sigma_priorY, mu_priorY) {
  
  minX <- min(dta$X)
  maxX <- max(dta$X)
  num_conf <- length(cov_cols)
  K <- length(current_cutoffs)
  exact_cuts <- c(minX, current_cutoffs, maxX)
  
  res <- rep(NA, K + 1)
  n_k <- table(dta$E)

  prior_var <- Sigma_priorY[1, 1]
  prior_mean <- mu_priorY[1]
  
  # Updating the intercept of the first experiment, and the rest follow.
  
  # Posterior variance.
  post_var <- 1 / prior_var + sum(n_k / current_vars[2, ])
  post_var <- 1 / post_var
  
  # Posterior mean.
  post_mean <- prior_mean / prior_var
  
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, E == ee)
    des_mat <- cbind(1, D$X - exact_cuts[ee], as.matrix(D[, cov_cols]))
    use_coefs <- c(0, current_coefs[2, ee, - 1])
    # The intercept as calculated below should be the same as the sum of beta
    # times the interval length.
    use_coefs[1] <- current_coefs[2, ee, 1] - current_coefs[2, 1, 1]
    
    resid <- D$Y - des_mat %*% matrix(use_coefs, ncol = 1)
    post_mean <- post_mean + sum(resid) / current_vars[2, ee]
  }
  res[1] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
  
  for (ee in 1 : K) {
    interval <- exact_cuts[ee + 1] - exact_cuts[ee]
    res[ee + 1] <- res[ee] + current_coefs[2, ee, 2] * interval
  }
  
  return(res)
}
