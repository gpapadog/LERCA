UpdateExpCoef <- function(dta, cov_cols, current_cutoffs, current_coefs,
                          current_vars, Sigma_priorY, mu_priorY) {
  
  minX <- min(dta$X)
  maxX <- max(dta$X)
  num_conf <- length(cov_cols)
  K <- length(current_cutoffs)
  
  # Coefficient of exposure at different experiments, set at current values.
  res <- current_coefs[2, , 2]
  
  exact_cuts <- c(minX, current_cutoffs, maxX)
  
  prior_var <- Sigma_priorY[2, 2]
  prior_mean <- mu_priorY[2]
  
  # Defining experiment based on current configuration.
  n_k <- table(dta$E)
  
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, E == ee)

    des_mat <- cbind(1, D[, cov_cols])
    resid <- D$Y - des_mat %*% current_coefs[2, ee, - 2, drop = FALSE]
    X_s <- D$X - exact_cuts[ee]
    
    
    # Posterior variance.
    
    post_var <- 1 / prior_var + sum(X_s ^ 2) / current_vars[2, ee]
    interval <- exact_cuts[ee + 1] - exact_cuts[ee]
    if (ee != K + 1) {
      for (ll in (ee + 1) : (K + 1)) {
        post_var <- post_var + n_k[ll] / current_vars[2, ll] * (interval ^ 2)
      }
    }
    post_var <- 1 / post_var
    
    
    # Posterior mean.
    
    post_mean <- prior_mean / prior_var
    post_mean <- post_mean + sum(X_s * resid) / current_vars[2, ee]
    if (ee != K + 1) {
      for (ll in (ee + 1) : (K + 1)) {
        D_ll <- subset(dta, E == ll)
        des_mat_ll <- cbind(1, D_ll$X - exact_cuts[ll], D_ll[, cov_cols])
        
        # What are the coefficients of the linear model for the residuals.
        coef_ll <- c(0, res[ll], current_coefs[2, ll, - c(1, 2)])
        # Intercept.
        coef_ll[1] <- current_coefs[2, ee, 1]
        if (ee + 1 <= ll - 1) {
          for (rr in (ee + 1) : (ll - 1)) {
            coef_ll[1] <- coef_ll[1] + res[rr] * (exact_cuts[rr + 1] -
                                                    exact_cuts[rr])
          }
        }
        coef_ll <- matrix(coef_ll, ncol = 1)
        resid_ll <- D_ll$Y - des_mat_ll %*% coef_ll
        
        # Adding to the posterior mean.
        post_mean <- post_mean + sum(resid_ll) * interval / current_vars[2, ll]
      }
    }
    post_mean <- post_var * post_mean
    res[ee] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
  }
  return(res)
}

