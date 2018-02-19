UpdateIntSlope <- function(dta, cov_cols, current_cutoffs, current_coefs,
                           current_vars, Sigma_priorY, mu_priorY) {
  
  num_conf <- ifelse(is.null(cov_cols), 0, length(cov_cols))
  K <- length(current_cutoffs)
  exact_cuts <- c(min(dta$X), current_cutoffs, max(dta$X))
  n_k <- table(dta$E)  # Number of observations per experiment.
  
  if (num_conf > 0) {
    stop('Code has to change to incorporate covariates.')
  }
  
  # Intercepts and slopes set at current values.
  r <- current_coefs[2, , 1 : 2]
  
  
  # -------- Updating the intercept of the first experiment. -------- #
  
  prior_var <- Sigma_priorY[1, 1]
  prior_mean <- mu_priorY[1]
  
  post_var <- 1 / prior_var + sum(n_k / current_vars[2, ])
  post_var <- 1 / post_var
  
  post_mean <- prior_mean / prior_var
  for (ee in 1 : (K + 1)) {
    D <- subset(dta, E == ee)
    resid <- D$Y - (r[ee, 1] - r[1, 1]) - r[ee, 2] * (D$X - exact_cuts[ee])
    post_mean <- post_mean + sum(resid) / current_vars[2, ee]
  }
  post_mean <- post_var * post_mean
  r[1, 1] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
  
  
  # -------- Updating the intercept of remaing experiments. -------- #
  
  for (ee in 1 : K) {
    r[ee + 1, 1] <- r[ee, 1] + r[ee, 2] * (exact_cuts[ee + 1] - exact_cuts[ee])
  }
  
  
  # -------- Updating the slope of each experiment. -------- #
  
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, E == ee)

    prior_var <- Sigma_priorY[2, 2]
    curr_var <- current_vars[2, ee]
    prior_mean <- mu_priorY[2]
    
    interval <- exact_cuts[ee + 1] - exact_cuts[ee]
    X_s <- D$X - exact_cuts[ee]
    
    # ---- Updating the slope of current experiment.
    
    # Posterior variance.

    post_var <- 1 / prior_var + sum(X_s ^ 2) / curr_var
    if (ee != K + 1) {
      wh <- (ee + 1) : (K + 1)
      sum_vars <- sum(n_k[wh] / current_vars[2, wh])
      post_var <- post_var + (interval ^ 2) * sum_vars
    }
    post_var <- 1 / post_var
    
    
    # Posterior mean.
    
    resid <- D$Y - r[ee, 1]
    post_mean <- prior_mean / prior_var + sum(X_s * resid) / curr_var
    if (ee != K + 1) {
      for (ll in (ee + 1) : (K + 1)) {
        D_ll <- subset(dta, E == ll)
        resid_mean <- r[ll, 1] - r[ee, 2] * interval 
        resid_mean <- resid_mean + r[ll, 2] * (D_ll$X - exact_cuts[ll])
        resid <- D_ll$Y - resid_mean
        post_mean <- post_mean + sum(interval * resid) / current_vars[2, ll]
      }
    }
    post_mean <- post_var * post_mean
    
    r[ee, 2] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
    
    
    # ---- Updating all intercepts that include that slope ------
    
    if (ee != K + 1) {
      for (ll in ee : K) {
        r[ll + 1, 1] <- r[ll, 1] + r[ll, 2] * (exact_cuts[ll + 1] - exact_cuts[ll])
      }
    }
  }
  return(r)
}

