#' Intercepts and coefficients of exposure
#' 
#' Updating the outcome model intercepts and coefficients for the exposure.
#' 
#' @param dta A data set including a column of the exposure of interest as X,
#' the outcome of interest as Y, and all potential confounders as C1, C2, ...
#' @param cov_cols The indices of the columns in dta corresponding to the
#' potential confounders.
#' @param current_cutoffs Numeric of length K. The current values for the
#' points in the experiment configuraiton.
#' @param current_coefs The current coefficients of the MCMC. Three dimensional
#' array with dimensions corresponding to exposure/outcome model, experiment,
#' and coefficients (intercept, slope, covariates).
#' @param current_vars Matrix. Rows correspond to exposure/outcome model, and
#' columns to experiments. Entries are the current variances.
#' @param mu_priorY The mean of the normal prior on the coefficients of the
#' outcome model. Numeric vector with entries corresponding to intercept, slope
#' of exposure, and potential covariates.  If left NULL, it is set to 0 for all
#' parameters.
#' @param Sigma_priorY The covariance matrix of the normal prior on the
#' regression coefficients of the outcome model. If left NULL, it is set to
#' diagonal with entries 100 ^ 2 for all parameters.
#' 
#' @return A matrix with dimensions that correspond to the exposure/outcome
#' model, and the intercept/slope coefficient.
#'
UpdateIntSlope <- function(dta, cov_cols, current_cutoffs, current_coefs,
                           current_vars, mu_priorY, Sigma_priorY) {
  
  num_conf <- ifelse(is.null(cov_cols), 0, length(cov_cols))
  K <- length(current_cutoffs)
  exact_cuts <- c(min(dta$X), current_cutoffs, max(dta$X))
  n_k <- table(dta$E)  # Number of observations per experiment.
  
  # Intercepts and slopes set at current values.
  r <- current_coefs[2, , 1 : 2]
  
  
  # -------- Updating the intercept of the first experiment. -------- #
  
  prior_var <- Sigma_priorY[1, 1]
  post_var <- 1 / prior_var + sum(n_k / current_vars[2, ])
  post_var <- 1 / post_var
  
  prior_mean <- mu_priorY[1]
  post_mean <- prior_mean / prior_var
  
  for (ee in 1 : (K + 1)) {

    D <- subset(dta, E == ee)
    des_mat <- cbind(1, D$X - exact_cuts[ee])
    des_coef <- c(r[ee, 1] - r[1, 1], r[ee, 2])
    
    if (num_conf > 0) {
      des_mat <- cbind(des_mat, as.matrix(D[, cov_cols]))
      des_coef <- c(des_coef, current_coefs[2, ee, - c(1, 2)])
    }
    
    des_coef <- matrix(des_coef, ncol = 1)
    resid <- D$Y - des_mat %*% des_coef
    post_mean <- post_mean + sum(resid) / current_vars[2, ee]
  }
  post_mean <- post_var * post_mean
  r[1, 1] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
  
  
  # -------- Updating the intercept of remaing experiments. -------- #
  
  # Here, we use that the covariates are centered.
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
    
    post_mean <- prior_mean / prior_var
    
    des_mat <- matrix(1, nrow = nrow(D), ncol = 1)
    des_coef <- r[ee, 1]
    if (num_conf > 0) {
      des_mat <- cbind(des_mat, as.matrix(D[, cov_cols]))
      des_coef <- c(des_coef, current_coefs[2, ee, - c(1, 2)])
    }
    des_coef <- matrix(des_coef, ncol = 1)
    resid <- D$Y - des_mat %*% des_coef
    post_mean <- post_mean + sum(X_s * resid) / curr_var
    
    if (ee != K + 1) {
      for (ll in (ee + 1) : (K + 1)) {
        
        D_ll <- subset(dta, E == ll)
        des_mat_ll <- cbind(1, D_ll$X - exact_cuts[ll])
        des_coef_ll <- c(r[ll, 1] - r[ee, 2] * interval, r[ll, 2])
        if (num_conf > 0) {
          des_mat_ll <- cbind(des_mat_ll, as.matrix(D_ll[, cov_cols]))
          des_coef_ll <- c(des_coef_ll, current_coefs[2, ll, - c(1, 2)])
        }
        des_coef_ll <- matrix(des_coef_ll, ncol = 1)
        resid <- D_ll$Y - des_mat_ll %*% des_coef_ll
        
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

