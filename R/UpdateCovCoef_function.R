#' Update within experiment coefficients
#' 
#' Updating the coefficients that use the within experiment likelihood:
#' intercept and coefficients of the exposure model, and the coefficients of
#' the covariates in the outcome model.
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
#' @param current_alphas Array of dimensions that correspond to the exposure or
#' outcome model, the experiment, and potential confounding. Entries are 0/1
#' corresponding to exclusion/inclusion of the covaraite in the corresponding
#' model of the experiment.
#' @param current_vars Matrix. Rows correspond to exposure/outcome model, and
#' columns to experiments. Entries are the current variances.
#' @param mu_priorX The mean of the normal prior on the coefficients of the
#' exposure model. Numeric vector of length equal to the number of potential
#' confounders + 1 with the first entry corresponding to the intercept.
#' @param Sigma_priorX Covariance matrix of the normal prior on the regression
#' coefficients of the exposure model.
#' @param mu_priorY The mean of the normal prior on the coefficients of the
#' outcome model. Numeric vector with entries corresponding to intercept, slope
#' of exposure, and potential covariates.
#' @param Sigma_priorY The covariance matrix of the normal prior on the
#' regression coefficients of the outcome model.
#' 
#' @return Array of dimensions that correspond to the exposure/outcome model,
#' experiment, and coefficients (intercept, covariates). The intercepts of the
#' outcome model are NA, since they are not updated with this function.
#' 
UpdateCovCoef <- function(dta, cov_cols, current_cutoffs, current_coefs,
                          current_alphas, current_vars, mu_priorX,
                          Sigma_priorX, mu_priorY, Sigma_priorY) {
  
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
      if (length(which_in) > 0) {
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
  }
  
  return(r)
}