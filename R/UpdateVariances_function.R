#' Exposure and outcome model residual variances
#' 
#' Updating the residual variance of the exposure and outcome models within
#' each experiment.
#' 
#' @param dta A data set including a column of the exposure of interest as X,
#' the outcome of interest as Y, and all potential confounders as C1, C2, ...
#' @param current_cutoffs Numeric of length K. The current values for the
#' points in the experiment configuraiton.
#' @param current_coefs The current coefficients of the MCMC. Three dimensional
#' array with dimensions corresponding to exposure/outcome model, experiment,
#' and coefficients (intercept, slope, covariates).
#' @param cov_cols The indices of the columns in dta corresponding to the
#' potential confounders.
#' @param alpha_priorX The shape parameter of the inverse gamma prior on the
#' residual variance of the exposure model.
#' @param beta_priorX The rate parameter of the inverse gamma prior on the
#' residual variance of the exposure model.
#' @param alpha_priorY The shape parameter of the inverse gamma prior on the
#' residual variance of the outcome model.
#' @param beta_priorY The rate parameter of the inverse gamma prior on the
#' residual variance of the outcome model.
#' 
#' @return Matrix with rows corresponding to the exposure/outcome model and
#' columns corresponding to the experiment.
#' 
UpdateVariances <- function(dta, current_cutoffs, current_coefs, cov_cols,
                            alpha_priorX, beta_priorX, alpha_priorY,
                            beta_priorY) {
  
  minX <- min(dta$X) - 0.001
  maxX <- max(dta$X) + 0.001
  K <- length(current_cutoffs)
  cuts <- c(minX, current_cutoffs, maxX)
  
  r <- array(0, dim = c(2, K + 1))
  
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, E == ee)
    n_k <- nrow(D)  # Number of observations in current experiment.

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
  
