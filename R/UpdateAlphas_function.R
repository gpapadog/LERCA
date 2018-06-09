#' Update the inclusion indicators from their full conditional.
#' 
#' Use Gibbs sampling to update the inclusion indicator of all covariates in
#' the exposure and outcome models.
#' 
#' @param dta A data set including a column of the exposure of interest as X,
#' the outcome of interest as Y, and all potential confounders as C1, C2, ...
#' @param cov_cols The indices of the columns in dta corresponding to the
#' potential confounders.
#' @param current_cutoffs Numeric of length K. The current values for the
#' points in the experiment configuraiton.
#' @param current_alphaY Matrix of 0/1 entries. Rows correspond to experiments,
#' and columns to covariate. Entries 0/1 represent exclusion/inclusion of the
#' covaraite in the corresponding outcome model.
#' @param current_coefs The current coefficients of the MCMC. Three dimensional
#' array with dimensions corresponding to exposure/outcome model, experiment,
#' and coefficients (intercept, slope, covariates).
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
#' @param omega The omega of the BAC prior.
#' 
#' @return Array of the inclusion indicators where the dimensions correspond to
#' the exposure/outcome model, experiment and covariates.
#' 
UpdateAlphas <- function(dta, cov_cols, current_cutoffs, current_alphaY,
                         current_coefs, current_vars, mu_priorX, Sigma_priorX,
                         mu_priorY, Sigma_priorY, omega) {
  
  K <- length(current_cutoffs)
  minX <- min(dta$X)
  maxX <- max(dta$X)
  num_conf <- length(cov_cols)
  cuts <- c(minX - 0.001, current_cutoffs, maxX + 0.001)
  exact_cuts <- c(minX, current_cutoffs, maxX)
  
  r <- array(NA, dim = c(2, dim(current_alphaY)))
  
  # ---- Update alphas.
  for (ee in 1 : (K + 1)) {
    
    D <- subset(dta, X > cuts[ee] & X <= cuts[ee + 1])

    for (jj in 1 : num_conf) {
      
      # Exposure model.
      corresp_alphaY <- current_alphaY[ee, jj]
      
      alpha0 <- ifelse(corresp_alphaY == 0, omega / (omega + 1), 1 / 2)
      
      curr_cov <- D[, cov_cols[jj]]
      other_cov <- cbind(Int = 1, as.matrix(D[, cov_cols[- jj]]))
      resid <- D$X - other_cov %*% current_coefs[1, ee, - c(2, jj + 2)]
      
      prior_var <- Sigma_priorX[jj + 1, jj + 1]
      prior_sd <- sqrt(prior_var)
      post_var <- 1 / prior_var + sum(curr_cov ^ 2) / current_vars[1, ee]
      post_var <- 1 / post_var
      
      post_mean <- mu_priorX[jj + 1] / prior_var
      post_mean <- post_mean + sum(curr_cov * resid) / current_vars[1, ee]
      post_mean <- post_var * post_mean
      
      alpha1 <- corresp_alphaY * log(1 / 2) +
        (1 - corresp_alphaY) * log(1 / (omega + 1)) +
        dnorm(x = 0, mean = mu_priorX[jj + 1], sd = prior_sd, log = TRUE) -
        dnorm(0, mean = post_mean, sd = sqrt(post_var), log = TRUE)
      alpha1 <- exp(alpha1)
      
      probs <- c(alpha0, alpha1)
      if (is.infinite(probs[2]) & is.finite(probs[1])) {
        probs <- c(0, 1)
      }
      r[1, ee, jj] <- sample(c(0, 1), size = 1, prob = probs)
      
      # Outcome model.
      corresp_alphaX <- r[1, ee, jj]
      
      alpha0 <- ifelse(corresp_alphaX == 0, 1 / 2, 1 / (omega + 1))
      
      curr_cov <- D[, cov_cols[jj]]
      other_cov <- cbind(Int = 1, X_s = D$X - exact_cuts[ee])
      other_cov <- cbind(other_cov, as.matrix(D[, cov_cols[- jj]]))
      resid <- D$Y - other_cov %*% current_coefs[2, ee, - (jj + 2)]
      
      prior_var <- Sigma_priorY[jj + 2, jj + 2]
      prior_sd <- sqrt(prior_var)
      post_var <- 1 / prior_var + sum(curr_cov ^ 2) / current_vars[2, ee]
      post_var <- 1 / post_var
      
      post_mean <- mu_priorY[jj + 2] / prior_var
      post_mean <- post_mean + sum(curr_cov * resid) / current_vars[2, ee]
      post_mean <- post_var * post_mean
      
      alpha1 <- corresp_alphaX * log(omega / (omega + 1)) +
        (1 - corresp_alphaX) * log(1 / 2) +
        dnorm(x = 0, mean = mu_priorY[jj + 2], sd = prior_sd, log = TRUE) -
        dnorm(0, mean = post_mean, sd = sqrt(post_var), log = TRUE)
      alpha1 <- exp(alpha1)
      
      probs <- c(alpha0, alpha1)
      if (is.infinite(probs[2]) & is.finite(probs[1])) {
        probs <- c(0, 1)
      }
      r[2, ee, jj] <- sample(c(0, 1), size = 1, prob = probs)
    }
  }
  return(r)
}