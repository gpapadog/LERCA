# current_cutoffs <- cutoffs[cc, ii, ]
# current_alphaY <- alphas[1, cc, ii - 1, , ]
# current_coefs <- coefs[, cc, ii, , ]
# current_vars <- variances[, cc, ii, ]

UpdateAlphas <- function(dta, cov_cols, current_cutoffs, current_alphaY,
                         current_coefs, current_vars, Sigma_priorX, mu_priorX,
                         Sigma_priorY, mu_priorY,
                         omega) {
  
  minX <- min(dta$X) - 0.00001
  maxX <- max(dta$X) + 0.00001
  num_conf <- length(cov_cols)
  
  r <- array(NA, dim = c(2, dim(current_alphaY)))
  
  # ---- Update alphas.
  for (ee in 1 : (K + 1)) {
    
    cuts <- c(minX, current_cutoffs, maxX)
    D <- subset(dta, X > cuts[ee] & X <= cuts[ee + 1])

    for (jj in 1 : num_conf) {
      
      # Exposure model.
      corresp_alphaY <- current_alphaY[ee, jj]
      
      alpha0 <- ifelse(corresp_alphaY == 0, omega / (omega + 1), 1 / 2)
      
      curr_cov <- D[, cov_cols[jj]]
      other_cov <- cbind(1, as.matrix(D[, cov_cols[- jj]]))
      resid <- D$X - other_cov %*% current_coefs[1, ee, - c(2, jj + 2)]
      
      prior_var <- Sigma_priorX[jj + 1, jj + 1]
      prior_sd <- sqrt(prior_var)
      post_var <- sum(curr_cov ^ 2) / current_vars[1, ee]
      post_var <- 1 / (post_var + 1 / prior_var)
      
      post_mean <- sum(curr_cov * resid) / current_vars[1, ee]
      post_mean <- post_mean + mu_priorX[jj + 1] / prior_var
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
      other_cov <- cbind(1, D$X, as.matrix(D[, cov_cols[- jj]]))
      resid <- D$Y - other_cov %*% current_coefs[2, ee, - (jj + 2)]
      
      prior_var <- Sigma_priorY[jj + 2, jj + 2]
      prior_sd <- sqrt(prior_var)
      post_var <- sum(curr_cov ^ 2) / current_vars[2, ee]
      post_var <- 1 / (post_var + 1 / prior_var)
      
      post_mean <- sum(curr_cov * resid) / current_vars[2, ee]
      post_mean <- post_mean + mu_priorY[jj + 1] / prior_var
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