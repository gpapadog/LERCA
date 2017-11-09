LERCAfullcond <- function(dta, chains, Nsims, K, cov_cols, omega = 5000,
                          mu_priorX = NULL, Sigma_priorX = NULL,
                          mu_priorY = NULL, Sigma_priorY = NULL,
                          alpha_priorX = 0.001, beta_priorX = 0.001,
                          alpha_priorY = 0.001, beta_priorY = 0.001,
                          starting_cutoffs = NULL,
                          prop_distribution = c('Normal', 'Uniform'),
                          normal_percent = 1, print_every = 1000) {
  
  prop_distribution <- match.arg(prop_distribution)
  
  num_exper <- K + 1
  num_conf <- length(cov_cols)
  minX <- min(dta$X)
  maxX <- max(dta$X)
  
  
  if (is.null(mu_priorX)) {
    mu_priorX <- rep(0, num_conf + 1)
  }
  if (is.null(mu_priorY)) {
    mu_priorY <- rep(0, num_conf + 2)
  }
  if (is.null(Sigma_priorX)) {
    Sigma_priorX <- diag(rep(100 ^ 2, num_conf + 1))
  }
  if (is.null(Sigma_priorY)) {
    Sigma_priorY <- diag(rep(100 ^ 2, num_conf + 2))
  }
  
  
  # -------- Where to save --------- #
  arrays <- MakeArrays(chains = chains, Nsims = Nsims, num_exper = num_exper,
                       num_conf = num_conf, omega = omega, minX = minX,
                       maxX = maxX, starting_cutoffs = starting_cutoffs)
  alphas <- arrays$alphas
  cutoffs <- arrays$cutoffs
  coefs <- arrays$coefs
  variances <- arrays$variances
  
  acc <- rep(0, chains)
  
  
  # -------- STEP 3. MCMC. ---------- #
  
  for (cc in 1 : chains) {
    for (ii in 2 : Nsims) {
      
      if (ii %% 100 == 0) {
        print(paste('Iteration', ii))
      }
      
      # ----- Update experiment configuration.
      current_cutoffs <- cutoffs[cc, ii - 1, ]
      current_coefs <- coefs[, cc, ii - 1, , ]
      current_vars <- variances[, cc, ii - 1, ]
      
      exp_upd <- UpdateExperiments(dta = dta, cov_cols = cov_cols,
                                   current_cutoffs = current_cutoffs,
                                   current_coefs = current_coefs,
                                   current_vars = current_vars,
                                   prop_distribution = prop_distribution,
                                   normal_percent = normal_percent)
      cutoffs[cc, ii, ] <- exp_upd$cutoffs
      acc[cc] <- acc[cc] + exp_upd$acc
      
      
      # ----- Updating the coefficients.
      current_cutoffs <- cutoffs[cc, ii, ]
      current_alphas <- alphas[, cc, ii - 1, , ]
      current_vars <- variances[, cc, ii - 1, ]
      
      coef_upd <- UpdateCoefficients(dta = dta, cov_cols = cov_cols,
                                     current_cutoffs = current_cutoffs,
                                     current_alphas = current_alphas,
                                     current_vars = current_vars,
                                     Sigma_priorX = Sigma_priorX,
                                     mu_priorX = mu_priorX,
                                     Sigma_priorY = Sigma_priorY,
                                     mu_priorY = mu_priorY)
      coefs[, cc, ii, , ] <- coef_upd

      # ------ Updating the variances.
      
      current_coefs <- coefs[, cc, ii, , ]
      var_upd <- UpdateVariances(dta = dta, current_cutoffs = current_cutoffs,
                                 current_coefs = current_coefs,
                                 alpha_priorX = alpha_priorX,
                                 beta_priorX = beta_priorX,
                                 alpha_priorY = alpha_priorY,
                                 beta_priorY = beta_priorY)
      variances[, cc, ii, ] <- var_upd

      # ---- Update alphas.
      current_cutoffs <- cutoffs[cc, ii, ]
      current_alphaY <- alphas[1, cc, ii - 1, , ]
      current_coefs <- coefs[, cc, ii, , ]
      current_vars <- variances[, cc, ii, ]
      
      alphas_upd <- UpdateAlphas(dta = dta, cov_cols = cov_cols,
                                 current_cutoffs = current_cutoffs,
                                 current_alphaY = current_alphaY,
                                 current_coefs = current_coefs,
                                 current_vars = current_vars,
                                 Sigma_priorX = Sigma_priorX,
                                 mu_priorX = mu_priorX,
                                 Sigma_priorY = Sigma_priorY,
                                 mu_priorY = mu_priorY, omega = omega)
      alphas[, cc, ii, , ] <- alphas_upd
      
      if (print_every > 0) {
        if (ii %% print_every == 0) {
          par(mfrow = c(2, ceiling(K / 2)), mar = rep(2, 4))
          
          for (kk in 1 : K) {
            plot(cutoffs[cc, 1 : ii, kk], type = 'l', col = cc)
          }
          
        }
      }
    }
  }
  return(list(cutoffs = cutoffs, alphas = alphas, coefs = coefs,
              variances = variances, acc = acc,
              acc_percent = acc / (Nsims - 1)))
}