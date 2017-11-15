#' Local Exposure-Response Confounding Adjustment.
#' 
#' Function that calculates the mean exposure response curve allowing for
#' differential confounding at different exposure levels.
#' 
#' @export
LERCA <- function(dta, chains, Nsims, K, cov_cols, omega = 5000,
                  mu_priorX = NULL, Sigma_priorX = NULL,
                  mu_priorY = NULL, Sigma_priorY = NULL,
                  alpha_priorX = 0.001, beta_priorX = 0.001,
                  alpha_priorY = 0.001, beta_priorY = 0.001,
                  starting_cutoffs = NULL,
                  prop_distribution = c('Uniform', 'Normal'),
                  normal_percent = 1, plot_every = 0,
                  comb_probs = c(0.01, 0.5, 0.99),
                  split_probs = c(0.2, 0.95),
                  s_upd_probs = c(99 / 100, 1 / 100)) {
  
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
  
  acc <- array(0, dim = c(3, 2, chains))
  dimnames(acc) <- list(kind = c('separate', 'jumpOver', 'jumpWithin'),
                        num = c('attempt', 'success'), chain = 1 : chains)  
  
  # -------- STEP 3. MCMC. ---------- #
  
  for (cc in 1 : chains) {
    for (ii in 2 : Nsims) {
      
      if (ii %% 100 == 0) {
        print(paste('Iteration', ii))
      }
      
      # ----- Update experiment configuration and alphas.
      current_cutoffs <- cutoffs[cc, ii - 1, ]
      current_coefs <- coefs[, cc, ii - 1, , ]
      current_vars <- variances[, cc, ii - 1, ]
      current_alphas <- alphas[, cc, ii - 1, , ]
      
      wh_s_upd <- sample(c(1, 2, 3), 1, prob = s_upd_probs)
      acc[wh_s_upd, 1, cc] <- acc[wh_s_upd, 1, cc] + 1
      
      if (wh_s_upd == 1) {
        
        exp_upd <- UpdateExperiments(dta = dta, cov_cols = cov_cols,
                                     current_cutoffs = current_cutoffs,
                                     current_coefs = current_coefs,
                                     current_vars = current_vars,
                                     prop_distribution = prop_distribution,
                                     normal_percent = normal_percent)
        cutoffs[cc, ii, ] <- exp_upd$cutoffs
        acc[1, 2, cc] <- acc[1, 2, cc] + exp_upd$acc
        
        # ---- Update alphas.
        current_cutoffs <- cutoffs[cc, ii, ]
        current_alphaY <- alphas[2, cc, ii - 1, , ]
        
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
        
      } else if (wh_s_upd == 2) {
        
        jump_upd <- JumpOver(dta = dta, current_cutoffs = current_cutoffs,
                             current_alphas, approximate = TRUE,
                             cov_cols = cov_cols, omega = omega,
                             comb_probs = comb_probs, split_probs = split_probs,
                             min_exper_sample = min_exper_sample)
        
        acc[2, 2, cc] <- acc[2, 2, cc] + jump_upd$acc
        cutoffs[cc, ii, ] <- jump_upd$new_cutoffs
        alphas[, cc, ii, , ] <- jump_upd$new_alphas
        
      } else if (wh_s_upd == 3) {  # Within experiment simultaneous update.
        
        jump_upd <- JumpWithin(dta = dta, current_cutoffs = current_cutoffs,
                               current_alphas = current_alphas,
                               cov_cols = cov_cols, approximate = TRUE,
                               omega = omega, alpha_probs = alpha_probs,
                               min_exper_sample = min_exper_sample)
        
        acc[3, 1, cc] <- acc[3, 1, cc] + jump_upd$acc
        cutoffs[cc, ii, ] <- jump_upd$new_cutoffs
        alphas[, cc, ii, , ] <- jump_upd$new_alphas
        
      }
      
      # ----- Updating the coefficients.
      current_cutoffs <- cutoffs[cc, ii, ]
      current_alphas <- alphas[, cc, ii, , ]
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
      
      if (plot_every > 0) {
        if (ii %% plot_every == 0) {
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
              acc_percent = acc[, 2, ] / acc[, 1, ]))
}
