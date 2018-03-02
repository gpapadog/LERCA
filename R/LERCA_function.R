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
                  starting_alphas = NULL,
                  prop_distribution = c('Uniform', 'Normal'),
                  normal_percent = 1, plot_every = 0,
                  comb_probs = c(0.01, 0.5, 0.99),
                  split_probs = c(0.2, 0.95),
                  s_upd_probs = c(94 / 100, 1 / 100, 5 / 100),
                  alpha_probs = c(0.01, 0.5, 0.99),
                  min_exper_sample = 20) {
  
  progress <- floor(seq(2, Nsims, length.out = 11))[- 1]
  dta <- as.data.frame(dta)
  prop_distribution <- match.arg(prop_distribution)
  
  num_exper <- K + 1
  num_conf <- ifelse(is.null(cov_cols), 0, length(cov_cols))
  minX <- min(dta$X)
  maxX <- max(dta$X)
  
  if (num_conf == 0) {
    s_upd_probs <- c(1, 0, 0)
    warning('Update for experiment configuration without alphas.')
  } else {
    warning('JumpWithin can be improved to propose slopes.')
    warning('JumpOver is not performed with ensuring continuous ER.')
    if (any(abs(colMeans(dta[, cov_cols])) > 1e-10)) {
      stop('Covariates need to be centered.')
    }
  }
  
  
  if (is.null(mu_priorX)) {
    mu_priorX <- rep(0, num_conf + 1)
  }
  if (is.null(mu_priorY)) {
    # For a continuous ER, the intercepts of experiments higher than 1 are not
    # random, and the corresponding prior on mu_priorY is ignored.
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
                       maxX = maxX, starting_cutoffs = starting_cutoffs,
                       starting_alphas = starting_alphas)
  cutoffs <- arrays$cutoffs
  coefs <- arrays$coefs
  variances <- arrays$variances
  alphas <- arrays$alphas
  
  acc <- array(0, dim = c(3, 2, chains))
  dimnames(acc) <- list(kind = c('separate', 'jumpOver', 'jumpWithin'),
                        num = c('attempt', 'success'), chain = 1 : chains)
  
  # -------- STEP 3. MCMC. ---------- #
  
  for (cc in 1 : chains) {
    for (ii in 2 : Nsims) {
      
      if (ii %in% progress) {
        print(paste0('Chain ', cc, ' - ', 10 * which(progress == ii), '% done.'))
      }
      
      current_coefs <- coefs[, cc, ii - 1, , ]
      current_vars <- variances[, cc, ii - 1, ]
      current_alphas <- alphas[, cc, ii - 1, , ]  # NULL if no covariates.
      current_cutoffs <- cutoffs[cc, ii - 1, ]

      # ----- MCMC 1 : Update experiment configuration and alphas. ------ #
      
      wh_s_upd <- sample(c(1, 2, 3), 1, prob = s_upd_probs)
      acc[wh_s_upd, 1, cc] <- acc[wh_s_upd, 1, cc] + 1
      
      # Updating the experiment configuration and alphas separately.
      if (wh_s_upd == 1) {
        
        exp_upd <- UpdateExperiments(dta = dta, cov_cols = cov_cols,
                                     current_cutoffs = current_cutoffs,
                                     current_coefs = current_coefs,
                                     current_vars = current_vars,
                                     min_exper_sample = min_exper_sample,
                                     prop_distribution = prop_distribution,
                                     normal_percent = normal_percent,
                                     mu_priorY = mu_priorY,
                                     Sigma_priorY = Sigma_priorY)
        cutoffs[cc, ii, ] <- exp_upd$cutoffs
        current_cutoffs <- cutoffs[cc, ii, ]
        coefs[, cc, ii, , ] <- exp_upd$coefs  # Current or proposed.
        current_coefs <- coefs[, cc, ii, , ]
        acc[1, 2, cc] <- acc[1, 2, cc] + exp_upd$acc
        
        if (num_conf > 0) {
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
        }
        
      # Simulatinuous update: Jumping over the current experiment.
      } else if (wh_s_upd == 2) {
        
        jump_upd <- JumpOver(dta = dta, current_cutoffs = current_cutoffs,
                             current_alphas, approximate = TRUE,
                             cov_cols = cov_cols, omega = omega,
                             comb_probs = comb_probs, split_probs = split_probs,
                             min_exper_sample = min_exper_sample)
        
        acc[2, 2, cc] <- acc[2, 2, cc] + jump_upd$acc
        cutoffs[cc, ii, ] <- jump_upd$new_cutoffs
        alphas[, cc, ii, , ] <- jump_upd$new_alphas
       
      # Simulatinuous update: Jumping within the current experiment.
      } else if (wh_s_upd == 3) {
        
        jump_upd <- JumpWithin(dta = dta, current_cutoffs = current_cutoffs,
                               current_alphas = current_alphas,
                               current_coefs = current_coefs,
                               cov_cols = cov_cols, approximate = TRUE,
                               omega = omega, alpha_probs = alpha_probs,
                               min_exper_sample = min_exper_sample)
        
        acc[3, 1, cc] <- acc[3, 1, cc] + jump_upd$acc
        cutoffs[cc, ii, ] <- jump_upd$new_cutoffs
        alphas[, cc, ii, , ] <- jump_upd$new_alphas
        coefs[, cc, ii, , ] <- jump_upd$new_coefs
      }
      
      
      # Defining experiment based on current configuration.
      cuts <- c(minX - 0.001, current_cutoffs, maxX + 0.001)
      dta$E <- sapply(dta$X, function(x) sum(x >= cuts))
      
      
      # ------ Updating the covariates' coefficients.
      
      # If no covariates exist, we only update the intercepts of the exposure.
      # When covariates are centered, the outcome coefficients can be updated
      # separately from the remaining parameters.
      
      coefs[, cc, ii, , ] <- coefs[, cc, ii - 1, , ]
      coef_upd <- UpdateCovCoef(dta = dta, cov_cols = cov_cols,
                                current_cutoffs = current_cutoffs,
                                current_coefs = current_coefs,
                                current_alphas = current_alphas,
                                current_vars = current_vars,
                                Sigma_priorX = Sigma_priorX,
                                mu_priorX = mu_priorX,
                                Sigma_priorY = Sigma_priorY,
                                mu_priorY = mu_priorY)
      coefs[1, cc, ii, , - 2] <- coef_upd[1, , ]
      coefs[2, cc, ii, , - c(1, 2)] <- coef_upd[2, , - 1]
      current_coefs <- coefs[, cc, ii, , ]
      
      
      # ------ Updating the variances.
      
      var_upd <- UpdateVariances(dta = dta, current_cutoffs = current_cutoffs,
                                 current_coefs = current_coefs,
                                 cov_cols = cov_cols,
                                 alpha_priorX = alpha_priorX,
                                 beta_priorX = beta_priorX,
                                 alpha_priorY = alpha_priorY,
                                 beta_priorY = beta_priorY)
      variances[, cc, ii, ] <- var_upd
      current_vars <- variances[, cc, ii, ]
      
      
      # ----- Updating the intercept and slope of outcome models.
      
      int_slope_upd <- UpdateIntSlope(dta = dta, cov_cols = cov_cols,
                                      current_cutoffs = current_cutoffs,
                                      current_coefs = current_coefs,
                                      current_vars = current_vars,
                                      Sigma_priorY = Sigma_priorY,
                                      mu_priorY = mu_priorY)
      coefs[2, cc, ii, , c(1, 2)] <- int_slope_upd
      
      
      
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


# predict_at <- seq(min(dta$X), max(dta$X), length.out = 100)
# y_predict <- rep(NA, length(predict_at))
# for (xx in 1 : 100) {
#   cuts <- c(min(dta$X), cutoffs[cc, ii, ])
#   wh_exper <- sum(cutoffs[cc, ii, ] < predict_at[xx]) + 1
#   pred_vec <- c(1, predict_at[xx] - cuts[wh_exper])
#   y_predict[xx] <- sum(coefs[2, cc, ii, wh_exper, 1 : 2] * pred_vec)
# }
# plot(predict_at, y_predict, type = 'l', main = 'ER - one iteration')
# abline(v = cuts[- 1])
