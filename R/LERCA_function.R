#' Local Exposure-Response Confounding Adjustment.
#' 
#' Function that calculates the mean exposure response curve allowing for
#' differential confounding at different exposure levels.
#' 
#' @param dta A data set including a column of the exposure of interest as X,
#' the outcome of interest as Y, and all potential confounders as C1, C2, ...
#' @param chains The number of MCMC chains.
#' @param Nsims The number of posterior samples per chain.
#' @param K The number of points in the experiment configuration.
#' @param cov_cols The indices of the columns in dta corresponding to the
#' potential confounders.
#' @param omega The parameter of the BAC prior on the inclusion indicators.
#' @param mu_priorX The mean of the normal prior on the coefficients of the
#' exposure model. Numeric vector of length equal to the number of potential
#' confounders + 1 with the first entry corresponding to the intercept. If left
#' NULL, it is set to 0 for all parameters.
#' @param Sigma_priorX Covariance matrix of the normal prior on the regression
#' coefficients of the exposure model. If left NULL, it is set to diagonal with
#' entries 100 ^ 2 for all parameters.
#' @param mu_priorY The mean of the normal prior on the coefficients of the
#' outcome model. Numeric vector with entries corresponding to intercept, slope
#' of exposure, and potential covariates.  If left NULL, it is set to 0 for all
#' parameters.
#' @param Sigma_priorY The covariance matrix of the normal prior on the
#' regression coefficients of the outcome model. If left NULL, it is set to
#' diagonal with entries 100 ^ 2 for all parameters.
#' @param alpha_priorX The shape parameter of the inverse gamma prior on the
#' residual variance of the exposure model.
#' @param beta_priorX The rate parameter of the inverse gamma prior on the
#' residual variance of the exposure model.
#' @param alpha_priorY The shape parameter of the inverse gamma prior on the
#' residual variance of the outcome model.
#' @param beta_priorY The rate parameter of the inverse gamma prior on the
#' residual variance of the outcome model.
#' @param starting cutoffs Matrix with rows corresponding to different chains.
#' Each row includes K ordered values of MCMC starting cutoffs. If left NULL,
#' random started values are used.
#' @param starting_alphas Array with dimensions corresponding to the model
#' (exposure / outcome), the experiment, and the potential confounders. Entries
#' 0/1 represent exclusion/inclusion of the covariate in the corresponding
#' model.
#' @param approx_likelihood Logical. If set to TRUE the likelihood of the data
#' in the jump over and jump within moves will be calculated based on the BIC
#' approximation. Defaults to TRUE. Option FALSE not supported for now.
#' @param prop_distribution Character vector. Options include 'Uniform' or
#' 'Normal' representing the type of distribution that will be used to propose
#' a move of the cutoffs in the separate update. Defaults to uniform.
#' @param normal_percent Numeric. Parameter controling the width of a normal
#' proposal for the experiment configuration. Used only when prop_distribution
#' is set to Normal. Smaller values represent smaller variance of the truncated
#' normal proposal distribution.
#' @param plot_every Integer. Plot the locations of the experiment
#' configuration every plot_every iteration. Defaults to 0 leading to no
#' plotting.
#' @param comb_probs When two experiments are combined, comb_probs represents
#' the probability of alpha = 1 when 0, 1, and 2 corresponding alphas are equal
#' to 1. Vector of length 3. Defaults to (0.01, 0.5, 0.99).
#' @param split_probs When one experiment is split, split_probs describes the
#' probability that the alpha of a new experiment is equal to 1, when the alpha
#' of the current experiment is 0, and when it is 1. Vector of length 2.
#' Defaults to (0.2, 0.95).
#' @param s_upd_probs Numeric of length three. The probability that each of the
#' separate, jump over, and jump within moves is proposed. Defaults to (0.8,
#' 0.1, 0.1).
#' @param alpha_probs The probability that a proposed alpha is equal to 1, when
#' 0, 1, and 2 alphas of the surrounding experiments are equal to 1. Vector of
#' length 3. Defaults to (0.01, 0.5, 0.99).
#' @param min_exper_sample The minimum number of observations within an
#' experiment. Defaults to 20.
#' @param jump_slope_tune The standard deviation of the proposal on the slopes
#' for the jump over move. Defaults to 0.05.
#' 
#' @export
LERCA <- function(dta, chains, Nsims, K, cov_cols, omega = 5000,
                  mu_priorX = NULL, Sigma_priorX = NULL,
                  mu_priorY = NULL, Sigma_priorY = NULL,
                  alpha_priorX = 0.001, beta_priorX = 0.001,
                  alpha_priorY = 0.001, beta_priorY = 0.001,
                  starting_cutoffs = NULL,
                  starting_alphas = NULL,
                  approx_likelihood = TRUE,
                  prop_distribution = c('Uniform', 'Normal'),
                  normal_percent = 1, plot_every = 0,
                  comb_probs = c(0.01, 0.5, 0.99),
                  split_probs = c(0.2, 0.95),
                  s_upd_probs = c(0.8, 0.1, 0.1),
                  alpha_probs = c(0.01, 0.5, 0.99),
                  min_exper_sample = 20, jump_slope_tune = 0.05) {
  
  
  prop_distribution <- match.arg(prop_distribution)
  if (!approx_likelihood) {
    stop('approx_likelihood can only be set to TRUE for now.')
  }
  
  progress <- floor(seq(2, Nsims, length.out = 11))[- 1]
  dta <- as.data.frame(dta)
  
  num_exper <- K + 1
  num_conf <- length(cov_cols)
  minX <- min(dta$X)
  maxX <- max(dta$X)
  
  if (num_conf == 0) {
    s_upd_probs <- c(1, 0, 0)
    warning('Update for experiment configuration without alphas.')
  } else {
    if (any(abs(colMeans(dta[, cov_cols])) > 1e-10)) {
      stop('Covariates need to be centered.')
    }
  }
  
  
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
  arrays <- MakeArrays(X = dta$X, chains = chains, Nsims = Nsims,
                       num_exper = num_exper, num_conf = num_conf,
                       omega = omega, minX = minX, maxX = maxX,
                       starting_cutoffs = starting_cutoffs,
                       starting_alphas = starting_alphas,
                       min_exper_sample = min_exper_sample)
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
                             current_coefs = current_coefs,
                             current_alphas = current_alphas,
                             approx_likelihood = approx_likelihood,
                             cov_cols = cov_cols, omega = omega,
                             Sigma_priorY = Sigma_priorY,
                             mu_priorY = mu_priorY, comb_probs = comb_probs,
                             split_probs = split_probs,
                             min_exper_sample = min_exper_sample,
                             jump_slope_tune = jump_slope_tune)
        
        acc[2, 2, cc] <- acc[2, 2, cc] + jump_upd$acc
        cutoffs[cc, ii, ] <- jump_upd$new_cutoffs
        alphas[, cc, ii, , ] <- jump_upd$new_alphas
        coefs[, cc, ii, , ] <- jump_upd$new_coefs
        
      # Simulatinuous update: Jumping within the current experiment.
      } else if (wh_s_upd == 3) {
        
        jump_upd <- JumpWithin(dta = dta, current_cutoffs = current_cutoffs,
                               current_alphas = current_alphas,
                               current_coefs = current_coefs,
                               cov_cols = cov_cols,
                               approx_likelihood = approx_likelihood,
                               omega = omega, mu_priorY = mu_priorY,
                               Sigma_priorY = Sigma_priorY,
                               alpha_probs = alpha_probs,
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
