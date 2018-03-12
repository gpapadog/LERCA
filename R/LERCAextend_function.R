#' Extending the LERCA fit.
#' 
#' Acquiring more posterior samples from a LERCA fit.
#' 
#' @param dta The data set including the exposure as X, the outcome as Y, and
#' all potential confounders as C1, C2, ...
#' @param extend_Nsims The number of additional posterior sampels.
#' @param lerca The previous LERCA fit.
#' @param plot_every At how many iterations we want to plot the posterior
#' samples of the experiment configuration. Defaults to 0.
#' 
#' @return Returns a list with the same elements as LERCA including posterior
#' samples of the initial and the extended fit.
#' 
LERCAextend <- function(dta, extend_Nsims, lerca, plot_every = 0) {
  
  chains <- dim(lerca$cutoffs)[1]
  K <- dim(lerca$cutoffs)[3]
  prev_Nsims <- dim(lerca$cutoffs)[2]
  
  starting_cutoffs <- lerca$cutoffs[, prev_Nsims, , drop = FALSE]
  starting_cutoffs <- abind::adrop(starting_cutoffs, drop = 2)
  starting_alphas <- lerca$alphas[, , prev_Nsims, , , drop = FALSE]
  starting_alphas <- abind::adrop(starting_alphas, drop = 3)
  starting_coefs <- lerca$coefs[, , prev_Nsims, , , drop = FALSE]
  starting_coefs <- abind::adrop(starting_coefs, drop = 3)
  
  starting_vars <- lerca$variances[, , prev_Nsims, , drop = FALSE]
  
  cov_cols <- lerca$lerca_specs$cov_cols
  omega <- lerca$lerca_specs$omega
  mu_priorX <- lerca$lerca_specs$mu_priorX
  Sigma_priorX <- lerca$lerca_specs$Sigma_priorX
  mu_priorY <- lerca$lerca_specs$mu_priorY
  Sigma_priorY <- lerca$lerca_specs$Sigma_priorY
  alpha_priorX <- lerca$lerca_specs$alpha_priorX
  beta_priorX <- lerca$lerca_specs$beta_priorX
  alpha_priorY <- lerca$lerca_specs$alpha_priorY
  beta_priorY <- lerca$lerca_specs$beta_priorY
  approx_likelihood <- lerca$lerca_specs$approx_likelihood
  prop_distribution <- lerca$lerca_specs$prop_distribution
  normal_percent <- lerca$lerca_specs$normal_percent
  comb_probs <- lerca$lerca_specs$comb_probs
  split_probs <- lerca$lerca_specs$split_probs
  s_upd_probs <- lerca$lerca_specs$s_upd_probs
  alpha_probs <- lerca$lerca_specs$alpha_probs
  min_exper_sample <- lerca$lerca_specs$min_exper_sample
  

  lerca_ext <- LERCA(dta = dta, chains = chains, Nsims = extend_Nsims, K = K,
                     cov_cols = cov_cols, omega = omega, mu_priorX = mu_priorX,
                     mu_priorY = mu_priorY, Sigma_priorX = Sigma_priorX,
                     Sigma_priorY = Sigma_priorY, alpha_priorX = alpha_priorX,
                     beta_priorX = beta_priorX, alpha_priorY = alpha_priorY,
                     beta_priorY = beta_priorY,
                     starting_cutoffs = starting_cutoffs,
                     starting_alphas = starting_alphas,
                     starting_coefs = starting_coefs,
                     starting_vars = starting_vars,
                     approx_likelihood = approx_likelihood,
                     prop_distribution = prop_distribution,
                     normal_percent = normal_percent, plot_every = plot_every,
                     comb_probs = comb_probs, split_probs = split_probs,
                     s_upd_probs = s_upd_probs, alpha_probs = alpha_probs,
                     min_exper_sample = min_exper_sample)

  
  cutoffs <- abind(lerca$cutoffs, lerca_ext$cutoffs, along = 2)
  alphas <- abind(lerca$alphas, lerca_ext$alphas, along = 3)
  coefs <- abind(lerca$coefs, lerca_ext$coefs, along = 3)
  variances <- abind(lerca$variances, lerca_ext$variances, along = 3)
  acc <- lerca$acc + lerca_ext$acc
  acc_percent <- acc[, 2, ] / acc[, 1, ]

  return(list(cutoffs = cutoffs, alphas = alphas, coefs = coefs,
              variances = variances, acc = acc, acc_percent = acc_percent,
              lerca_specs = lerca$lerca_specs))
  
}