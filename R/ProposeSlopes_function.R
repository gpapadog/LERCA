#' Proposed slope values of Jump Over move.
#' 
#' Propose values for the slopes of all experiments in the Jump Over move
#' ensuring continous ER and unchanged likelihood in experiments that are not
#' affected.
#' 
#' @param wh_cut Integer. Which one of the K points in the experiment
#' configuration was selected to be moved.
#' @param current_coefs The current coefficients in an array format, with
#' dimensions corresponding to the exposure/outcome model, the experiments, and
#' the coefficient (intercept, slope, covariates).
#' @param cuts The points in the experiment configuration includin the minimum
#' and maximum values.
#' @param sj1star The proposed point in the experiment configuration.
#' @param sj1 The current point that is proposed to be moved.
#' @param sj The point in the current experiment configuration before sj1.
#' @param sj2 The point in the current experiment configuration after sj1.
#' 
ProposeSlopes <- function(wh_cut, current_coefs, cuts, sj1star, sj1, sj, sj2) {
  
  proposed_coefs <- current_coefs
  
  # Finding the intercept at experiment wh_cut + 2.
  if (wh_cut < K) {
    delta_k2 <- current_coefs[2, wh_cut + 2, 1]  # Intercept of exp K + 2.
  } else if (wh_cut == K) {
    add_on <- current_coefs[2, wh_cut + 1, 2] * (cuts[K + 2] - cuts[K + 1])
    delta_k2 <- current_coefs[2, wh_cut + 1, 1] + add_on
  }
  
  if (sj1star > sj1) {
    
    # We will first set the slope for experiment k.
    beta_tilde <- current_coefs[2, wh_cut, 2] * (sj1 - sj)
    beta_tilde <- beta_tilde + current_coefs[2, wh_cut + 1, 2] * (sj1star - sj1)
    beta_tilde <- beta_tilde / (sj1star - sj)
    
    unif_range <- range(c(beta_tilde, current_coefs[2, wh_cut, 2]))
    u <- runif(1, min = unif_range[1], max = unif_range[2])
    proposed_coefs[2, wh_cut, 2] <- u
    
    # The slope for experiment k + 1 is now deterministic.
    beta_k1 <- delta_k2 - current_coefs[2, wh_cut, 1] - u * (sj1star - sj)
    beta_k1 <- beta_k1 / (sj2 - sj1star)
    proposed_coefs[2, wh_cut + 1, 2] <- beta_k1
    
    # Setting the intercept of experiment k + 1.
    delta_k1 <- proposed_coefs[2, wh_cut, 1]  # Intercept of previous exp.
    delta_k1 <- delta_k1 + proposed_coefs[2, wh_cut, 2] * (sj1star - sj)
    proposed_coefs[2, wh_cut + 1, 1] <- delta_k1
    
    # Calculating the range of the reverse move (for proposal ratio).
    beta_tilde_rev <- delta_k2 - proposed_coefs[2, wh_cut, 1]
    beta_tilde_rev <- beta_tilde_rev - proposed_coefs[2, wh_cut, 2] * (sj1 - sj)
    beta_tilde_rev <- beta_tilde_rev / (sj2 - sj1)
    unif_range_rev <- range(c(beta_tilde_rev, current_coefs[2, wh_cut + 1, 2]))
    
  } else if (sj1star < sj1) {
    
    # We will first set to slope for experiment k + 1.
    beta_tilde <- delta_k2 - current_coefs[2, wh_cut, 1]
    beta_tilde <- beta_tilde - current_coefs[2, wh_cut, 2] * (sj1star - sj)
    beta_tilde <- beta_tilde / (sj2 - sj1star)
    
    unif_range <- range(c(beta_tilde, current_coefs[2, wh_cut + 1, 2]))
    u <- runif(1, min = unif_range[1], max = unif_range[2])
    proposed_coefs[2, wh_cut + 1, 2] <- u
    
    beta_k <- delta_k2 - current_coefs[2, wh_cut, 1]
    beta_k <- (beta_k - u * (sj2 - sj1star)) / (sj1star - sj)
    proposed_coefs[2, wh_cut, 2] <- beta_k
    
    # Setting the intercept of experiment k + 1.
    delta_k1 <- proposed_coefs[2, wh_cut, 1]  # Intercept of previous exp.
    delta_k1 <- delta_k1 + proposed_coefs[2, wh_cut, 2] * (sj1star - sj)
    proposed_coefs[2, wh_cut + 1, 1] <- delta_k1
    
    # Calculating the range of the reverse move (for proposal ratio).
    beta_tilde_rev <- (proposed_coefs[2, wh_cut, 2] * (sj1star - sj) +
                         proposed_coefs[2, wh_cut + 1, 2] * (sj1 - sj1star))
    beta_tilde_rev <- beta_tilde_rev / (sj1 - sj)
    unif_range_rev <- range(c(beta_tilde_rev, proposed_coefs[2, wh_cut, 2]))
  }
  
  return(list(proposed_coefs = proposed_coefs, beta_tilde = beta_tilde,
              unif_range = unif_range, unif_range_rev = unif_range_rev,
              beta_tilde_rev = beta_tilde_rev))
}