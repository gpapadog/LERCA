#' @param current_coefs The current coefficients in an array format, with
#' dimensions corresponding to the exposure/outcome model, the experiments, and
#' the coefficient (intercept, slope, covariates).
#' @param approximate Logical. If set to true the BIC will be used to calculate
#' the marginal likelihood. FALSE not supported yet.
#' @param cov_cols The indices of the columns including the covariates.
#' @param comb_probs When two experiments are combined, comb_probs represents
#' the probability of alpha = 1 when 0, 1, and 2 corresponding alphas are equal
#' to 1. Vector of length 3.
#' @param split_probs When one experiment is split, split_probs describes the
#' probability that the alpha of a new experiment is equal to 1, when the alpha
#' of the current experiment is 0, and when it is 1. Vector of length 2.
#' @param min_exper_sample The minimum number of observations within an
#' experiment.
JumpOver <- function(dta, current_cutoffs, current_alphas, current_coefs,
                     approximate = TRUE, cov_cols, omega = 5000,
                     mu_priorY, Sigma_priorY,
                     comb_probs = c(0.01, 0.5, 0.99), tune = 0.05,
                     split_probs = c(0.2, 0.95), min_exper_sample = 20) {
  
  r <- list(new_cutoffs = current_cutoffs, new_alphas = current_alphas,
            new_coefs = current_coefs, acc = FALSE,
            current_cutoffs = current_cutoffs)
  
  if (!approximate) {
    stop('approximate FALSE not supported yet.')
  }
  
  minX <- min(dta$X)
  maxX <- max(dta$X)
  K <- length(current_cutoffs)
  
  # ------ STEP 1 : Choosing the cutoff that will be moved. ------- #
  
  wh_cut <- sample(K, 1)
  cuts <- c(minX, current_cutoffs, maxX)
  
  # --- Which experiment the value will get moved to.
  
  # Experiment choices with probability.
  exper_choice <- (1 : (K + 1))[- c(wh_cut, wh_cut + 1)]
  prob_exper <- sapply(2 : length(cuts), function(x) cuts[x] - cuts[x - 1])
  prob_exper <- prob_exper[- c(wh_cut, wh_cut + 1)]
  # Cutoff will be moved to experiment:
  split_exper <- sample(exper_choice, 1, prob = prob_exper)
  new_s <- runif(1, min = cuts[split_exper], max = cuts[split_exper + 1])
  proposed_cutoffs <- sort(c(new_s, current_cutoffs[- wh_cut]))
  prop_cuts <- c(minX, proposed_cutoffs, maxX)
  

  # --- If the proposed value creates experiment with no data, do not accept.

  dta$new_exper <- sapply(dta$X, function(x) sum(prop_cuts < x))
  dta$prev_exper <- sapply(dta$X, function(x) sum(cuts < x))
  dta$new_exper[dta$new_exper == 0] <- 1
  dta$prev_exper[dta$prev_exper == 0] <- 1
  
  if (length(unique(dta$new_exper)) < K + 1 |
      any(table(dta$new_exper) < min_exper_sample)) {
    return(r)
  }
  
  
  # ------- STEP 2: Proposing new alphas for the experiments. --------- #
  
  
  # --- What will the proposed alphas be.
  proposed_alphas <- array(NA, dim = dim(current_alphas))
  dimnames(proposed_alphas) <- dimnames(current_alphas)
  
  # The alphas that remain the same.
  curr_exper_same <- setdiff(exper_choice, split_exper) # Prev exper number.
  start_same <- cuts[curr_exper_same] # Starting cutoff.
  prop_exper_same <- which(prop_cuts %in% start_same) # Proposed experiment.
  proposed_alphas[, prop_exper_same, ] <- current_alphas[, curr_exper_same, ]
  
  # The alphas that get split.
  which_new_s <- which(prop_cuts == new_s)
  prop_exper_split <- c(which_new_s - 1, which_new_s)
  curr_exper_split <- which(cuts == prop_cuts[which_new_s - 1])
  
  for (ee in prop_exper_split) {
    for (mm in 1 : 2) {
      for (jj in 1 : dim(current_alphas)[3]) {
        probs <- split_probs[current_alphas[mm, curr_exper_split, jj] + 1]
        probs <- c(1 - probs, probs)
        proposed_alphas[mm, ee, jj] <- sample(c(0, 1), 1, prob = probs)
      }
    }
  }
  
  # The alphas that get combined.
  start_comb <- cuts[wh_cut]
  curr_exper_comb <- c(wh_cut, wh_cut + 1)
  prop_exper_comb <- which(prop_cuts == start_comb) # Which one is the combined
  combine_alphas <- current_alphas[, curr_exper_comb, ]
  sum_combine <- apply(combine_alphas, c(1, 3), sum)
  for (mm in 1 : 2) {
    for (jj in 1 : ncol(sum_combine)) {
      probs <- comb_probs[sum_combine[mm, jj] + 1]
      probs <- c(1 - probs, probs)
      proposed_alphas[mm, prop_exper_comb, jj] <- sample(c(0, 1), 1,
                                                         prob = probs)
    }
  }
  
  
  # ------- STEP 3: Proposing coefficients for the new experiments. ------- #
  
  coef_prop <- JumpOverCoef(current_coefs = current_coefs,
                            prop_cuts = prop_cuts, cuts = cuts,
                            prop_exper_same = prop_exper_same,
                            curr_exper_same = curr_exper_same,
                            curr_exper_comb = curr_exper_comb,
                            prop_exper_comb = prop_exper_comb,
                            curr_exper_split = curr_exper_split,
                            prop_exper_split = prop_exper_split,
                            tune = tune)
  proposed_coefs <- coef_prop$proposed_coefs
  slope_u <- coef_prop$u
  slope_u_rev <- coef_prop$u_rev
  
  
  # ------- STEP 4. Calculating the probability of acceptance ------ #
  
  AR <- 0
  
  # Step 4a. The log-Likelihood difference.
  
  # The new experiments are prop_exper_split and prop_exper_comb.
  for (ee in c(prop_exper_split, prop_exper_comb)) {
    D <- subset(dta, new_exper == ee)
    AR <- AR + LogLike(D = D, curr_exper_alphas = proposed_alphas[, ee, ],
                       curr_coefsY = proposed_coefs[2, ee, 1 : 2],
                       X_s_cut = prop_cuts[ee], cov_cols = cov_cols)
  }
  for (ee in c(curr_exper_split, curr_exper_comb)) {
    D <- subset(dta, prev_exper == ee)
    AR <- AR - LogLike(D = D, curr_exper_alphas = current_alphas[, ee, ],
                       curr_coefsY = current_coefs[2, ee, 1 : 2],
                       X_s_cut = cuts[ee], cov_cols = cov_cols)
  }
  
  
  # Step 4b. The log prior difference.
  
  # For the cutoffs.
  diff_prop <- sapply(2 : (K + 2), function(x) prop_cuts[x] - prop_cuts[x - 1])
  AR <- AR + sum(log(diff_prop))
  
  diff_curr <- sapply(2 : (K + 2), function(x) cuts[x] - cuts[x - 1])
  AR <- AR - sum(log(diff_curr))
  
  
  # For the alphas.
  prob_alphas <- matrix(c(omega, omega, 1, omega) / (3 * omega + 1),
                        nrow = 2, ncol = 2, byrow = TRUE)
  
  # For the proposed values.
  prop_changed_exper <- c(prop_exper_comb, prop_exper_split)
  tab_alpha_prop <- table(proposed_alphas[1, prop_changed_exper, ],
                          proposed_alphas[2, prop_changed_exper, ])
  wh_rows <- as.numeric(rownames(tab_alpha_prop)) + 1
  wh_cols <- as.numeric(colnames(tab_alpha_prop)) + 1
  AR <- AR + sum(tab_alpha_prop * log(prob_alphas)[wh_rows, wh_cols])
  
  # For the current values.
  curr_changed_exper <- c(curr_exper_comb, curr_exper_split)
  tab_alphas_curr <- table(current_alphas[1, curr_changed_exper, ],
                           current_alphas[2, curr_changed_exper, ])
  wh_rows <- as.numeric(rownames(tab_alphas_curr)) + 1
  wh_cols <- as.numeric(colnames(tab_alphas_curr)) + 1
  AR <- AR - sum(tab_alphas_curr * log(prob_alphas)[wh_rows, wh_cols])
  
  
  # For the slopes that changed.
  prior_sd <- sqrt(Sigma_priorY[2, 2])
  AR <- AR +
    sum(dnorm(proposed_coefs[2, , 2], mean = mu_priorY[2],
              sd = prior_sd, log = TRUE)) -
    sum(dnorm(current_coefs[2, , 2], mean = mu_priorY[2],
              sd = prior_sd, log = TRUE))
  
  
  
  # Step 4c. The log proposal difference.
  
  # For the cutoffs.
  
  # Probability of proposing current state from the proposed.
  prop_exper_size <- sapply(2 : length(prop_cuts),
                            function(x) prop_cuts[x] - prop_cuts[x - 1])
  # Proposal would be uniform over all but prop_exper_split experiments.
  AR <- AR + log(1 / sum(prop_exper_size[- prop_exper_split]))
  
  # Probability of proposing proposal from current.
  curr_exper_size <- sapply(2 : length(cuts),
                            function(x) cuts[x] - cuts[x - 1])
  AR <- AR - log(1 / sum(curr_exper_size[- curr_exper_comb]))
  
  
  # For the alphas.
  
  # Probability of proposing current alphas from the proposed.
  # For the experiment that got split.
  prop_alpha_split <- apply(proposed_alphas[, prop_exper_split, ], c(1, 3), sum)
  curr_alpha_split <- current_alphas[, curr_exper_split, ]
  tab_alpha_split <- table(curr_alpha_split, prop_alpha_split)
  wh_rows <- as.numeric(rownames(tab_alpha_split)) + 1
  wh_cols <- as.numeric(colnames(tab_alpha_split)) + 1
  
  probs_alpha_split <- rbind(1 - comb_probs, comb_probs)
  AR <- AR + sum(tab_alpha_split * log(probs_alpha_split)[wh_rows, wh_cols])
  
  # For the experiments that got combined.
  curr_alpha_comb <- apply(current_alphas[, curr_exper_comb, ], c(1, 3), sum)
  prop_alpha_comb <- proposed_alphas[, prop_exper_comb, ]
  tab_alpha_comb <- table(curr_alpha_comb, prop_alpha_comb)
  wh_rows <- as.numeric(rownames(tab_alpha_comb)) + 1
  wh_cols <- as.numeric(colnames(tab_alpha_comb)) + 1
  
  probs_alpha_comb <-
    sapply(split_probs, function(x) c((1 - x) ^ 2, x * (1 - x), x ^ 2))
  AR <- AR + sum(tab_alpha_comb * log(probs_alpha_comb)[wh_rows, wh_cols])
  
  # Probability of proposing proposed from current.
  wh_rows <- as.numeric(rownames(t(tab_alpha_split))) + 1
  wh_cols <- as.numeric(colnames(t(tab_alpha_split))) + 1
  AR <- AR - sum(t(tab_alpha_split) * log(probs_alpha_comb)[wh_rows, wh_cols])
  
  wh_rows <- as.numeric(rownames(t(tab_alpha_comb))) + 1
  wh_cols <- as.numeric(colnames(t(tab_alpha_comb))) + 1
  AR <- AR - sum(t(tab_alpha_comb) * log(probs_alpha_split)[wh_rows, wh_cols])
  
  
  # For the slopes.
  AR <- AR + dnorm(slope_u_rev, mean = 0, sd = tune, log = TRUE)
  AR <- AR - dnorm(slope_u, mean = 0, sd = tune, log = TRUE)
  
  
  # ------- Accepting or rejecting the move. -------- #
  acc <- (log(runif(1)) < AR)
  if (!acc) {
    return(r)
  }
  
  # Else we accepted the move.
  
  r$new_cutoffs <- proposed_cutoffs
  r$new_alphas <- proposed_alphas
  r$acc <- TRUE
  r$new_coefs <- proposed_coefs
  return(r)
}


