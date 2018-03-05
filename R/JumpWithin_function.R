#' Simultaneous move of experiment, alphas and slopes within experiment.
#' 
#' Proposing a simultaneous move of the experiment configuration and inclusion
#' indicators maintaining the order of the experiment configuration. New values
#' for the slopes are proposed ensuring continuous ER. We use the likelihood
#' integrating out the coefficients of the covariates and variance terms.
#' 
#' @param current_cutoffs The current values of the experiment configuration.
#' Vector of length K.
#' @param current_alphas The current values of the inclusion indicators. Array
#' of dimensions 2 (exposure & outcome model), experiments, covariates.
#' @param current_coefs The current values of the coefficients. Array of
#' dimensions 2 (exposure & outcome model), experiments, and coefficients
#' (intercept, slope, covariates).
#' @param cov_cols The indices of the columns including the covariates.
#' @param approx_likelihood Logical. If set to true the BIC will be used to
#' calculate the marginal likelihood. FALSE not supported yet.
#' @param omega The omega parameter of the BAC prior.
#' @param alpha_probs The probability that a proposed alpha is equal to 1, when
#' 0, 1, and 2 alphas of the surrounding experiments are equal to 1. Vector of
#' length 3. Defaults to (0.01, 0.5, 0.99).
#' @param min_exper_sample The minimum number of observations within an
#' experiment. Defaults to 20.
#' 
JumpWithin <- function(dta, current_cutoffs, current_alphas, current_coefs,
                       cov_cols, approx_likelihood = TRUE,
                       approx_jumps = FALSE, omega = 5000, Sigma_priorY,
                       mu_priorY, alpha_probs = c(0.01, 0.5, 0.99),
                       min_exper_sample = 20) {
  
  r <- list(new_cutoffs = current_cutoffs, new_alphas = current_alphas,
            new_coefs = current_coefs, acc = FALSE)
  
  if (!approx_likelihood) {
    stop('approx_likelihood FALSE not supported yet.')
  }
  
  minX <- min(dta$X)
  maxX <- max(dta$X)
  K <- length(current_cutoffs)
  
  # ------ STEP 1. Proposing the cutoffs. -------- #
  
  # Which cutoff we will try to change simultaneously with alphas.
  wh_cut <- sample(K, 1)
  cuts <- c(minX - 0.001, current_cutoffs, maxX + 0.001)
  exact_cuts <- c(minX, current_cutoffs, maxX)

  # Proposed cutoff values.
  new_s <- runif(1, min = exact_cuts[wh_cut], max = exact_cuts[wh_cut + 2])
  proposed_cutoffs <- current_cutoffs
  proposed_cutoffs[wh_cut] <- new_s
  exact_prop_cuts <- c(minX, proposed_cutoffs, maxX)
  prop_cuts <- c(minX - 0.001, proposed_cutoffs, maxX + 0.001)
  
  
  # ------ STEP 2. Proposing the alphas. -------- #
  
  # --- What will the proposed alphas be.
  proposed_alphas <- array(NA, dim = dim(current_alphas))
  dimnames(proposed_alphas) <- dimnames(current_alphas)
  
  # The alphas that remain the same.
  exper_change <- c(wh_cut, wh_cut + 1)
  proposed_alphas[, - exper_change, ] <- current_alphas[, - exper_change, ]
  
  # The alphas that are generated.
  for (mm in 1 : 2) {
    for (jj in 1 : dim(current_alphas)[3]) {
      number_ones <- sum(current_alphas[mm, exper_change, jj])
      probs <- alpha_probs[number_ones + 1]
      probs <- c(1 - probs, probs)
      proposed_alphas[mm, exper_change, jj] <- sample(c(0, 1), 2, prob = probs,
                                                      replace = TRUE)
    }
  }
  
  
  # ------ STEP 3. Proposing the coefficients. -------- #
  
  # Changing the slopes of the new experiments, and the intercept between them.
  
  sj <- exact_cuts[wh_cut]
  sj1 <- exact_cuts[wh_cut + 1]
  sj1star <- exact_prop_cuts[wh_cut + 1]
  sj2 <- exact_cuts[wh_cut + 2]
  
  prop_slopes <- ProposeSlopes(wh_cut = wh_cut, current_coefs = current_coefs,
                               cuts = exact_cuts, sj1star = sj1star, sj1 = sj1,
                               sj = sj, sj2 = sj2)
  proposed_coefs <- prop_slopes$proposed_coefs
  unif_range <- prop_slopes$unif_range
  unif_range_rev <- prop_slopes$unif_range_rev
  
  
  # If the proposed value creates experiment with not sufficient data, reject.
  
  dta$new_exper <- sapply(dta$X, function(x) sum(prop_cuts < x))
  dta$prev_exper <- sapply(dta$X, function(x) sum(cuts < x))
  
  if (length(unique(dta$new_exper)) < K + 1 |
      any(table(dta$new_exper) < min_exper_sample)) {
    return(r)
  }
  
  
  # ----- STEP 4. Calculating the probability of acceptance ------ #
  
  AR <- 0

  # Step 4a. The log Likelihood difference.
  
  for (ee in exper_change) {
    D <- subset(dta, new_exper == ee)
    AR <- AR + LogLike(D = D, curr_exper_alphas = proposed_alphas[, ee, ],
                       curr_coefsY = proposed_coefs[2, ee, 1 : 2],
                       X_s_cut = exact_prop_cuts[ee], cov_cols = cov_cols,
                       approx_jumps = approx_jumps)
  }
  for (ee in exper_change) {
    D <- subset(dta, prev_exper == ee)
    AR <- AR - LogLike(D = D, curr_exper_alphas = current_alphas[, ee, ],
                       curr_coefsY = current_coefs[2, ee, 1 : 2],
                       X_s_cut = exact_cuts[ee], cov_cols = cov_cols,
                       approx_jumps = approx_jumps)
  }
  
  
  # Step 4b. The log prior difference.
  
  # For the cutoffs.
  diff_prop <- sapply(2 : (K + 2), function(x) exact_prop_cuts[x] - exact_prop_cuts[x - 1])
  AR <- AR + sum(log(diff_prop))
  
  diff_curr <- sapply(2 : (K + 2), function(x) exact_cuts[x] - exact_cuts[x - 1])
  AR <- AR - sum(log(diff_curr))
  
  
  # For the alphas.
  prob_alphas <- matrix(c(omega, omega, 1, omega) / (3 * omega + 1),
                        nrow = 2, ncol = 2, byrow = TRUE)
  
  # For the proposed alpha values.
  tab_alpha_prop <- table(proposed_alphas[1, exper_change, ],
                          proposed_alphas[2, exper_change, ])
  wh_rows <- as.numeric(rownames(tab_alpha_prop)) + 1
  wh_cols <- as.numeric(colnames(tab_alpha_prop)) + 1
  AR <- AR + sum(tab_alpha_prop * log(prob_alphas)[wh_rows, wh_cols])
  
  # For the current alpha values.
  tab_alphas_curr <- table(current_alphas[1, exper_change, ],
                           current_alphas[2, exper_change, ])
  wh_rows <- as.numeric(rownames(tab_alphas_curr)) + 1
  wh_cols <- as.numeric(colnames(tab_alphas_curr)) + 1
  AR <- AR - sum(tab_alphas_curr * log(prob_alphas)[wh_rows, wh_cols])
  
  
  # For the slopes that changed.
  if (!approx_jumps) {
    prior_sd <- sqrt(Sigma_priorY[2, 2])
    AR <- AR +
      sum(dnorm(proposed_coefs[2, c(wh_cut, wh_cut + 1), 2], mean = mu_priorY[2],
                sd = prior_sd, log = TRUE)) -
      sum(dnorm(current_coefs[2, c(wh_cut, wh_cut + 1), 2], mean = mu_priorY[2],
                sd = prior_sd, log = TRUE))
  }

  
  # Step 4c. The log proposal difference. (The cutoffs' proposal is symmetric.)

  # For the alphas.
  proposal_probs <- sapply(alpha_probs,
                           function(x) c((1 - x) ^ 2, x * (1 - x), x ^ 2))
  proposal_probs <- t(proposal_probs)
  dimnames(proposal_probs) <- list(current = 0 : 2, proposed = 0 : 2)
  
  sum_curr_alphas <- apply(current_alphas[, exper_change, ], c(1, 3), sum)
  sum_prop_alphas <- apply(proposed_alphas[, exper_change, ], c(1, 3), sum)
  
  cross_sum <- table(sum_curr_alphas, sum_prop_alphas)
  wh_rows <- as.numeric(rownames(cross_sum)) + 1
  wh_cols <- as.numeric(colnames(cross_sum)) + 1
  
  
  # Probability of proposing current alphas from the proposed.
  AR <- AR - sum(t(cross_sum) * log(proposal_probs)[wh_cols, wh_rows])
 
  # Probability of proposing proposed from current.
  AR <- AR - sum(cross_sum * log(proposal_probs)[wh_rows, wh_cols])
  
  
  # For the slopes.
  if (!approx_jumps) {
    AR <- AR + log(unif_range[2] - unif_range[1])
    AR <- AR - log(unif_range_rev[2] - unif_range_rev[1])
  }

  
  # ------- Accepting or rejecting the move. -------- #

  if (log(runif(1)) < AR) {
    r$new_cutoffs <- proposed_cutoffs
    r$new_alphas <- proposed_alphas
    r$acc <- TRUE
    r$new_coefs <- proposed_coefs
  }
  
  return(r)
}


