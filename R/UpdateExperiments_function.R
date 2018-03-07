#' Separate move of the experiment configuration
#' 
#' Updating the experiment configuration separately from the inclusion
#' indicators.
#' 
#' @param dta Data frame including the covariates as C1, C2, ..., the exposure
#' as X and the outcome as Y.
#' @param cov_cols The indices of the columns including the covariates.
#' @param current_cutoffs Numeric of length K. The current values for the
#' points in the experiment configuraiton.
#' @param current_coefs The current coefficients in an array format, with
#' dimensions corresponding to the exposure/outcome model, the experiments, and
#' the coefficient (intercept, slope, covariates).
#' @param current_vars Matrix. Rows correspond to exposure/outcome model, and
#' columns to experiments. Entries are the current variances.
#' @param min_exper_sample The minimum number of observations within an
#' experiment. Defaults to 20.
#' @param prop_distribution Character vector. Options include 'Uniform' or
#' 'Normal' representing the type of distribution that will be used to propose
#' a move of the cutoffs in the separate update. Defaults to uniform.
#' @param normal_percent Numeric. Parameter controling the width of a normal
#' proposal for the experiment configuration. Used only when prop_distribution
#' is set to Normal. Smaller values represent smaller variance of the truncated
#' normal proposal distribution. Defaults to 1.
#' @param mu_priorY Vector of length equal to the number of covariates + 2 with
#' entries corresponding to the prior mean of the intercept, slope, coefficient
#' in the outcome model.
#' @param Sigma_priorY The normal prior covariance matrix of the parameters in
#' mu_priorY.
#' 
#' @return List. Entries are the new, current and proposed experiment
#' configuration, the new coefficients and the indicator of acceptance of the
#' proposed move.
#' 
UpdateExperiments <- function(dta, cov_cols, current_cutoffs, current_coefs,
                              current_vars, min_exper_sample = 20,
                              prop_distribution = c('Uniform', 'Normal'),
                              normal_percent = 1, mu_priorY, Sigma_priorY) {
  
  prop_distribution <- match.arg(prop_distribution)
  minX <- min(dta$X)
  maxX <- max(dta$X)
  num_conf <- ifelse(is.null(cov_cols), 0, length(cov_cols))
  r <- list(cutoffs = current_cutoffs, current_cutoffs = current_cutoffs,
            acc = FALSE, coefs = current_coefs)
  if (min_exper_sample < 1) {
    stop('Set min_exper_sample to 1 or higher.')
  }

  K <- length(current_cutoffs)
  cuts <- c(minX, current_cutoffs, maxX)
  exact_cuts <- c(minX, current_cutoffs, maxX)
  
  
  # ---- STEP 1: Proposing a new value for the experiment configuration ---- #
    
  proposed_cutoffs <- current_cutoffs
  wh_cut <- sample(K, 1)
  choose_from <- c(cuts[wh_cut], cuts[wh_cut + 2])
  
  if (prop_distribution == 'Uniform') {
    proposed_cutoffs[wh_cut] <- runif(1, min = choose_from[1], max = choose_from[2])
  } else {
    # If from a truncated normal, we proposed with mean the current value.
    prop_sd <- normal_percent / 4 * (choose_from[2] - choose_from[1])
    val <- truncnorm::rtruncnorm(1, a = choose_from[1], b = choose_from[2],
                                 mean = proposed_cutoffs[wh_cut],
                                 sd = prop_sd)
    proposed_cutoffs[wh_cut] <- val
  }

  r$proposed_cutoffs <- proposed_cutoffs
  prop_cuts <- c(minX - 0.001, proposed_cutoffs, maxX + 0.001)
  exact_prop_cuts <- c(minX, proposed_cutoffs, maxX)
  
  # If we have less than min_exper_sample data, we reject.

  new_exper <- sapply(dta$X, function(x) sum(x < prop_cuts))
  if (length(table(new_exper)) < K + 1 |
      any(table(new_exper) < min_exper_sample)) {
    return(r)
  }
  
  
  # ---- STEP 2: Proposing slope values for adjacent experiments. ---- #
  
  # Change the slopes of experiments k, k + 1, and intercept of k + 1.
  
  sj <- cuts[wh_cut]
  sj1 <- cuts[wh_cut + 1]
  sj2 <- cuts[wh_cut + 2]
  sj1star <- proposed_cutoffs[wh_cut]
  
  prop_slopes <- ProposeSlopes(wh_cut = wh_cut, current_coefs = current_coefs,
                               cuts = cuts, sj1star = sj1star, sj1 = sj1,
                               sj = sj, sj2 = sj2)
  proposed_coefs <- prop_slopes$proposed_coefs
  unif_range <- prop_slopes$unif_range
  unif_range_rev <- prop_slopes$unif_range_rev
    
  # ---- STEP 3: Calculating the acceptance probability. ---- #
  
  # 3A. Calculating the likelihood ratio.
  
  # 3A - Exposure model: Only data between current and proposed value.
  
  D <- subset(dta, X >= sj1 & X <= sj1star | X >= sj1star & X <= sj1)
  
  # Experiment that data were in before and are proposed to be in.
  curr_exper <- sum(current_cutoffs < D$X[1]) + 1
  prop_exper <- sum(proposed_cutoffs < D$X[1]) + 1
  
  # Log likelihood difference for the exposure model.
  des_mat <- matrix(1, nrow = nrow(D), ncol = 1)
  if (num_conf > 0) {
    des_mat <- cbind(des_mat, as.matrix(D[, cov_cols]))
  }
  mean_prop_exp <- des_mat %*% current_coefs[1, prop_exper, - 2]
  var_prop_exp <- current_vars[1, prop_exper]
  
  mean_curr_exp <- des_mat %*% current_coefs[1, curr_exper, - 2]
  var_curr_exp <- current_vars[1, curr_exper]
  
  logAR <- mvnfast::dmvn(X = D$X, mu = mean_prop_exp, log = TRUE,
                         sigma = var_prop_exp * diag(nrow(D))) -
    mvnfast::dmvn(X = D$X, mu = mean_curr_exp, log = TRUE,
                  sigma = var_curr_exp * diag(nrow(D)))
  
  
  # 3A Outcome model: All data within the two experiments.
  
  D <- subset(dta, X >= sj & X <= sj2)
  D$curr_E <- sapply(D$X, function(x) sum(cuts < x))
  D$prop_E <- sapply(D$X, function(x) sum(prop_cuts < x))
  
  # Proposed likelihood of experiments wh_cut, wh_cut + 1.

  for (ee in c(wh_cut, wh_cut + 1)) {
    
    D_1 <- subset(D, prop_E == ee)
    des_mat <- cbind(1, D_1$X - prop_cuts[ee])
    if (num_conf > 0) {
      des_mat <- cbind(des_mat, as.matrix(D_1[, cov_cols]))
    }
    mean_prop_out <- des_mat %*% proposed_coefs[2, ee, ]
    var_prop_out <- current_vars[2, ee]
    
    logAR <- logAR + mvnfast::dmvn(D_1$Y, mu = mean_prop_out, log = TRUE,
                                   sigma = var_prop_out * diag(nrow(D_1)))
  }
  
  
  # Current_likelihood of experiments wh_cut, wh_cut + 1.
  
  for (ee in c(wh_cut, wh_cut + 1)) {
    
    D_1 <- subset(D, curr_E == ee)
    des_mat <- cbind(1, D_1$X - cuts[ee])
    if (num_conf > 0) {
      des_mat <- cbind(des_mat, as.matrix(D_1[, cov_cols]))
    }
    mean_curr_out <- des_mat %*% current_coefs[2, ee, ]
    var_curr_out <- current_vars[2, ee]
    
    logAR <- logAR - mvnfast::dmvn(D_1$Y, mu = mean_curr_out, log = TRUE,
                                   sigma = var_curr_out * diag(nrow(D_1)))
  }
  
  
  
  # 3B. Calculating the prior ratio.
  
  # For the experiment configuration.
  logAR <- logAR + log(choose_from[2] - proposed_cutoffs[wh_cut]) +
    log(proposed_cutoffs[wh_cut] - choose_from[1])
  logAR <- logAR - log(choose_from[2] - current_cutoffs[wh_cut]) -
    log(current_cutoffs[wh_cut] - choose_from[1])
  
  # For the coefficients that changed (slope of experiment wh_cut, wh_cut + 1).
  prior_sd <- sqrt(Sigma_priorY[2, 2])
  logAR <- logAR +
    sum(dnorm(proposed_coefs[2, c(wh_cut, wh_cut + 1), 2], mean = mu_priorY[2],
              sd = prior_sd, log = TRUE)) -
    sum(dnorm(current_coefs[2, c(wh_cut, wh_cut + 1), 2], mean = mu_priorY[2],
              sd = prior_sd, log = TRUE))
  
  
  # 3C. Calculating the proposal ratio.
  
  # For the cutoffs.
  if (prop_distribution == 'Normal') {
    logAR <- logAR +
      log(truncnorm::dtruncnorm(x = current_cutoffs[wh_cut],
                                a = choose_from[1], b = choose_from[2],
                                mean = proposed_cutoffs[wh_cut],
                                sd = prop_sd)) -
      log(truncnorm::dtruncnorm(x = proposed_cutoffs[wh_cut],
                                a = choose_from[1], b = choose_from[2],
                                mean = current_cutoffs[wh_cut],
                                sd = prop_sd))
  }
  
  # For the coefficients.
  logAR <- logAR + log(unif_range[2] - unif_range[1])
  logAR <- logAR - log(unif_range_rev[2] - unif_range_rev[1])
  
  if (log(runif(1)) < logAR) {
    r$acc <- TRUE
    r$cutoffs <- proposed_cutoffs
    r$coefs <- proposed_coefs
  }
  return(r)
}
