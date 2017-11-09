UpdateExperiments <- function(dta, cov_cols, current_cutoffs, current_coefs,
                              current_vars,
                              prop_distribution = c('Normal', 'Uniform'),
                              normal_percent = 1) {
  
  minX <- min(dta$X)
  maxX <- max(dta$X)
  
  prop_distribution <- match.arg(prop_distribution)
  K <- length(current_cutoffs)
  proposed_cutoffs <- current_cutoffs
  wh_cut <- sample(K, 1)
  cuts <- c(minX, current_cutoffs, maxX)
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
  
  # Calculating the AR.
  cuts <- c(minX - 0.0001, current_cutoffs, maxX + 0.0001)
  
  sj <- cuts[wh_cut]
  sj1 <- cuts[wh_cut + 1]
  sj2 <- cuts[wh_cut + 2]
  sj1star <- proposed_cutoffs[wh_cut]
  
  D <- subset(dta, X >= sj1 & X <= sj1star | X >= sj1star & X <= sj1)
  
  # If no observations changed experiment, do not accept.
  if (nrow(D) == 0) {
    return(list(cutoffs = current_cutoffs, acc = FALSE))
    cutoffs[cc, ii, ] <- current_cutoffs
  }
  
  # Otherwise, observations changed experiment, calcualting the AR.
  
  curr_exper <- sum(current_cutoffs < D$X[1]) + 1
  prop_exper <- sum(proposed_cutoffs < D$X[1]) + 1
  
  # --- Likelihood ratio.
  
  # For the outcome.
  des_mat <- cbind(1, D$X, as.matrix(D[, cov_cols]))
  mean_prop_out <- des_mat %*% current_coefs[2, prop_exper, ]
  var_prop_out <- current_vars[2, prop_exper]
  
  mean_curr_out <- des_mat %*% current_coefs[2, curr_exper, ]
  var_curr_out <- current_vars[2, curr_exper]
  
  # For the exposure.
  des_mat <- cbind(1, as.matrix(D[, cov_cols]))
  mean_prop_exp <- des_mat %*% current_coefs[1, prop_exper, - 2]
  var_prop_exp <- current_vars[1, prop_exper]
  
  mean_curr_exp <- des_mat %*% current_coefs[1, curr_exper, - 2]
  var_curr_exp <- current_vars[1, curr_exper]
  
  AR <- mvnfast::dmvn(D$Y, mean_prop_out, sigma = var_prop_out * diag(nrow(D)),
                      log = TRUE) -
    mvnfast::dmvn(D$Y, mean_curr_out, sigma = var_curr_out * diag(nrow(D)),
                  log = TRUE) +
    mvnfast::dmvn(D$X, mean_prop_exp, sigma = var_prop_exp * diag(nrow(D)),
                  log = TRUE) -
    mvnfast::dmvn(D$X, mean_curr_exp, sigma = var_curr_exp * diag(nrow(D)),
                  log = TRUE)
  
  # --- Prior ratio.
  
  AR <- AR + log(choose_from[2] - proposed_cutoffs[wh_cut]) +
    log(proposed_cutoffs[wh_cut] - choose_from[1])
  AR <- AR - log(choose_from[2] - current_cutoffs[wh_cut]) -
    log(current_cutoffs[wh_cut] - choose_from[1])
  
  
  # --- Proposal ratio.
  
  if (prop_distribution == 'Normal') {
    AR <- AR +
      log(truncnorm::dtruncnorm(x = current_cutoffs[wh_cut],
                                a = choose_from[1], b = choose_from[2],
                                mean = proposed_cutoffs[wh_cut],
                                sd = prop_sd)) -
      log(truncnorm::dtruncnorm(x = proposed_cutoffs[wh_cut],
                                a = choose_from[1], b = choose_from[2],
                                mean = current_cutoffs[wh_cut],
                                sd = prop_sd))
  }
  
  if (log(runif(1)) < AR) {
    return(list(cutoffs = proposed_cutoffs, acc = TRUE))
  } else {
    return(list(cutoffs = current_cutoffs, acc = FALSE))
  }
}
