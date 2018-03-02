JumpOverCoef <- function(current_coefs, prop_cuts, cuts, prop_exper_same,
                         curr_exper_same, curr_exper_comb, prop_exper_comb,
                         curr_exper_split, prop_exper_split, tune) {
  
  # We need to set new coefficients.
  proposed_coefs <- current_coefs
  proposed_coefs[2, , 1 : 2] <- NA
  
  # Intercept of first experiment.
  proposed_coefs[2, 1, 1] <- current_coefs[2, 1, 1]
  
  # Intercept and slope of unchanged experiments, intercept of next.
  for (ee in 1 : length(prop_exper_same)) {
    set_exper <- prop_exper_same[ee]
    proposed_coefs[2, set_exper, ] <- current_coefs[2, curr_exper_same[ee], ]
    
    if (set_exper <= K) {
      next_int <- proposed_coefs[2, set_exper, 1]
      interval <- prop_cuts[set_exper + 1] - prop_cuts[set_exper]
      next_int <- next_int + proposed_coefs[2, set_exper, 2] * interval
      proposed_coefs[2, set_exper + 1, 1] <- next_int
    }
  }
  
  # Intercept of combined experiment.
  set_exper <- prop_exper_comb
  proposed_coefs[2, set_exper, 1] <- current_coefs[2, curr_exper_comb[1], 1]
  
  # Intercept of next to combined.
  if (curr_exper_comb[2] < K + 1) {
    next_int <- current_coefs[2, curr_exper_comb[2] + 1, 1]
    proposed_coefs[2, set_exper + 1, 1] <- next_int
  } else {
    next_int <- proposed_coefs[2, set_exper, 1]
    interval <- cuts[curr_exper_comb[1] + 1] - cuts[curr_exper_comb[1]]
    next_int <- next_int + current_coefs[2, curr_exper_comb[1], 2] * interval
    interval <- cuts[curr_exper_comb[2] + 1] - cuts[curr_exper_comb[2]]
    next_int <- next_int + current_coefs[2, curr_exper_comb[2], 2] * interval
  }
  
  # Slope of the combined experiment.
  interval <- prop_cuts[set_exper + 1] - prop_cuts[set_exper]
  slope_set <- (next_int - proposed_coefs[2, set_exper, 1]) / interval
  proposed_coefs[2, set_exper, 2] <- slope_set
  
  # Slope of the first experiment that was split.
  u <- rnorm(1, mean = 0, sd = tune)
  slope_split <- current_coefs[2, curr_exper_split, 2] + u
  proposed_coefs[2, prop_exper_split[1], 2] <- slope_split
  
  # Intercept of the second experiment that was split.
  set_exper <- prop_exper_split[2]
  interval <- prop_cuts[set_exper] - prop_cuts[set_exper - 1]
  next_int <- proposed_coefs[2, set_exper - 1, 1]
  next_int <- next_int + proposed_coefs[2, set_exper - 1, 2] * interval
  proposed_coefs[2, set_exper, 1] <- next_int
  
  # Setting the slope of the second experiment that was split.
  next_int <- current_coefs[2, curr_exper_split, 1]
  interval <- cuts[curr_exper_split + 1] - cuts[curr_exper_split]
  next_int <- next_int + current_coefs[2, curr_exper_split, 2] * interval
  
  interval <- prop_cuts[set_exper + 1] - prop_cuts[set_exper]
  slope_split <- (next_int - proposed_coefs[2, set_exper, 1]) / interval
  proposed_coefs[2, set_exper, 2] <- slope_split
  
  return(list(proposed_coefs = proposed_coefs, u = u))
  
}