#' Acquire the coefficients as a function of the exposure.
#' 
#' LERCA coefficients in each experiment are translated to coefficients as a
#' function of the exposure.
#' 
#' @param lerca The LERCA fit.
#' @param exp_values A vector of exposure values.
#' 
#' @export
ExposureCoefs <- function(lerca, exp_values) {
  
  chains <- dim(lerca$coefs)[2]
  num_points <- length(exp_values)
  num_vars <- dim(lerca$coefs)[5]
  num_conf <- num_vars - 2
  post_samples <- dim(lerca$coefs)[3]
  
  exp_coef <- array(0, dim = c(chains, post_samples, 2, num_points, num_vars))
  dimnames(exp_coef) <- list(chain = 1 : chains, sample = 1 : post_samples,
                             model = c('Exp', 'Out'), exposure = exp_values,
                             var = c('Int', 'X', paste0('C', 1 : num_conf)))

  for (ii in 1 : length(exp_values)) { # For every point on the exposure range
    for (cc in 1 : chains) {
      for (kk in 1 : post_samples) {
        exper <- sum(exp_values[ii] > lerca$cutoffs[cc, kk, ]) + 1
        exp_coef[cc, kk, , ii, ] <- lerca$coefs[, cc, kk, exper, ]
      }
    }
  }
  
  return(exp_coef)
}