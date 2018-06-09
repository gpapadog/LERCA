#' Acquire the inclusions probabilities as a function of the exposure.
#' 
#' LERCA inclusion indicators in each experiment are translated to inclusion
#' probabilities as a function of the exposure.
#' 
#' @param lerca The LERCA fit.
#' @param exp_values A vector of exposure values where the inclusion
#' probabilities should be calculated.
#' 
#' @export
ExposureInclusion <- function(lerca, exp_values) {
  
  chains <- dim(lerca$alphas)[2]
  num_points <- length(exp_values)
  num_conf <- dim(lerca$alphas)[5]
  post_samples <- dim(lerca$alphas)[3]
  
  inclusion <- array(0, dim = c(chains, 2, num_points, num_conf))
  dimnames(inclusion) <- list(chain = 1 : chains, model = c('Exp', 'Out'),
                              exposure = exp_values,
                              conf = paste0('C', 1 : num_conf))

  for (ii in 1 : length(exp_values)) {
    for (cc in 1 : chains) {
      for (kk in 1 : post_samples) {
        exper <- sum(exp_values[ii] > lerca$cutoffs[cc, kk, ]) + 1
        this_sample <- lerca$alphas[, cc, kk, exper, ]
        inclusion[cc, , ii, ] <- inclusion[cc, , ii, ] + this_sample
      }
    }
  }
  inclusion <- inclusion / post_samples
  
  return(inclusion)
}