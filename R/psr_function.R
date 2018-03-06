#' Calculating the PSR of the LERCA fit.
#' 
#' Based on the posterior samples of the mean ER estimates, calculate the
#' potential scale reduction factor to assess MCMC convergence.
#' 
#' @param samples Matrix. Rows and columns correspond to iterations and chains
#' accordingly.
#' 
#' @return Numeric. Estimate of the potential scale reduction factor. Values
#' close to 1 indicate MCMC convergence.
#' 
#' @export
psr <- function(samples) {
  R <- nrow(samples)
  if (any(is.na(samples))) {
    warning('Some values are NA. na.rm is set to TRUE.') 
  }
  B <- R * var(colMeans(samples, na.rm = TRUE))
  W <- mean(apply(samples, 2, var, na.rm = TRUE))
  value <- sqrt((B + (R - 1) * W) / (R * W))
  return(value)
}
