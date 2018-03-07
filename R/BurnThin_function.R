#' Burn in and thinning of the LERCA fit.
#' 
#' @param lerca The LERCA fit.
#' @param burn The number of initial iterations to discard.
#' @param thin The distance of iterations we want to keep.
#' 
#' @export
#' 
#' @example
#' lerca_short <- BurnThin(lerca = lerca, burn = 1000, thin = 10)
#' 
BurnThin <- function(lerca, burn, thin) {
  
  if (burn == 0 & thin == 1) {
    return(lerca)
  }
  
  keep <- seq(burn + 1, dim(lerca$cutoffs)[2], by = thin)
  
  lerca$cutoffs <- lerca$cutoffs[, keep, , drop = FALSE]
  lerca$alphas <- lerca$alphas[, , keep, , , drop = FALSE]
  lerca$coefs <- lerca$coefs[, , keep, , , drop = FALSE]
  lerca$variances <- lerca$variances[, , keep, , drop = FALSE]
  
  return(lerca)
}