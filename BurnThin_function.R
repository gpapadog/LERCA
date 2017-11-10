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