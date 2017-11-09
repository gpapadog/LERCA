BurnThin <- function(lerca, burn, thin) {
  
  if (burn == 0 & thin == 1) {
    return(lerca)
  }
  
  keep <- seq(burn + 1, dim(lerca$cutoffs)[2], by = thin)
  
  lerca$cutoffs <- lerca$cutoffs[, keep, ]
  lerca$alphas <- lerca$alphas[, , keep, , ]
  lerca$coefs <- lerca$coefs[, , keep, , ]
  lerca$variances <- lerca$variances[, , keep, ]
  
  return(lerca)
}