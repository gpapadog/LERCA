#' WAIC of LERCA fit using the full conditionals.
#'
#' @param lerca The LERCA output including cutoffs, coefficients and variances
#' of the exposure and outcome models.
#' @param dta The full data set including the confounders as C1, C2, ..., the
#' outcome as Y and the exposure as X.
#' 
#' @export
WAIC <- function(lerca, dta) {
  
  dta <- as.data.frame(dta)
  num_conf <- dim(lerca$alphas)[5]
  covs_cols <- which((names(dta) %in% paste0('C', 1 : num_conf)))
  N <- nrow(dta)
  keep <- lerca$keep
  
  # Values of mu_prior and Sigma_prior if NULL.
  if (is.null(mu_priorX)) {
    mu_priorX <- rep(0, num_conf + 1)
  }
  if (is.null(mu_priorY)) {
    mu_priorY <- rep(0, num_conf + 2)
  }
  if (is.null(Sigma_priorX)) {
    Sigma_priorX <- diag(rep(100 ^ 2, num_conf + 1))
  }
  if (is.null(Sigma_priorY)) {
    Sigma_priorY <- diag(rep(100 ^ 2, num_conf + 2))
  }
  
  
  lppd <- 0
  Elog <- 0
  logE <- 0
  pwaic1 <- 0
  pwaic2 <- 0
  
  for (ii in 1 : N) {
    
    log_like <- matrix(NA, nrow = length(keep), ncol = chains)
    
    for (ss in 1 : length(keep)) {
      for (cc in 1 : chains) {
        
        curr_exper <- sum(lerca$cutoffs[cc, keep[ss], ] < dta$X[ii]) + 1
        curr_alphaX <- lerca$alphas[cc, 1, keep[ss], curr_exper, ]
        curr_alphaY <- lerca$alphas[cc, 2, keep[ss], curr_exper, ]
        
        
        log_like[ss, cc] <- CalcLogLikeXandY(D = dta[ii, ],
                                             alphaX = curr_alphaX,
                                             alphaY = curr_alphaY,
                                             covs_cols = covs_cols,
                                             nu_priorX = nu_priorX,
                                             nu_priorY = nu_priorY,
                                             lambda_priorX = lambda_priorX,
                                             lambda_priorY = lambda_priorY,
                                             mu_priorX = mu_priorX,
                                             mu_priorY = mu_priorY, 
                                             Sigma_priorX = Sigma_priorX,
                                             Sigma_priorY = Sigma_priorY,
                                             approximate = FALSE)
      }
    }
    
    logE_ii <- log(mean(exp(log_like)))
    Elog_ii <- mean(log_like)
    
    lppd <- lppd + logE_ii
    pwaic1 <- pwaic1 + 2 * (logE_ii - Elog_ii)
    pwaic2 <- pwaic2 + var(as.numeric(log_like))
  }
  r <- c(lppd = lppd, pwaic1 = pwaic1, pwaic2 = pwaic2,
         waic1 = - 2 * (lppd - pwaic1), waic2 = - 2 * (lppd - pwaic2))
  return(r)
}
