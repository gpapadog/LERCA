#' Approximate log likelihood of an experiment.
#' 
#' Calculating the log likelihood based on the BIC approximation.
#' 
#' @param D The data set of the current experiment includes covariates,
#' exposure as 'X' and outcome as 'Y'.
#' @param curr_exper_alphas Matrix. Dimensions correspond to exposure/outcome
#' model and potential confounders. Entries 0/1 represent exlusion/inclusion of
#' the covariate in each model.
#' @param curr_coefsY Numeric of length two. Intercept and slope of the outcome
#' model. If left NULL, the likelihood is calculated integrating out intercept
#' and slope along with remaining coefficients.
#' @param X_s_cut Numeric. The point in the experiment configuration that
#' corresponds to the beginning of the current experiment.
#' @param cov_cols The indices of the columns in D that correspond to the
#' potential confounders.
#' 
LogLike <- function(D, curr_exper_alphas, curr_coefsY = NULL, X_s_cut,
                    cov_cols) {
  
  approx_jumps <- is.null(curr_coefsY)
  
  cov_matrix <- as.matrix(D[, cov_cols])
  
  # For the exposure model.
  des_mat <- cbind(1, matrix(cov_matrix[, curr_exper_alphas[1, ] == 1],
                             nrow = nrow(D)))
  lmod <- lm(D$X ~ des_mat - 1)
  log_like <- - BIC(lmod) / 2
  
  # For the outcome model.
  des_mat <- matrix(cov_matrix[, curr_exper_alphas[2, ] == 1], nrow = nrow(D))
  if (approx_jumps) {
    des_mat <- cbind(1, D$X, des_mat)
    lmod <- lm(D$Y ~ des_mat)
  } else {
    resid <- D$Y - cbind(1, D$X - X_s_cut) %*% matrix(curr_coefsY, ncol = 1)
    if (ncol(des_mat) > 0) {
      lmod <- lm(resid ~ des_mat - 1)
    } else {
      lmod <- lm(resid ~ - 1)
    }
    
  }
  log_like <- log_like - BIC(lmod) / 2
  
  return(log_like)
}