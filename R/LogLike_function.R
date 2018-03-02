#' Approximate log likelihood of an experiment.
#' 
#' Calculating the log likelihood based on the BIC approximation. We use X
#' instead of X - s in the outcome model since it provides the same BIC value.
LogLike <- function(D, curr_exper_alphas, curr_coefsY, X_s_cut, cov_cols) {
  
  cov_matrix <- as.matrix(D[, cov_cols])
  
  # For the exposure model.
  des_mat <- cbind(1, matrix(cov_matrix[, curr_exper_alphas[1, ] == 1],
                             nrow = nrow(D)))
  lmod <- lm(D$X ~ des_mat - 1)
  log_like <- - BIC(lmod) / 2
  
  # For the outcome model.
  des_mat <- matrix(cov_matrix[, curr_exper_alphas[2, ] == 1], nrow = nrow(D))
  resid <- D$Y - cbind(1, D$X - X_s_cut) %*% matrix(curr_coefsY, ncol = 1)
  if (ncol(des_mat) > 0) {
    lmod <- lm(resid ~ des_mat - 1)
  } else {
    lmod <- lm(resid ~ - 1)
  }
  log_like <- log_like - BIC(lmod) / 2
  
  return(log_like)
}