LogLike <- function(D, curr_exper_alphas, cov_cols) {
  
  cov_matrix <- as.matrix(D[, cov_cols])
  
  # For the exposure model.
  design_mat <- cbind(1, cov_matrix[, curr_exper_alphas[1, ] == 1])
  lmod <- lm(D$X ~ design_mat - 1)
  log_like <- - BIC(lmod) / 2
  
  # For the outcome model.
  design_mat <- cbind(1, cov_matrix[, curr_exper_alphas[2, ] == 1])
  lmod <- lm(D$Y ~ D$X + design_mat - 1)
  log_like <- log_like - BIC(lmod) / 2
  
  return(log_like)
}