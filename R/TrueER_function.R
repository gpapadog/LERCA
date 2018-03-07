#' @export
TrueER <- function(dta, true_cutoffs, out_coefs, predict_at = NULL,
                   grid_length = 100, XY_function = c('linear', 'other'),
                   XY_spec = NULL) {
  
  XY_function <- match.arg(XY_function)
  
  minX <- min(dta$X)
  maxX <- max(dta$X)

  if (is.null(predict_at)) {
    predict_at <- seq(minX, maxX, length.out = grid_length)
  }
  
  num_conf <- dim(out_coefs)[2] - 2
  des_mat <- cbind(Int = 1, X = dta$X)
  if (num_conf != 0) {
    covs_cols <- which(names(dta) %in% paste0('C', 1 : num_conf))
    des_mat <- cbind(des_mat, dta[, covs_cols, with = FALSE])
  }
  des_mat <- as.matrix(des_mat)
  
  # Array where the predicted counterfactual values are saved.
  counter <- array(NA, dim = c(length(predict_at), nrow(dta)))
  dimnames(counter) <- list(predict = predict_at, obs = 1 : nrow(dta))
  
  for (pp in 1:length(predict_at)) {
    # Find the experiment that the predicting exposure is in.
    exper <- sum(true_cutoffs <= predict_at[pp])
    D <- des_mat
    if (XY_function == 'linear') {
      D[, 2] <- predict_at[pp]
    } else {
      D[, 2] <- XY_spec(predict_at[pp])
    }
    counter[pp, ] <- D %*% out_coefs[exper, ]
  }
  return(list(x = predict_at, y = counter))
}
