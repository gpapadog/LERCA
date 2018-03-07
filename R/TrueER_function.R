#' Calculating the true ER of simulated data set.
#' 
#' Based on a simulated data set and known experiment configuration and outcome
#' model coefficients calculate the exposure response curve.
#' 
#' @param dta Data frame. Includes exposure as 'X', outcome as 'Y' and
#' covariates as C1, C2, ...
#' @param true_cutoffs Numeric vector. The true points of the experiment
#' configuration.
#' @param out_coefs Matrix. Rows correspond to experiments and columns to
#' coefficients (intercept, slope, covariates) in the outcome model.
#' @param predict_at The values of the exposure we want to predict the response
#' at. If left NULL, specify grid_length.
#' @param grid_length The number of exposure points we want to estimate the
#' mean response at. If predict_at is left NULL, an equally-distanced grid of
#' values of length grid_length over the observed exposure range will be used.
#' Defaults to 100.
#' @param XY_function The true ER shape. Options are 'linear' and 'other'.
#' Defaults to 'linear'.
#' @param XY_spec Function. If XY_function is set to 'other' specify the true
#' ER shape in XY_spec. Leave NULL otherwise.
#' 
#' @return List of two elements. The first one named 'x' is a vector of the 
#' exposure values at which we evaluated the true ER. The second one named 'y'
#' is a matrix of rows equal to the number of exposure values in 'x', and
#' columns equal to the number of observations, including the expected response
#' of an observation at a specific exposure value.
#' 
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
