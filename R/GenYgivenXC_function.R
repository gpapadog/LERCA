#' Generating the outcome given exposure and covariates.
#' 
#' In the simulated data, we specify the outcome model and residual variance
#' for every experiment and the outcome is generated.
#' 
#' @param dataset A data frame containing X, C1, C2, ..., and E.
#' @param out_coef A matrix where the ij element is the coefficient of
#' covariate Ci in the outcome model of experiment j.
#' @param bY Vector of length equal to the number of experiments including the
#' intercept of the outcome model in each experiment.
#' @param bYX If XY_function is set to linear, bYX includes the coefficient of
#' the exposure within each experiment and is of length num_exper. If
#' XY_function is set to other, bYX is of length 1 and corresponds to the
#' coefficient in front of the exposure term.
#' @param Ysd Vector of length number of experiments. Standard deviation of
#' outcome model residual.
#' @param XY_function String specifying whether the XY relationship is piece-
#' wise linear (set 'linear'), or a continuous function supplied by the XY_spec
#' arguement (set 'other').
#' @param XY_spec Needs to be specified if XY_function is set to 'other'. It is
#' the function that specifies the true ER relationship. Defaults to NULL.
#' 
GenYgivenXC <- function(dataset, out_coef, bY, bYX, Ysd, XY_function,
                        XY_spec = NULL) {
  
  dta <- data.table::copy(dataset)

  N <- nrow(dta)
  num_conf <- dim(out_coef)[1]
  num_exper <- length(unique(dta$E))
  
  Ymean <- bY[dta$E]
  if (!is.null(num_conf)) {
    covs <- which(names(dta) %in% paste0('C', 1 : num_conf))
    Ymean <- Ymean + colSums(out_coef[, dta$E] * t(dta[, covs, with = FALSE]))
  }

  if (XY_function == 'linear') {
    Ymean <- Ymean + with(dta, X * bYX[E])
  } else {
    Ymean <- Ymean + with(dta, XY_spec(X) * bYX)
  }

  Y <- Ymean + rnorm(N, mean = 0, sd = Ysd[dta$E])
  return(Y)
}