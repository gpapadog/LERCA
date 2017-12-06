#' CER from MCMC of multiple chains.
#' 
#' When multiple chains are ran, we use GetER to get the posterior CER samples.
#' 
#' @param dta The dataset including columns X (treatment), Y (outcome), and
#' potential confounders named as C1, C2, ...
#' @param cutoffs A 3-dimensional array with dimensions corresponding to chain,
#' iteration, and number of points in the experiment configuration.
#' @param coefs A 4-dimensional array with dimensions corresponding to chain,
#' iteration, experiment, and variable. Variables are intercept, exposure,
#' and the potential confounders.
#' @param predict_at A vector of the exposure values we wish to predict ER at.
#' If left NULL, a grid of length grid_length will be created.
#' @param grid_length The number of locations we will predict the ER. Defaults
#' to 100.
#' @param mean_only Logical. Set to FALSE if we want the individual ER
#' predictions. Set to TRUE if we are only interested in the posterior samples
#' of the mean ER. Defaults to FALSE.
#' @param other_function Whether we want to apply a different function to the
#' predictions. Defaults to NULL. Examples include exp(x), log(x).
#' 
#' @export
GetER <- function(dta, cutoffs, coefs, predict_at = NULL, grid_length = 100,
                  mean_only = FALSE, other_function = NULL) {
  
  minX <- min(dta$X)
  maxX <- max(dta$X)
  
  if (any(dim(cutoffs) + c(0, 0, 1) != dim(coefs)[- 4])) {
    stop('cutoffs and coefs are not of appropriate dimension.')
  }
  chains <- dim(cutoffs)[1]
  Nsims <- dim(cutoffs)[2]
  K <- dim(cutoffs)[3]
  
  if (is.null(predict_at)) {
    predict_at <- seq(minX, maxX, length.out = grid_length)
  }
  if (any(predict_at < minX) | any(predict_at > maxX)) {
    warning(paste0('Exposure values to predict at outside the exposure range',
                   'will be ignored.'))
  }
  
  if (mean_only) {  # We keep post samples of the mean ER.
    counter <- array(NA, dim = c(length(predict_at), chains, Nsims))
    dimnames(counter) <- list(predict = predict_at, chain = 1 : chains,
                              sim = 1 : Nsims)
  } else {
    counter <- array(NA, dim = c(length(predict_at), nrow(dta), chains, Nsims))
    dimnames(counter) <- list(predict = predict_at, obs = 1 : nrow(dta),
                              chain = 1 : chains, sim = 1 : Nsims)
  }
  if (!is.null(other_function)) {
    counter_other <- counter
  }
  
  for (cc in 1 : chains) {
    counter_chain <- GetER_1chain(dta, cutoffs = cutoffs[cc, , ],
                                  coefs = coefs[cc, , , ],
                                  predict_at = predict_at,
                                  mean_only = mean_only,
                                  other_function = other_function)
    if (mean_only) {
      counter[, cc, ] <- counter_chain$y
      if (!is.null(other_function)) {
        counter_other[, cc, ] <- counter_chain$y_other
      }
    } else {
      counter[, , cc, ] <- counter_chain$y
      if (!is.null(other_function)) {
        counter_other[, , cc, ] <- counter_chain$y_other
      }
    }
  }
  if (is.null(other_function)) {
    return(list(x = counter_chain$x, y = counter))
  }
  return(list(x = counter_chain$x, y = counter, y_other = counter_other))
}



