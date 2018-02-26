#' Calculate the CER
#' 
#' Calculate the average ER function from posterior cutoff and coefficient values from
#' one MCMC chain.
#' 
#' Function that takes the posterior samples of cutoffs, inclusion indicators and
#' coefficients of the outcome model and returns the posterior mean response for
#' multiple values of the exposure.
#' 
#' @param dta The dataset including columns X (treatment), Y (outcome), and potential
#' confounders named as C1, C2, ...
#' @param cutoffs A list of length equal to the number of posterior samples. Each
#' element of the list is a vector including the current set of cutoff values.
#' @param coefs A list of length equal to the number of posterior samples. Each
#' element of the list is a matrix of the outcome coefficients. The dimension of each
#' matrix is (number of experiments at current iteration) x (number of potential
#' confounders). Confounders with inclusion indicators equal to 0, have coefficients
#' equal to 0. The first two elements of each row of the matrix correspond to intercept
#' and slope of the exposure.
#' @param predict_at A vector of the exposure values we wish to predict ER at. If left
#' NULL, a grid of length grid_length over the exposure values will be created.
#' @param grid_length The number of locations we will predict the ER. Defaults to 100.
#' @param mean_only Logical. Set to FALSE if we want the individual ER predictions. Set
#' to TRUE if we are only interested in the posterior samples of the mean ER. Defaults
#' to FALSE.
#' @param other_function Whether we want to apply a different function to the
#' predictions. Defaults to NULL. Examples include exp(x), log(x).
GetER_1chain <- function(dta, cutoffs, coefs, predict_at = NULL, grid_length = 100,
                         mean_only = FALSE, other_function = NULL) {
  
  dta <- as.data.frame(dta)
  minX <- min(dta$X)
  maxX <- max(dta$X)
  Nsims <- dim(coefs)[1]
  
  if (is.null(predict_at)) {
    predict_at <- seq(minX, maxX, length.out = grid_length)
  }
  
  num_conf <- dim(coefs)[3] - 2
  des_mat <- cbind(Int = 1, X = dta$X)
  if (num_conf > 0) {
    cov_cols <- which(names(dta) %in% paste0('C', 1 : num_conf))
    des_mat <- cbind(des_mat, as.matrix(dta[, cov_cols]))
    meanC <- colMeans(dta[, cov_cols])
  }

  # Array where the predicted counterfactual values are saved.
  if (mean_only) {
    counter <- array(NA, dim = c(length(predict_at), Nsims))
    dimnames(counter) <- list(predict = predict_at, sim = 1:Nsims)
  } else {
    counter <- array(NA, dim = c(length(predict_at), nrow(dta), Nsims))
    dimnames(counter) <- list(predict = predict_at, obs = 1:nrow(dta), sim = 1:Nsims)
  }
  if (!is.null(other_function)) {
    counter_other <- counter
  }
  
  predict_in_range <- which(predict_at >= minX & predict_at <= maxX)
  
  
  # If we only want ONLY the mean of the linear function, there is simple way.
  if (mean_only & is.null(other_function)) {

    for (ii in 1 : Nsims) {
      current_cutoffs <- cutoffs[ii, ]
      exact_cuts <- c(minX, current_cutoffs, maxX)
      current_betas <- coefs[ii, , ]
      for (pp in 1 : length(predict_in_range)) {
        curr_loc <- predict_at[predict_in_range[pp]]
        exper <- sum(current_cutoffs <= curr_loc) + 1
        curr_pred <- c(1, curr_loc - exact_cuts[exper])
        if (num_conf > 0) {
          curr_pred <- c(curr_pred, meanC)
        }
        counter[predict_in_range[pp], ii] <- sum(curr_pred *
                                                   current_betas[exper, ])
      }
    }
    return(list(x = predict_at, y = counter))
  }
  
  
  # If we want posterior samples for every observation, or we want a different
  # function, we must run the following.
  
  for (ii in 1 : Nsims) {

    current_cutoffs <- cutoffs[ii, ]
    exact_cuts <- c(minX, current_cutoffs, maxX)
    current_betas <- coefs[ii, , ]
    k <- length(current_cutoffs)
    
    for (pp in 1 : length(predict_in_range)) {
      
      # Find the experiment that the predicting exposure is in.
      exper <- sum(current_cutoffs <= predict_at[predict_in_range[pp]]) + 1
      
      D <- des_mat
      D[, 2] <- predict_at[predict_in_range[pp]] - exact_cuts[exper]
      predictions <- D %*% current_betas[exper, ]
      
      if (mean_only & !is.null(other_function)) {
        counter[predict_in_range[pp], ii] <- mean(predictions)
        counter_other[predict_in_range[pp], ii] <- mean(other_function(predictions))
        res <- list(x = predict_at, y = counter, y_other = counter_other)
        
      } else if (!mean_only & is.null(other_function)) {
        counter[predict_in_range[pp], , ii] <- predictions
        res <- list(x = predict_at, y = counter)
        
      } else if (!mean_only & !is.null(other_function)) {
        counter[predict_in_range[pp], , ii] <- predictions
        counter_other[predict_in_range[pp], , ii] <- other_function(predictions)
        res <- list(x = predict_at, y = counter, y_other = counter_other)
        
      } else {
        stop('Error. Should never happen.')
      }
    }
  }
  return(res)
}


