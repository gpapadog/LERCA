#' Marginal mean of C across experiments to have E(C|X =x) continuous in x.
#' 
#' @param meanCexp1 A vector of length equal to the number of covariates
#' including the marginal mean of each covariate in experiment 1.
#' @param XCcov A matrix of two dimensions. Element ij describes the covariance
#' of variable Ci with X in experiment j.
#' @param exper_change A vector including the points sbar that specify the
#' experiments. exper_change includes the minimum and maximum of the exposure
#' range.
#' @param meanX A vector including the mean of X in each experiment.
#' @param varX A vector including the variance of X in each experiment.
XCcontinuous <- function(meanCexp1, XCcov, exper_change, meanX, varX) {
  
  num_conf <- length(meanCexp1)
  num_exper <- dim(XCcov)[2]
  ret_means <- matrix(NA, nrow = num_conf, ncol = num_exper)
  ret_means[, 1] <- meanCexp1
  
  if (num_exper > 1) {
    for (cc in 1:num_conf) {
      for (ee in 2:num_exper) {
        
        term1 <- XCcov[cc, ee - 1] / varX[ee - 1]
        term1 <- term1 * (exper_change[ee] - meanX[ee - 1])
        
        term2 <- - XCcov[cc, ee] / varX[ee] * (exper_change[ee] - meanX[ee])
        ret_means[cc, ee] <- ret_means[cc, ee - 1] + term1 + term2
      }
    }
  }
  return(ret_means)
}