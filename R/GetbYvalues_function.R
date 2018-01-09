#' Get the bY values that make the ER curve continuous in simulated data when
#' the X-Y relationship is truly piecewise linear.
#'
#' @param num_exper The number of experiments. Should be > 1.
#' @param Xrange The range of the exposure range. Specify if exper_change is
#' NULL.
#' @param exper_change The points of experiment change. Should be including the
#' edges of the exposure range. If left NULL, exper_change will be specified as
#' equally distanced points on Xrange such that num_exper experiments are
#' defined.
#' @param bYX If XY_function is set to linear, bYX includes the coefficient of
#' the exposure within each experiment and is of length num_exper. If
#' XY_function is set to other, bYX is of length 1 and corresponds to the
#' coefficient in front of the exposure term.
#' @param interYexp1 The intercept of the outcome model in experiment 1. If not
#' specified, it defaults to 0.
#' @param out_coef A matrix with number of rows equal to the number of
#' covariates and number of columns equal to the number of experiments
#' including the outcome model coefficients.
#' @param meanC A matrix of equal dimensions to out_coef including the mean of
#' each covariate in each experiment.
#' @param weights The weights to be given to the mean of each experiment in the
#' calculation of the overall mean. Vector of length equal to the number of
#' experiments. If left NULL, each experiment is given equal weight.
#' @param XY_function String specifying whether the XY relationship is piece-
#' wise linear (set 'linear'), or a continuous function supplied by the XY_spec
#' arguement (set 'other').
#'
#' @export
GetbYvalues <- function(num_exper, Xrange = NULL, exper_change = NULL, bYX,
                        interYexp1 = 0, out_coef, meanC, weights = NULL,
                        XY_function = c('linear', 'other')) {
  
  num_conf <- nrow(out_coef)
  XY_function <- match.arg(XY_function)
  
  if (num_exper == 1) {
    return(interYexp1)
  }
  if (num_exper < 1) {
    stop('Set num_exper equal to or greater than 1.')
  }
  if (ncol(out_coef) != num_exper) {
    stop('out_coef not compatible with num_exper.')
  }
  if (is.null(exper_change)) {
    if (is.null(Xrange)) {
      stop('Xrange has to be specified if exper_change is NULL.')
    }
    exper_change <- seq(Xrange[1], Xrange[2], length.out = num_exper + 1)
  }
  
  bY <- c(interYexp1, rep(NA, num_exper - 1))
  
  if (is.null(weights)) {
    weights <- rep(1, num_exper)
  }
  weights <- weights / sum(weights)
  meanC <- sweep(meanC, 2, weights, FUN = '*')
  overall_meanC <- rowSums(meanC)

  # Calculating the intercepts.
  for (ii in 2:num_exper) {
    bY[ii] <- bY[ii - 1]
    coef_diff <- out_coef[, ii - 1] - out_coef[, ii]
    bY[ii] <- bY[ii] + sum(coef_diff * overall_meanC)
    if (XY_function == 'linear') {
      bY[ii] <- bY[ii] + (bYX[ii - 1] - bYX[ii]) * exper_change[ii]
    }
  }
  
  return(bY)
}
