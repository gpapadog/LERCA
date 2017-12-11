#' Simulating data with differential confounding.
#' 
#' Funcition that generates data where different variables are confounders of
#' the exposure-response relationship based on the exposure range we are
#' looking at.
#' 
#' @param N Sample size.
#' @param num_exper Number of experiments.
#' @param XCcorr A matrix of 2 dimensions. Element ij corresponds to the
#' correlation of variable Ci with exposure X in experiment j.
#' @param varC A vector of length equal to the number of covariates, describing
#' the variance of each of the C. We assume constant variance of C across
#' experiments.
#' @param Xrange Vector of length 2. These values correspond to the minimum and
#' the maximum of the uniform distribution from which X is generated.
#' @param exper_range The points of the exposure range Xrange at which
#' differential confounding might occur.
#' @param meanCexp1 A vector of length equal to the number of covariates
#' including the marginal mean of each covariate in experiment 1. The mean of 
#' the covariates in the rest experiments is specified by the algorithm such
#' that E(C|X=x) is continuous in x. If equal to NULL, 0 will be specified for
#' all variables.
#' @param out_coef A matrix of two dimensions. Element ij describes the
#' coefficient of Ci in the outcome model of experiment j.
#' @param interYexp1 The intercept of the outcome model in experiment 1. The
#' remaining intercepts are set such that the function E[Y | X = x] is
#' continuous in x.
#' @param bYX If XY_function is set to linear, bYX includes the coefficient of
#' the exposure within each experiment and is of length num_exper. If
#' XY_function is set to other, bYX is of length 1 and corresponds to the
#' coefficient in front of the exposure term.
#' @param Ysd A vector including the residual variance of the outcome model in
#' each experiment. Defaults to 1.
#' @param overall_meanC A string equal to 'true' or 'observed'. Defaults to
#' 'true'. This argument controls whether we will use the analytical or the
#' observed mean of C to ensure continuous ER at the points of the experiment
#' configuration. If set to 'true', individual simulated data sets do not
#' necessarily have a continuous ER, but they do in general. If set equal to
#' the 'observed', each simulated data set has a continuous ER, that might be
#' slightly different for every dataset. Both choices average to the same true
#' ER, for a large sample size, or number of replicated data sets.
#' @param XY_function String specifying whether the XY relationship is piece-
#' wise linear (set 'linear'), or a continuous function supplied by the XY_spec
#' arguement (set 'other').
#' @param XY_spec Needs to be specified if XY_function is set to 'other'. It is
#' the function that specifies the true ER relationship. Defaults to NULL.
#' 
#' @export
#' @examples
#' XCcorr <- matrix(c(0.4, 0.2, 0.2, 0), 2, 2)
#' out_coef <- matrix(c(0.2, 0.3, 0, 0.4), 2, 2)
#' bYX <- c(0, 0, 1, 2, 3)
#' sim <- SimDifferentialConfounding(N = 1000, num_exper = 2, XCcorr = XCcorr,
#'                                   varC = c(1, 1), Xrange = c(0, 10),
#'                                   exper_change = c(0, 4, 10),
#'                                   out_coef = out_coef, bYX = bYX)
#' f <- function(x) {
#'   return(x ^ 2)
#' }
#' sim <- SimDifferentialConfounding(N = 1000, num_exper = 2, XCcorr = XCcorr,
#'                                   varC = c(1, 1), Xrange = c(0, 10),
#'                                   exper_change = c(0, 4, 10),
#'                                   out_coef = out_coef, bYX = bYX,
#'                                   XY_function = 'other', XY_spec = f)

SimDifferentialConfounding <- function(N, num_exper, XCcorr, varC, Xrange,
                                       exper_change = NULL, meanCexp1 = NULL,
                                       out_coef, interYexp1 = 0, bYX, Ysd = 1,
                                       overall_meanC = c('true', 'observed'),
                                       XY_function = c('linear', 'other'),
                                       XY_spec = NULL) {
  
  overall_meanC <- match.arg(overall_meanC)
  XY_function <- match.arg(XY_function)
  
  num_conf <- nrow(XCcorr)
  if (length(Ysd) == 1) {
    Ysd <- ScalarToVector(Ysd, num_exper)
  }
  
  if (is.null(meanCexp1)) {
    meanCexp1 <- rep(0, num_conf)
  }
  if (is.null(exper_change)) {
    exper_change <- seq(Xrange[1], Xrange[2], length.out = num_exper + 1)
  }
  meanC_weights <- exper_change[- 1] - exper_change[- length(exper_change)]
  
  SimDiffConfChecks(num_exper, XCcorr, exper_change, Xrange, varC, meanCexp1)
  
  print('X is generated from a uniform distribution.')
  # If X is not from uniform, adjust GetbYvalues to correctly calculate the
  # mean of C. Also adjust meanX, varX.
  
  X <- runif(N, min = Xrange[1], max = Xrange[2])
  E <- as.numeric(cut(X, exper_change, include.lowest = TRUE, right = TRUE))
  dta <- data.table::data.table(X = X, E = E)
  
  print('Data per experiment:')
  print(with(dta, table(E)))
  
  # In order to generate C, we will create the "effective experiment" of each
  # observation. Even if exper_change defines a point s, if the XCcorr values
  # do not change, the effective experiment for the exposure is the same.
  new_exp_exper <- sapply(2 : num_exper, function(x) any(XCcorr[, x] !=
                                                           XCcorr[, x - 1]))
  new_exp_exper <- as.numeric(c(TRUE, new_exp_exper))
  dta$E_eff_exp <- sapply(dta$E, function(x) sum(new_exp_exper[1 : x]))
  
  # Then, we will have the effective experiment change representing the points
  # where the exposure-covariate correlation actually changes.
  eff_exper_change <- exper_change[new_exp_exper == 1]
  eff_num_exper <- length(eff_exper_change) - 1
  eff_XCcorr <- XCcorr[, new_exp_exper == 1, drop = FALSE]
  
  # Consider if I want to make this equal to the observed.
  meanX <- sapply(2 : length(eff_exper_change), function(x)
    mean(eff_exper_change[c(x - 1, x)]))
  varX <- sapply(2 : length(eff_exper_change), function(x)
    (eff_exper_change[x] - eff_exper_change[x - 1]) ^ 2 / 12)
  XCcov <- sapply(1 : eff_num_exper, function(x) eff_XCcorr[, x] *
                    sqrt(varC * varX[x]))
  
  meanC <- XCcontinuous(meanCexp1, XCcov, eff_exper_change, meanX, varX)
  
  Cfull <- NULL
  data.table::setkey(dta, E)
  
  for (ee in 1 : eff_num_exper) {
    
    wh <- with(dta, which(E_eff_exp == ee))
    X_meanX <- matrix(dta$X[wh] - meanX[ee], ncol = 1)
    covXC_ee <- XCcov[, ee, drop = FALSE]
    
    mu_bar <- 1 / varX[ee] * covXC_ee %*% t(X_meanX)
    mu_bar <- sweep(mu_bar, 1, meanC[, ee], FUN = '+')
    sigma_bar <- diag(varC) - 1 / varX[ee] * covXC_ee %*% t(covXC_ee)
    
    C <- matrix(nrow = length(wh), ncol = num_conf)
    for (ii in 1:length(wh)) {
      C[ii, ] <- mvnfast::rmvn(1, mu = mu_bar[, ii], sigma = sigma_bar)
    }
    Cfull <- rbind(Cfull, C)
  }
  
  dta <- dta[, c('X', 'E')]
  dta <- cbind(dta, Cfull)
  data.table::setnames(dta, names(dta)[- c(1, 2)], paste0('C', 1 : num_conf))
  
  # Whether we want to use the true overall mean or the observed.
  over_meanC <- meanC
  if (eff_num_exper < num_exper) {
    if (overall_meanC == 'true') {
      warning('overall_meanC set to observed when XCcorr has consecutive columns.')
    }
    overall_meanC <- 'observed'
  }
  
  if (overall_meanC == 'observed') {
    over_meanC <- colMeans(dta[, 3 : (2 + num_conf), with = FALSE])
    over_meanC <- matrix(over_meanC, nrow = num_conf, ncol = num_exper)
  }
  
  # Generating the outcome 'Y'.
  bY <- GetbYvalues(num_exper, exper_change = exper_change, bYX = bYX,
                    interYexp1 = interYexp1, out_coef = out_coef,
                    meanC = over_meanC, weights = meanC_weights,
                    XY_function = XY_function)
  dta$Y <- GenYgivenXC(dta, out_coef, bY, bYX, Ysd, XY_function, XY_spec)
  
  # Calculating covariance and correlation in simulated data.
  cov_dta <- data.frame(Y = dta$Y, X = dta$X, E = dta$E)
  cov_dta <- cbind(cov_dta, Cfull)
  
  covariances <- MakeCovArray(num_conf = num_conf, num_exper = num_exper)
  correlations <- covariances
  
  for (ee in 1:num_exper) {
    covariances[, , ee] <- cov(cov_dta[cov_dta$E == ee, - 3])
    correlations[, , ee] <- cor(cov_dta[cov_dta$E == ee, - 3])
  }
  
  return(list(data = dta, covar = covariances, corr = correlations, bY = bY,
              meanC = meanC, XCcov = XCcov, meanX = meanX, varX = varX))
}



SimDiffConfChecks <- function(num_exper, XCcorr, exper_change, Xrange,
                              varC, meanCexp1) {
  if (num_exper != ncol(XCcorr)) {
    stop('Number of columns in XCcor should be equal to num_exper.')
  }
  if (any(exper_change < Xrange[1]) | any(exper_change > Xrange[2])) {
    stop('exper_change should be inside Xrange.')
  }
  if (length(varC) != nrow(XCcorr)) {
    step('VarC is a vector of length equal to the number of confounders.')
  }
  if (length(meanCexp1) != nrow(XCcorr)) {
    stop('meanCexp1 should be a vector of length number of covariates.')
  }
}

MakeCovArray <- function(num_conf, num_exper) {
  covariances <- array(NA, dim = c(num_conf + 2, num_conf + 2, num_exper))
  dimnames(covariances)[[1]] <- c('Y', 'X', paste0('C', 1:num_conf))
  dimnames(covariances)[[2]] <- dimnames(covariances)[[1]]
  dimnames(covariances)[3] <- list(exper = 1:num_exper)
  return(covariances)
}
