#' Create arrays to save MCMC posterior samples
#' 
#' @param X Numeric vector. Observed exposure values. Can be left NULL if
#' min_exper_sample is set to 0.
#' @param chains The number of separate MCMC chains.
#' @param Nsims The number of posterior samples per chain.
#' @param num_exper The number of experiments we are allowing.
#' @param num_conf The number of potential confounders.
#' @param omega The omega of the BAC prior on inclusion indicators.
#' @param minX The minimum observed exposure value.
#' @param maxX The maximum observed exposure value.
#' @param starting cutoffs Matrix with rows corresponding to different chains.
#' Each row includes K ordered values of MCMC starting cutoffs. If left NULL,
#' random started values are used.
#' @param starting_alphas Array with dimensions corresponding to the model
#' (exposure / outcome), the experiment, and the potential confounders. Entries
#' 0/1 represent exclusion/inclusion of the covariate in the corresponding
#' model.
#' @param min_exper_sample The minimum number of observations within an
#' experiment. It will be used to ensure that starting cutoffs are allowed
#' under the prior specification.
#' 
MakeArrays <- function(X = NULL, chains, Nsims, num_exper, num_conf, omega,
                       minX, maxX, starting_cutoffs, starting_alphas,
                       min_exper_sample = 0) {
  
  # Alphas. Starting values are from the prior.
  
  alphas <- NULL
  if (num_conf > 0) {
    alphas <- array(NA, dim = c(2, chains, Nsims, num_exper, num_conf))
    dimnames(alphas) <- list(model = c('Exposure', 'Outcome'),
                             chain = 1 : chains, sample = 1 : Nsims,
                             exper = 1 : num_exper, conf = 1 : num_conf)
    
    if (is.null(starting_alphas)) {
      starting_alphas <- array(NA, dim = c(2, chains, num_exper, num_conf))
      for (cc in 1 : chains) {
        starting_alphas[2, cc, , ] <- sample(c(0, 1), num_exper * num_conf,
                                             prob = c(omega + 1, 2 * omega),
                                             replace = TRUE)
        for (ee in 1 : (K + 1)) {
          for (jj in 1 : num_conf) {
            probs <- ifelse(starting_alphas[2, cc, ee, jj] == 0, omega, 1)
            starting_alphas[1, cc, ee, jj] <- sample(c(0, 1), 1,
                                                     prob = c(probs, 1))
          }
        }
      }
    }
    alphas[, , 1, , ] <- starting_alphas
  }
  
  # Coefficient and variance values.
  
  variances <- array(1, dim = c(2, chains, Nsims, num_exper))
  dimnames(variances) <- list(model = c('Exposure', 'Outcome'),
                              chain = 1 : chains, sample = 1 : Nsims,
                              exper = 1 : num_exper)
  
  cov_names <- c('Int', 'X')
  if (num_conf > 0) {
    cov_names <- c(cov_names, paste0('C', 1 : num_conf))
  }
  coefs <- array(0, dim = c(2, chains, Nsims, num_exper, num_conf + 2))
  dimnames(coefs) <- list(model = c('Exposure', 'Outcome'),
                          chain = 1 : chains, sample = 1 : Nsims,
                          exper = 1 : num_exper, covar = cov_names)
  coefs[1, , , , 2] <- NA  # No exposure coefficient for exposure model.
  
  
  # Experiment configuration with starting values from the prior if NULL.
  
  cutoffs <- array(0, dim = c(chains, Nsims, K))
  dimnames(cutoffs) <- list(chain = 1 : chains, sample = 1 : Nsims,
                            point = 1 : K)
  
  if (is.null(starting_cutoffs)) {
    starting_cutoffs <- matrix(NA, nrow = chains, ncol = K)
    
    for (cc in 1 : chains) {
      
      # Ensuring that the starting cutoffs are allowed by the prior.
      sample_new <- TRUE
      while (sample_new) {
        # Get the even ordered statistics from uniform sample.
        x <- sort(runif(2 * (K + 1), minX, maxX))
        start_cuts <- x[seq(2, 2 * K, by = 2)]
        # If we have no constraint on the minimum sample size, accept.
        if (min_exper_sample == 0) {
          sample_new <- FALSE
        } else {  # Check that minimum sample size is satisfied.
          cuts_x <- c(minX - 0.001, start_cuts, maxX + 0.001)
          expers <- table(sapply(X, function(x) sum(x < cuts_x)))
          if (length(expers) == K + 1 & all(expers >= min_exper_sample)) {
            sample_new <- FALSE
          }
        }
      }
      starting_cutoffs[cc, ] <- start_cuts
    }
  }
  cutoffs[, 1, ] <- starting_cutoffs[1 : chains, ]
  
  
  
  # Starting values for coefficients and variances.
  for (cc in 1 : chains) {
    for (ee in 1 : (K + 1)) {
      coefs[1, cc, 1, ee, - 2] <- rnorm(num_conf + 1, mean = 0, sd = 10)
      coefs[2, cc, 1, ee, - 1] <- rnorm(num_conf + 1, mean = 0, sd = 10)
    }
    variances[1, cc, 1, ] <- invgamma::rinvgamma(num_exper, 10, 20)
    variances[2, cc, 1, ] <- invgamma::rinvgamma(num_exper, 10, 20)
  }
  # Intercept starting values.
  for (cc in 1 : chains) {
    exact_cuts <- c(minX, cutoffs[cc, 1, ], maxX)
    coefs[2, cc, 1, ee, 1] <- rnorm(1, mean = 0, sd = 10)
    
    for (ee in 1 : K) {
      interval <- exact_cuts[ee + 1] - exact_cuts[ee]
      coefs[2, cc, 1, ee + 1, 1] <- coefs[2, cc, 1, ee, 1] +
        coefs[2, cc, 1, ee, 2] * interval
    }
  }
  
  return(list(alphas = alphas, cutoffs = cutoffs, coefs = coefs,
              variances = variances))
}