MakeArrays <- function(chains, Nsims, num_exper, num_conf, omega, minX, maxX,
                       starting_cutoffs, starting_alphas) {
  
  # Alphas. Starting values are from the prior.
  
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
  
  # Coefficient and variance values.
  
  variances <- array(1, dim = c(2, chains, Nsims, num_exper))
  dimnames(variances) <- list(model = c('Exposure', 'Outcome'),
                              chain = 1 : chains, sample = 1 : Nsims,
                              exper = 1 : num_exper)
  
  coefs <- array(0, dim = c(2, chains, Nsims, num_exper, num_conf + 2))
  dimnames(coefs) <- list(model = c('Exposure', 'Outcome'),
                          chain = 1 : chains, sample = 1 : Nsims,
                          exper = 1 : num_exper,
                          covar = c('Int', 'X-s', 1 : num_conf))
  
  # Starting values for coefficients and variances.
  for (cc in 1 : chains) {
    for (ee in 1 : (K + 1)) {
      coefs[1, cc, 1, ee, - 2] <- rnorm(num_conf + 1, mean = 0, sd = 10)
      coefs[2, cc, 1, ee, ] <- rnorm(num_conf + 2, mean = 0, sd = 10)
    }
    variances[1, cc, 1, ] <- invgamma::rinvgamma(num_exper, 10, 20)
    variances[2, cc, 1, ] <- invgamma::rinvgamma(num_exper, 10, 20)
  }
  
  
  # Experiment configuration with starting values from the prior if NULL.
  
  cutoffs <- array(0, dim = c(chains, Nsims, K))
  dimnames(cutoffs) <- list(chain = 1 : chains, sample = 1 : Nsims,
                            point = 1 : K)
  
  if (is.null(starting_cutoffs)) {
    starting_cutoffs <- matrix(NA, nrow = chains, ncol = K)
    for (cc in 1 : chains) {
      x <- sort(runif(2 * (K + 1), minX, maxX))
      starting_cutoffs[cc, ] <- x[seq(2, 2 * K, by = 2)]
    }
  }
  cutoffs[, 1, ] <- starting_cutoffs[1 : chains, ]

  return(list(alphas = alphas, cutoffs = cutoffs, coefs = coefs,
              variances = variances))
}