#' Posterior weights of all models in MCMC.
#' 
#' Function that calculates and orders the models that were chosen in the MCMC
#' at least once.
#' 
#' @param mcmc The LERCA fit.
#' @param model Options are 'Outcome' or 'Exposure'. Defines which model we 
#' want to look at.
#' @param experiment Which experiment to look at. Defaults to 1. Cannot be more
#' than K + 1. Results will be interpretable for the exposure range in the
#' intersection of experiments across MCMC iterations.
#' 
#' @export
ModelWeights <- function(mcmc, model = c('Outcome', 'Exposure'),
                         experiment = 1) {
  
  num_cov <- dim(mcmc$alphas)[5]
  model <- match.arg(model)
  model <- ifelse(model == 'Outcome', 2, 1)
  if (experiment > dim(mcmc$alphas)[4]) {
    stop('Specify experiment within the experiment range.')
  }
  
  alphas <- mcmc$alphas[model, , , experiment, ]
  alphas <- plyr::adply(alphas, 1)
  alphas <- alphas[, - 1]
  names(alphas) <- paste0('C', 1 : num_cov)
  alphas$entry <- 1 : nrow(alphas)
  
  alphas <- data.table::as.data.table(alphas)
  unique_alpha <- alphas[, list(number_times = length(entry)),
                         by = eval(names(alphas)[1 : num_cov])]
  unique_alpha[, proportion := number_times / nrow(alphas)]
  data.table::setorder(unique_alpha, - number_times)
  
  return(unique_alpha)
}