#' Plotting LERCA results.
#' 
#' Estimated ER, cutoffs, coefficients, inclusion probabilities.
#' 
#' @param dta The data that we are using with column X for exposure.
#' @param lerca The LERCA fit.
#' @param ER The estimated ER as acquired using GetER().
#' @param variable Which variable's results to be plotted. Intercept is 1,
#' Exposure is 2, and the covariates are 3 onwards.
#' @param wh_model Which model's coefficients and inclusion probabilities to
#' be plotted. Options are 1, 2 for exposure and outcome model accordingly.
#' @param probs The quantiles of the distribution to be plotted as intervals.
#' @param inclusion The inclusion probabilities as a function of the exposure
#' and as acquired using ExposureInclusion().
#' @param coefs The coefficients as a function of the exposure as acquired
#' using ExposureCoefs().
PlotLERCA <- function(dta, lerca, ER, variable = NULL, wh_model = NULL,
                      probs, inclusion = NULL, coefs = NULL) {
  
  dta <- as.data.frame(dta)
  
  exp_values <- ER$x
  ER_results <- t(apply(ER$y, 1, get_stats, probs = probs))
  ER_results <- as.data.frame(ER_results)
  ER_results$x <- exp_values
  
  # Plotting the ER.
  
  plot_ER <- ggplot() +
    geom_ribbon(data = ER_results, aes(x = x, ymin = LB, ymax = UB),
                fill = 'grey80') +
    geom_line(data = ER_results, aes(x = x, y = mean)) +
    theme(panel.background = element_blank()) +
    ylab('Response') + xlab('Exposure')
  
  
  # Plotting the cutoffs.
  
  cutoffs_results <- data.frame(x = as.numeric(lerca_short$cutoffs))
  plot_cutoffs <- ggplot(cutoffs_results, aes(x = x)) + 
    geom_histogram(aes(y = ..density..), binwidth = 0.06) +
    geom_density(alpha=.2, adjust = 1 / 2, fill="#FF6666") +
    xlab('') + ylab('Experiment Configuration') + xlim(range(dta$X)) +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  
  # Plotting the observed distribution.
  
  observed_exposure <- data.frame(x = dta$X)
  plot_exposures <- ggplot(observed_exposure, aes(x = x)) + 
    geom_histogram(aes(y = ..density..), binwidth = 0.1) +
    geom_density(alpha=.2, adjust = 1, fill="#FF6666") +
    xlab(expression(PM[2.5])) +
    ylab(expression(paste('Observed ', PM[2.5]))) +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  r <- list(plot_ER = plot_ER, plot_cutoffs = plot_cutoffs,
            plot_exposures = plot_exposures)
    

  # Plotting the coefficient.

  if (!is.null(coefs) & !is.null(variable)) {
    
    slope_stat <- apply(coefs[, , wh_model, , variable], 3, get_stats,
                        probs = probs)
    slope_stat <- as.data.frame(t(slope_stat))
    slope_stat$x <- exp_values
    
    plot_coef <- ggplot() +
      geom_ribbon(data = slope_stat, aes(x = x, ymin = LB, ymax = UB),
                  fill = 'grey80') +
      geom_line(data = slope_stat, aes(x = x, y = mean)) +
      theme(panel.background = element_blank()) +
      geom_hline(yintercept = 0, linetype = 2) +
      ylab('Coefficient') + xlab('')
  
    r <- append(r, list(plot_coef = plot_coef))
  }
  
  # Plotting the inclusion indicators.
  
  if (!is.null(inclusion) & !is.null(variable) & (variable > 2)) {
    
    inclusion_plot <- apply(inclusion[, , , variable - 2], 2 : 3, mean)
    inclusion_plot <- adply(inclusion_plot, 1 : 2)
    names(inclusion_plot)[3] <- 'value'
    inclusion_plot$exposure <- as.numeric(as.character(inclusion_plot$exposure))

    plot_inclusion <- ggplot(data = inclusion_plot,
                             aes(x = exposure, y = value, linetype = model)) +
      geom_line() + xlab('Exposure') + ylab('Inclusion probability')
    
    r <- append(r, list(plot_inclusion = plot_inclusion))
  }
  return(r)
}


