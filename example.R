library(LERCA)
library(ggplot2)
data('toyData')

K <- 2
omega <- 5000

chains <- 2
Nsims <- 6000
burn <- 4000
thin <- 4

set.seed(1234)
lerca <- LERCA(dta = toyData, chains = chains, Nsims = Nsims, K = K,
               cov_cols = 3 : 6, omega = omega)

lerca_short <- BurnThin(lerca = lerca, burn = burn, thin = thin)

waic <- WAIC(lerca = lerca_short, dta = toyData)
waic

# Get the ER estimates over a set of exposure values.
ER <- GetER(dta = toyData, cutoffs = lerca_short$cutoffs,
            coefs = lerca_short$coefs[2, , , , ], mean_only = TRUE)

# Acquire the inclusion probabilities as a function of the exposure for both models.
inclusion <- ExposureInclusion(lerca = lerca_short, exp_values = ER$x)

# Acquire the coefficients as a function of the exposure.
coefs <- ExposureCoefs(lerca = lerca_short, exp_values = ER$x)

# Acquire models' posterior weights.
post_weights <- ModelWeights(lerca_short, model = 'Outcome', experiment = 1)
post_weights[1 : 4, ]

probs <- c(0.1, 0.9)
plots <- PlotLERCA(dta = toyData, lerca = lerca_short, ER = ER, probs = probs,
                   coefs = coefs, variable = 2, wh_model = 2)
plots$plot_ER
plots$plot_cutoffs
plots$plot_exposures
plots$plot_coef

# Plot the inclusion probability of C2 as a function of the exposure:
PlotLERCA(dta = toyData, lerca = lerca_short, ER = ER, probs = probs,
          coefs = coefs, inclusion = inclusion, variable = 2 + 3,
          wh_model = 2)$plot_inclusion
