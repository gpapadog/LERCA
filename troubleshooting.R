# Simulation of fixed number of points in the experiment configuration.

library(data.table)
library(invgamma)
library(truncnorm)
library(mvnfast)

setwd('~/Documents/Causal_ER_curve/LERCA/')
load_path <- 'Data_specs/data_specs4.Rdata'

source_path <- '~/Github/LERCA/R/'
file_source <- c('SimDifferentialConfounding',
                 'XCcontinuous', 'GetbYvalues', 'GenYgivenXC',
                 'WAIC', 'GetER', 'GetER_1chain')
sapply(paste0(source_path, file_source, '_function.R'), source, .GlobalEnv)

source_path <- '~/Github/LERCA_FullCond/'
file_source <- list.files(source_path, pattern = '*_function.R$')
sapply(paste0(source_path, file_source), source, .GlobalEnv)



index <- 3

load(load_path)
attach(data_specs)

# --------

N <- 200 * num_exper
overall_meanC <- 'observed'

chains <- 3
Nsims <- 1000
burn <- 30000
thin <- 70
coverged_psr_diff <- 0.05
plot_every <- 500

prop_distribution <- 'Uniform'
# normal_percent <- 1 / 2
BIC_approximation <- TRUE
omega <- 5000
comb_probs <- c(0.01, 0.5, 0.99)
split_probs <- c(0.2, 0.95)
s_upd_probs <- c(99 / 100, 1 / 100)
K <- num_exper - 1


print(index)
set.seed(index * 2)

sim <- SimDifferentialConfounding(N = N, num_exper = num_exper,
                                  XCcorr = XCcorr, varC = varC,
                                  Xrange = Xrange, bYX = bYX,
                                  exper_change = exper_change,
                                  meanCexp1 = meanCexp1, out_coef = out_coef,
                                  interYexp1 = interYexp1, Ysd = Ysd,
                                  overall_meanC = overall_meanC,
                                  XY_function = XY_function,
                                  XY_spec = XY_spec)
dta <- as.data.frame(sim$data)
cov_cols <- which(names(dta) %in% paste0('C', 1 : num_conf))


# ------- STEP 1. Priors -------- #

# Priors on the regression coefficients.
mu_priorX <- rep(0, num_conf + 1)
mu_priorY <- rep(0, num_conf + 2)
Sigma_priorX <- diag(rep(100 ^ 2, num_conf + 1))
Sigma_priorY <- diag(rep(100 ^ 2, num_conf + 2))

# Priors on the variance components.
alpha_priorX <- 0.001
alpha_priorY <- 0.001
beta_priorX <- 0.001
beta_priorY <- 0.001

# ------ LERCA code ------ #

num_exper <- K + 1
num_conf <- length(cov_cols)
minX <- min(dta$X)
maxX <- max(dta$X)


# -------- Where to save --------- #
arrays <- MakeArrays(chains = chains, Nsims = Nsims, num_exper = num_exper,
                     num_conf = num_conf, omega = omega, minX = minX,
                     maxX = maxX, starting_cutoffs = NULL)
alphas <- arrays$alphas
cutoffs <- arrays$cutoffs
coefs <- arrays$coefs
variances <- arrays$variances

acc <- array(0, dim = c(2, 2, chains))
dimnames(acc) <- list(kind = c('separate', 'jumpOver'),
                      num = c('attempt', 'success'),
                      chain = 1 : chains)


# -------- STEP 3. MCMC. ---------- #

for (cc in 1 : chains) {
  for (ii in 2 : Nsims) {
    
    if (ii %% 100 == 0) {
      print(paste('Iteration', ii))
    }
    
    # ----- Update experiment configuration.
    current_cutoffs <- cutoffs[cc, ii - 1, ]
    current_coefs <- coefs[, cc, ii - 1, , ]
    current_vars <- variances[, cc, ii - 1, ]
    current_alphas <- alphas[, cc, ii - 1, , ]

    wh_s_upd <- sample(c(1, 2), 1, prob = s_upd_probs)
    if (wh_s_upd == 1) {

      acc[1, 1, cc] <- acc[1, 1, cc] + 1
      
      exp_upd <- UpdateExperiments(dta = dta, cov_cols = cov_cols,
                                   current_cutoffs = current_cutoffs,
                                   current_coefs = current_coefs,
                                   current_vars = current_vars,
                                   prop_distribution = prop_distribution,
                                   normal_percent = normal_percent)
      cutoffs[cc, ii, ] <- exp_upd$cutoffs
      acc[1, 2, cc] <- acc[1, 2, cc] + exp_upd$acc
      
      # ---- Update alphas.
      current_cutoffs <- cutoffs[cc, ii, ]
      current_alphaY <- alphas[2, cc, ii - 1, , ]
      current_coefs <- coefs[, cc, ii - 1, , ]
      current_vars <- variances[, cc, ii - 1, ]
      
      alphas_upd <- UpdateAlphas(dta = dta, cov_cols = cov_cols,
                                 current_cutoffs = current_cutoffs,
                                 current_alphaY = current_alphaY,
                                 current_coefs = current_coefs,
                                 current_vars = current_vars,
                                 Sigma_priorX = Sigma_priorX,
                                 mu_priorX = mu_priorX,
                                 Sigma_priorY = Sigma_priorY,
                                 mu_priorY = mu_priorY, omega = omega)
      alphas[, cc, ii, , ] <- alphas_upd      

    } else if (wh_s_upd == 2) {
      
      acc[2, 1, cc] <- acc[2, 1, cc] + 1
      
      jump_upd <- JumpOver(dta = dta, current_cutoffs = current_cutoffs,
                           current_alphas, approximate = TRUE,
                           cov_cols = cov_cols, omega = omega,
                           comb_probs = comb_probs, split_probs = split_probs)
      
      acc[2, 2, cc] <- acc[2, 2, cc] + jump_upd$acc
      cutoffs[cc, ii, ] <- jump_upd$new_cutoffs
      alphas[, cc, ii, , ] <- jump_upd$new_alphas
    }
    
    
    # ----- Updating the coefficients.
    current_cutoffs <- cutoffs[cc, ii, ]
    current_alphas <- alphas[, cc, ii, , ]
    current_vars <- variances[, cc, ii - 1, ]
    
    coef_upd <- UpdateCoefficients(dta = dta, cov_cols = cov_cols,
                                   current_cutoffs = current_cutoffs,
                                   current_alphas = current_alphas,
                                   current_vars = current_vars,
                                   Sigma_priorX = Sigma_priorX,
                                   mu_priorX = mu_priorX,
                                   Sigma_priorY = Sigma_priorY,
                                   mu_priorY = mu_priorY)
    coefs[, cc, ii, , ] <- coef_upd
    
    # ------ Updating the variances.
    
    current_coefs <- coefs[, cc, ii, , ]
    var_upd <- UpdateVariances(dta = dta, current_cutoffs = current_cutoffs,
                               current_coefs = current_coefs,
                               alpha_priorX = alpha_priorX,
                               beta_priorX = beta_priorX,
                               alpha_priorY = alpha_priorY,
                               beta_priorY = beta_priorY)
    variances[, cc, ii, ] <- var_upd
    
    
    if (plot_every > 0) {
      if (ii %% plot_every == 0) {
        par(mfrow = c(2, ceiling(K / 2)), mar = rep(2, 4))
        
        for (kk in 1 : K) {
          plot(cutoffs[cc, 1 : ii, kk], type = 'l', col = cc)
        }
        
      }
    }
  }
}

# ----- END OF LERCA CODE ----- #

lerca <- list(cutoffs = cutoffs[, 1 : ii, , drop = FALSE],
              alphas = alphas[, , 1 : ii, , , drop = FALSE],
              coefs = coefs[, , 1 : ii, , , drop = FALSE],
              variances = variances[, , 1 : ii, , drop = FALSE])

burn <- 10000
thin <- 10
lerca_short <- BurnThin(lerca, burn = burn, thin = thin)

cutoffs_keep <- lerca_short$cutoffs
coefs_keep <- lerca_short$coefs
variances_keep <- lerca_short$variances
alphas_keep <- lerca_short$alphas


# --------------- STEP 4. Diagnostics. ---------------- #

# Diagnostics for the cutoffs.
par(mfrow = c(2, 2), mar = rep(2, 4))

for (kk in 1 : K) {
  plot(1, type = 'n', ylim = range(cutoffs_keep[, , kk]),
       xlim = c(0, dim(cutoffs_keep)[2]))
  for (cc in 1 : chains) {
    lines(cutoffs_keep[cc, , kk], col = cc)
  }
}


# Diagnostics for alpha.
model <- 1
chain <- 2
round(apply(alphas_keep[model, chain, , , ], c(2, 3), mean), 3)
round(t(out_coef), 4)
round(t(XCcorr), 4)

par(mfrow = c(2, 2), mar = rep(2, 4))

jj <- 3
wh_model <- 1
for (kk in 1 : (K + 1)) {
  plot(1, xlim = c(0, dim(cutoffs_keep)[2]), ylim = c(0, 1), type = 'n')
  for (cc in 1 : chains) {
    lines(alphas_keep[wh_model, cc, , kk, jj], col = cc)
  }
}


# Diagnostics for the variances.
par(mfrow = c(2, 2), mar = rep(2, 4))

wh_model <- 1
for (ee in 1 : num_exper) {
  plot(1, xlim = c(0, dim(cutoffs_keep)[2]), type = 'n',
       ylim = range(variances_keep[wh_model, , , ee]))
  for (cc in 1 : chains) {
    lines(variances_keep[wh_model, cc, , ee], col = cc)
  }
}


# Diagnostics for the coefficients.
wh_coef <- 2
wh_model <- 2
for (ee in 1 : num_exper) {
  plot(1, xlim = c(0, dim(cutoffs_keep)[2]), type = 'n',
       ylim = range(coefs_keep[wh_model, , , ee, wh_coef], na.rm = TRUE))
  for (cc in 1 : chains) {
    lines(coefs_keep[wh_model, cc, , ee, wh_coef], col = cc)
  }
}



ER <- GetER(dta = dta, cutoffs = cutoffs_keep,
            coefs = coefs_keep[2, , , , ], mean_only = TRUE)



plot_ER <- apply(ER$y, c(1, 2), mean)
plot(1, type = 'n', ylim = range(plot_ER), xlim = range(ER$x))
for (cc in 1 : chains) {
  lines(ER$x, plot_ER[, cc], type = 'l', col = cc)
}


point_x <- c(10, 30, 50)
for (pp in 1 : length(point_x)) {
  plot(1, type = 'n', ylim = range(ER$y[point_x[pp], , ]), xlim = c(1, dim(ER$y)[3]))
  for (cc in 1 : chains) {
    lines(ER$y[point_x[pp], cc, ], col = cc)
  }
}

