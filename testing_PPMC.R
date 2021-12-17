# Testing Modified PPMC 

# The purpose of this script is to test the different functions developed to 
# implement the posterior predictive model checking methods available in the
# PPC.R file.

# Contents ----
# 0. Prepare Environment
# 1. Simulate data
# 2. Fit the AR-PCM model
# 3. Compute the PPMC methods

# 0. Prepare Environment ----
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)

source("R/IRT_models.R")
source("R/PPC.R")

# 1. Simulate data ----

# Simulate data based on an AR-PCM.
# Assumptions: Unidimensionality, local independence, stationarity, measurement
# invariance, AR process of order 1, and constant discrimination.
set.seed(123)

nT      <- 200 # Number of time points.  
lambda  <- 0.5 # Autoregressive effect.
I       <- 6   # Number of items.
K       <- 5   # Number of categories per items.
M       <- K - 1

theta <- rep(NA, nT)

theta[1] <- rnorm(1)

for (i in 2:nT) {
  theta[i] <- lambda * theta[i - 1] + rnorm(1, 0, 1) 
}
rm(i)

# Create item parameters
thresholds <- t(apply(matrix(runif(I * M, .3, 1), I), 1, cumsum))
if (M == 1) {thresholds <- t(thresholds)}
thresholds <- -(thresholds - rowMeans(thresholds))
thresholds <- thresholds + rnorm(I)
thresholds <- thresholds * -1

# Location
delta <- rowMeans(thresholds)

# Step parameters
taus <- thresholds - delta

# Generate responses
probs.array <- array(NA, dim = c(length(theta), I, K))

for (y in 0:M) {
  probs.array[, , y + 1] <- P.GPCM(y     = y, 
                                   alpha = rep(1, I), 
                                   delta = delta, 
                                   taus  = taus, 
                                   theta = theta, 
                                   M     = M)
}
responses   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
responses   <- responses + 1 # To fit the PCM in stan, items should be coded starting from 1. 
rm(probs.array, y)

# Generate missing values
na_prop <- 0 # proportion of missing values
na_ind  <- rbinom(dim(responses)[1], 1, na_prop)
responses[na_ind == 1, ] <- NA
rm(na_prop)

arpcm_data <- list(I  = I,
                   nT = nT,
                   K  = K,
                   M  = M,
                   na_ind    = na_ind,
                   responses = responses)
arpcm_par  <- list(theta      = theta,
                   taus       = taus,
                   delta      = delta,
                   thresholds = thresholds,
                   lambda     = lambda)
rm(theta, responses)

# Simulate data based on a bidimensional AR-PCM model

# Now, theta is a matrix with two factors that are weakly correlated.
theta <- matrix(NA, nrow = nT, ncol = 2)
mu    <- c(0, 0)
Sigma <- matrix(c(1, 0.3, 0.3, 1), 2)

theta[1, ] <- MASS::mvrnorm(1, mu = mu, Sigma = Sigma)

for (i in 2:nT) {
  theta[i, ] <- lambda * theta[i - 1] + MASS::mvrnorm(1, mu = mu, Sigma = Sigma) 
}
rm(i, mu, Sigma)

# Now, let's use the same item parameters to generate the data but using three
# items for one factor and three items for the other factor.

# Responses factor 1
probs.array <- array(NA, dim = c(dim(theta)[1], I / 2, K))

for (y in 0:M) {
  probs.array[, , y + 1] <- P.GPCM(y     = y, 
                                   alpha = rep(1, I / 2), 
                                   delta = delta[1:3], 
                                   taus  = taus[1:3, ], 
                                   theta = theta[, 1], 
                                   M     = M)
}
responses1   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
responses1   <- responses1 + 1 # To fit the PCM in stan, items should be coded starting from 1. 
rm(probs.array, y)

# Responses factor 2
probs.array <- array(NA, dim = c(dim(theta)[1], I / 2, K))

for (y in 0:M) {
  probs.array[, , y + 1] <- P.GPCM(y     = y, 
                                   alpha = rep(1, I / 2), 
                                   delta = delta[4:6], 
                                   taus  = taus[4:6, ], 
                                   theta = theta[, 2], 
                                   M     = M)
}
responses2   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
responses2   <- responses2 + 1 # To fit the PCM in stan, items should be coded starting from 1. 
rm(probs.array, y)

responses <- cbind(responses1, responses2)
rm(responses1, responses2)

# Generate missing values
responses[na_ind == 1, ] <- NA

bidim_data <- list(I  = I,
                   nT = nT,
                   K  = K,
                   M  = M,
                   na_ind    = na_ind,
                   responses = responses)
bidim_par  <- list(theta      = theta,
                   taus       = taus,
                   delta      = delta,
                   thresholds = thresholds,
                   lambda     = lambda)
rm(theta, responses, taus, thresholds, delta, lambda, na_ind)

#plot (apply(arpcm_data$responses[1:50, 1:3], 1, mean), type = "l")
#lines(apply(arpcm_data$responses[1:50, 4:6], 1, mean), col = 2)

# 2. Fit the AR-PCM model ----

# First we load the stan model that we are fitting.
# This is an old version of the ar-pcm model that does not account for non-linear
# trends and that assumes that variance of the innovations is fixed at 1.
model <- stan_model(file = "Stan/ar_irt_pcm_na.stan", verbose = TRUE)

# function to generate random initial values

arpcm_inits <- function() {
  list(lambda = runif(1, 0, 1),
       beta   = array(rnorm(24, 0, 3), dim = c(6, 4)),
       inno   = rnorm(200, 0, 3))
} 

arpcm_stan <- list(
  nT = arpcm_data$nT,
  I  = arpcm_data$I,
  K  = arpcm_data$K,
  N  = arpcm_data$nT * arpcm_data$I,
  N_obs = sum(!is.na(c(arpcm_data$responses))),
  tt = rep(1:arpcm_data$nT, arpcm_data$I),
  ii = rep(1:arpcm_data$I, each = arpcm_data$nT),
  tt_obs = rep(1:arpcm_data$nT, arpcm_data$I)[!is.na(c(arpcm_data$responses))],
  ii_obs = rep(1:arpcm_data$I, each = arpcm_data$nT)[!is.na(c(arpcm_data$responses))],
  y_obs  = c(arpcm_data$responses)[!is.na(c(arpcm_data$responses))]
  )

begin.time <- proc.time()
arpcm_fit <- sampling(
  model,                            # Stan model. 
  data = arpcm_stan,                  # Data.
  iter = 2000,                      # Number of iterations.
  chains  = 3,                      # Number of chains.
  warmup  = 1000,                   # Burn-in samples.
  init    = arpcm_inits,
  pars    = c("inno", "betasum", "log_lik", "imp_y"),
  include = FALSE,
  control = list(adapt_delta=0.95) # Other parameters to control sampling behavior.
  )
run.time <- proc.time() - begin.time
rm(begin.time)

bidim_stan <- list(
  nT = bidim_data$nT,
  I  = bidim_data$I,
  K  = bidim_data$K,
  N  = bidim_data$nT * bidim_data$I,
  N_obs = sum(!is.na(c(bidim_data$responses))),
  tt = rep(1:bidim_data$nT, bidim_data$I),
  ii = rep(1:bidim_data$I, each = bidim_data$nT),
  tt_obs = rep(1:bidim_data$nT, bidim_data$I)[!is.na(c(bidim_data$responses))],
  ii_obs = rep(1:bidim_data$I, each = bidim_data$nT)[!is.na(c(bidim_data$responses))],
  y_obs  = c(bidim_data$responses)[!is.na(c(bidim_data$responses))]
)

begin.time <- proc.time()
bidim_fit <- sampling(
  model,                            # Stan model. 
  data = bidim_stan,                  # Data.
  iter = 2000,                      # Number of iterations.
  chains  = 3,                      # Number of chains.
  warmup  = 1000,                   # Burn-in samples.
  init    = arpcm_inits,
  pars    = c("inno", "betasum", "log_lik", "imp_y"),
  include = FALSE,
  control = list(adapt_delta=0.95) # Other parameters to control sampling behavior.
)
run.time <- proc.time() - begin.time
rm(begin.time)

# Check convergence and parameter recovery.

# Replace here the model to check 
fit   <- bidim_fit
tpars <- bidim_par
tmp   <- list()

tmp$beta    <- summary(fit, pars = "beta")$summary
tmp$theta   <- summary(fit, pars = "theta")$summary
tmp$lambda  <- summary(fit, pars = "lambda")$summary

betapars  <- paste0("beta[", rep(c(1, ceiling(I / 2), I), each = 3), 
                    ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", ceiling(unname(quantile(1:nT, 
                                                      probs = seq(0, 1, length.out = 9)))), 
                    "]")

mcmc_rhat(rhat(fit))

traceplot(fit, pars = betapars, inc_warmup = FALSE)
traceplot(fit, pars = thetapars, inc_warmup = FALSE)
traceplot(fit, pars = "lambda", inc_warmup = FALSE)

plot(c(t(tpars$thresholds[, 1:M])), tmp$beta[, 1], pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(t(tpars$thresholds[, 1:M])), tmp$beta[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0  = c(t(tpars$thresholds)), 
         y0  = tmp$beta[, 4], 
         y1  = tmp$beta[, 8],
         col = rgb(0, 0, 0, 0.25))

plot(tpars$lambda, tmp$lambda[, 1], xlim = c(0, 1), ylim = c(0,1), pch = 16,
     ylab = "Estimated Autoregressive Effect", xlab = "True Autoregressive Effect")
abline(0, 1, col = 2, lwd = 2)
segments(x0  = tpars$lambda,
         y0  = tmp$lambda[, 4],
         y1  = tmp$lambda[, 8],
         col = rgb(0, 0, 0, 0.25))

# for bidim use tpars$theta[, 1], tpars$theta[, 2], or apply(tpars$theta, 1, mean)
plot(tpars$theta[, 2], type = "l", ylim = c(-3, 3), ylab = "Theta")
polygon(c(1:nT, rev(1:nT)),
        c(tmp$theta[, 4], rev(tmp$theta[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
lines(tmp$theta[, 1], col = "red", lwd = 1)

rm(tmp, tpars, fit, betapars, thetapars)

# 3. Compute the PPMC methods ----
par(mfrow = c(2, 3))
ppc.itcor(arpcm_fit, arpcm_stan, method = "polyserial", quiet =TRUE)
ppc.itcor(bidim_fit, bidim_stan, method = "polyserial", quiet =TRUE)

ppc.itcor2(arpcm_fit, arpcm_stan, method = "polyserial", quiet =TRUE)
ppc.itcor2(bidim_fit, bidim_stan, method = "polyserial", quiet =TRUE)

ppc.itcor3(arpcm_fit, arpcm_stan, quiet =TRUE)
ppc.itcor3(bidim_fit, bidim_stan, quiet =TRUE)

dev.off()
# END ----