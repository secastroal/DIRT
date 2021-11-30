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
rm(theta, responses, taus, thresholds, delta, I, K, lambda, M, na_ind, nT)

plot(apply(arpcm_data$responses[1:50,], 1, mean), type = "l")
lines(apply(bidim_data$responses[1:50,], 1, mean), col = 2)

# 2. Fit the AR-PCM model ----
# 3. Compute the PPMC methods ----

ppc.itcor(fit, standata, method = "polyserial", quiet =TRUE)
ppc.itcor2(fit, standata, method = "polyserial", quiet =TRUE)
ppc.itcor3(fit, standata, quiet =TRUE)

# END ----