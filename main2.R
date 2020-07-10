# Jorge: Trying out a new Stan model (vandemeulebroecke2.stan).
# 
# This document fit a variation of the longitudinal IRT model proposed 
# by Vandemeulebroecke(2017). The new stan model aims to optimize the 
# original one by not estimating the theta parameters repeatedly. For
# this, the model was rewritten such that it uses arrays. Data need 
# to be in array.

# 0. Prepare enviroment ----
rm(list = ls())
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# rfiles <- list.files("R/", full.names = TRUE)
# sapply(rfiles,source,.GlobalEnv)
# rm(rfiles)
source("R/IRT_models.R")

# Jorge:
# I tested the functions in IRT_models.R, just in case.
# It all seems ok.
# See R/testing_functions.R

# 1. Vandemeulebroecke et al (2017) ----
# We will generate date based on the model proposed by Vandemeulebroecke et al (2017), 
# and fit the model to the data. 

## Optional ##
# If this script is executed through a batch file. One can specify
# additional arguments to change the number of subjects, the number
# of time points, and the number of items.

#arg <- commandArgs(trailingOnly = TRUE)
#arg <- as.numeric(arg)

# Simulating data
# For this, we need: 
#   The parameters of the model,
#   A dataset of covariates,
#   A dataset that includes responses of persons to items in long format.

set.seed(123)
# Let's define the conditions of the data
n  <- 2000  # Number of subjects. 
p  <- 10    # Number of items. 
K  <- 5    # Number of categories per items.
nC <- 2    # Number of covariates.
nt <- 20   # Number of time points. #! In the paper not all persons were measured the same number of times.

## Optional ## # When batched.
# n  <- arg[1]
# nt <- arg[2]
# p  <- arg[3]

# Now, we generate the thetas given the longitudinal model proposed by Vandemeulebroecke

gamma0    <- rnorm(n)                # Random intercept, which varies over persons.
gammastar <- rnorm(n, 0, 0.03)       # Random intercept of the random slope, which varies over persons.
#J# Let's closely follow the priors (please check, Sebas!!):
# gammastar_mu <- rnorm(1, 0, 2.5)
# gammastar_pr <- rgamma(1, shape=.2, rate=.2)
# gammastar_sd <- 1 / sqrt(gammastar_pr)
# gammastar    <- gammastar_mu + gammastar_sd * rnorm(n)

beta      <- runif(nC, -0.05, 0.05)  # Coefficients of the covariates.
#J# Probably it won't matter, but let's draw these from N(0,1), its prior:
#beta <- rnorm(nC)

#J# Z is not standardized. Better like this:
Z      <- replicate(nC, scale(rnorm(n))[1:n, ])
gamma1 <- gammastar + Z %*% beta  # Compute random slope.

# Create matrix to store thetas over time
theta <- matrix(NA, n, nt)
for (i in 1:nt) {
  theta[, i] <- gamma0 + gamma1 * i 
}
rm(i)

theta <- reshape(as.data.frame(theta), 
                 varying       = list(1:nt), 
                 new.row.names = 1:(n * nt),
                 direction     = "long")      # Turn theta matrix into long format.
theta <- theta[order(theta[, 3]), c(3, 1, 2)] # Order theta matrix by id.    

# Generate GRM responses given the generated thetas.

# Create item parameters
alpha      <- rlnorm(p, 0, 0.25)    # Discrimination parameters.
delta      <- matrix(NA, p, K - 1)  # Matrix to store difficulty parameters.
delta[, 1] <- rnorm(p, -2)
for (i in 2:(K - 1)){
  delta[, i] <- delta[, i - 1] + runif(p, 0.4, 1.2)  #!#!# edited!
}
delta <- delta - rowMeans(delta)                     #!#!# edited!
delta <- delta + rnorm(p)                            #!#!# edited!
rm(i)

IP          <- cbind(delta, alpha) # Matrix of item parameters to input in P.GRM

probs.array <- P.GRM(K - 1, IP, theta[, 3])
probs.array[probs.array < 0] <- 0
responses   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
rm(probs.array, delta, alpha)

# Restructure data to add id, time, and item factors.
simdata <- data.frame(theta[, 1:2], responses)

#!!!#
simdata2 <- array(NA, c(n, nt, p))
for (i1 in 1:n) 
{
  for (i2 in 1:nt)
  {
    tmp <- which((simdata$id == i1) & (simdata$time == i2))
    simdata2[i1, i2, ] <- as.numeric(simdata[tmp, 3:(p+2)])
  }
}
#!!!#

# Run Vandemeulebroecke's model

# In the paper, Vandemeulebroecke run 10000 iterations in stan (half burn-in) with 10 chains and no thining. 
# Nine out of the ten chains converge to the same results, one chain got stuck.The running time was 24 hours
# per chain. The reported effective size was 3700 on average using the function effectiveSize of the coda 
# package.

# Prepare data for stan

standata <- list(n            = n,                      # Number of persons.
                 nt           = nt,                     # Number of time points
                 p            = p,                      # Number of items.
                 K            = K,                      # Number of categories per item.
                 nC           = nC,                     # Number of standardized covariates.
                 Y            = simdata2 + 1,           # Array of responses.
                 Z            = Z,                      # Matrix of standardized covariates.
                 m_mu_gamma1  = 0,                      # Mean of the prior of the mean of gamma 1.
                 sd_mu_gamma1 = 2.5,                    # SD of the prior of the mean gamma 1.
                 m_alpha      = 1,                      # Mean of the prior of alpha.
                 sd_alpha     = 2.5,                    # SD of the prior of alpha.
                 m_kappa      = 0,                      # Mean of the prior of kappa.
                 sd_kappa     = 2.5,                    # SD of the prior of kappa.
                 a_pr_gamma1  = 0.2,                    # shape of the prior of the sd of gamma 1
                 b_pr_gamma1  = 0.2                     # rate of the prior of the sd of gamma 1
                 )
t0 <- proc.time()

# There are two parameterizations of this stan model. The first one is 
# the classical IRT parameterization: alpha_j * (theta_i - delta_j). The
# second one is as follows: alpha_j * theta_i - kappa_j.
# Use the stan models "vandemeulebroecke2.stan" and 
# "vandemeulebroecke_oriarray.stan", respectively.

fit.vande <- stan(file   = "Stan/vandemeulebroecke_oriarray.stan", 
                  data   = standata,
                  iter   = 5250, 
                  warmup = 250, # it seems enough from previous runs
                  chains = 10, 
                  thin   = 5, 
                  cores  = 10, 
                  pars   = c("alpha", "kappa",  "beta_i", "theta"))
time.vande <- proc.time() - t0
rm(t0)

sum.vande <- list()

sum.vande$alpha    <- summary(fit.vande, pars = "alpha")$summary
sum.vande$kappa    <- summary(fit.vande, pars = "kappa")$summary
sum.vande$beta_i   <- summary(fit.vande, pars = "beta_i")$summary
sum.vande$theta    <- summary(fit.vande, pars = "theta")$summary

betapars <- paste0("beta_i[", rep(c(1, ceiling(p / 2), p), each = 3), 
                   ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", rep(c(1, ceiling(n / 2), n), each = 3), 
                    ",", rep(c(1, ceiling(nt / 2), nt), times = 3), "]")

pdf(file = paste0("Vandeplots_oriarray_N", n, "_nT", nt, "_I", p, ".pdf"))
if (length(warnings()) != 0) {
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  mtext(names(warnings()), side = 3, line = (1:length(warnings())) * -5, adj = 0)
}

mcmc_rhat(rhat(fit.vande))

traceplot(fit.vande, pars = "alpha", inc_warmup = TRUE)
traceplot(fit.vande, pars = betapars, inc_warmup = TRUE)
traceplot(fit.vande, pars = thetapars, inc_warmup = TRUE)

fit.array <- as.array(fit.vande)

mcmc_acf(fit.array[,c(1, 5, 10),], 
         regex_pars = "alpha", 
         lags = 20)

mcmc_acf(fit.array[,c(1, 5, 10),], 
         pars = betapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array[,c(1, 5, 10),], 
         pars = thetapars[c(1, 5, 9)], 
         lags = 20)

plot(IP[,5], sum.vande$alpha[,1], pch = 20,
     xlab = "True alpha",
     ylab = "Estimated alpha",
     xlim = c(0, 2),
     ylim = c(0, 2),
     main = paste0("Discrimination; cor = ", round(cor(IP[,5], sum.vande$alpha[,1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = IP[, 5], 
         y0 = sum.vande$alpha[, 4], 
         y1 = sum.vande$alpha[, 8],
         col = rgb(0, 0, 0, 0.25))

plot(c(t(IP[,1:4])), sum.vande$beta_i[,1], pch = 20,
     xlab = "True Locations",
     ylab = "Estimated Locations",
     xlim = c(-4, 4),
     ylim = c(-4, 4),
     main = paste0("Locations; cor = ", 
                   round(cor(c(t(IP[,1:4])), sum.vande$beta_i[,1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = c(t(IP[, 1:4])), 
         y0 = sum.vande$beta_i[, 4], 
         y1 = sum.vande$beta_i[, 8],
         col = rgb(0, 0, 0, 0.25))

plot(theta[, 3], sum.vande$theta[,1], pch = 20,
     xlab = "True theta",
     ylab = "Estimated theta",
     xlim = c(-10, 10),
     ylim = c(-10, 10),
     main = paste0("Theta; cor = ", 
                   round(cor(theta[,3], sum.vande$theta[,1]), 3)))
abline(0, 1, col = 2, lwd = 2)
dev.off()
rm(list = ls())
