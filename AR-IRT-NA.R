# AR-IRT GRM ----
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)

source("R/IRT_models.R")

set.seed(123)

# This files simulates a simple AR-IRT model for a single individual.
# The Thetas over time are described by a lag 1 autoregressive model.
# The relation between the items and the person are modeled through a GRM model.

nT      <- 200 # Number of time points.  
lambda  <- 0.5 # Autoregressive effect.
#inn_var <- 0.7 # Variance of the innovation.


theta <- rep(NA, nT)

theta[1] <- rnorm(1)

for (i in 2:nT) {
  theta[i] <- lambda * theta[i - 1] + rnorm(1, 0, 1) 
}

# Next, we generate data based on the GRM and the thetas we just created.

I <- 6  # Number of items.
K <- 5  # Number of categories per items.

# Create item parameters
alpha      <- rlnorm(I, 0, 0.25)    # Discrimination parameters.
delta      <- matrix(NA, I, K - 1)  # Matrix to store difficulty parameters.
delta[, 1] <- rnorm(I, -2)
for (i in 2:(K - 1)){
  delta[, i] <- delta[, i - 1] + runif(I, 0.4, 1.2)  #!#!# edited!
}
delta <- delta - rowMeans(delta)                     #!#!# edited!
delta <- delta + rnorm(I)                            #!#!# edited!
rm(i)

IP <- cbind(delta, alpha) # Matrix of item parameters to input in P.GRM

probs.array <- P.GRM(K - 1, IP, theta)
probs.array[probs.array < 0] <- 0
responses   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
responses   <- responses + 1 # To fit the GRM in stan, items should be coded starting from 1. 
rm(probs.array, delta, alpha)

# Introduce some missing values. Usually, in ESM data, the missing values happen
# for all observations at certain time points. For example, when participants
# miss a beep or because additional missing values are incorporated to account
# for unequal spaced observations. Hence, we make observations missing per
# row.

na_prop <- 0.1 # proportion of missing values
na_ind  <- rbinom(dim(responses)[1], 1, na_prop)
responses[na_ind == 1, ] <- NA

# Fitting an AR-IRT (GRM) model in stan when there are missing values.

# This model uses the available data to estimate all the parameters. 
# Even the thetas for which there were not observed data (NA rows).
# In the posterior simulations, the missing values can be simulated
# just as the observed values are simulated.

standata <- list(nT = nT,
                 I  = I,
                 K  = K,
                 N  = nT * I,
                 N_obs = sum(!is.na(c(responses))),
                 tt = rep(1:nT, I),
                 ii = rep(1:I, each = nT),
                 tt_obs = rep(1:nT, I)[!is.na(c(responses))],
                 ii_obs = rep(1:I, each = nT)[!is.na(c(responses))],
                 y_obs  = c(responses)[!is.na(c(responses))])

fit <- stan(file = "Stan/ar_irt_na.stan",   # Stan model. 
            data = standata,                  # Data.
            iter = 2000,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 1000,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

sum.fit <- list()

sum.fit$alpha   <- summary(fit, pars = "alpha")$summary
sum.fit$beta    <- summary(fit, pars = "beta")$summary
sum.fit$theta   <- summary(fit, pars = "theta")$summary

betapars  <- paste0("beta[", rep(c(1, ceiling(I / 2), I), each = 3), 
                    ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", ceiling(unname(quantile(1:nT, 
                                                      probs = seq(0, 1, length.out = 9)))), 
                    "]")

pdf(file = paste0("AR_IRT_nT", nT, 
                  "_I", I, ".pdf"))
if (length(warnings()) != 0) {
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  mtext(names(warnings()), side = 3, line = (1:length(warnings())) * -5, adj = 0)
}

mcmc_rhat(rhat(fit))

traceplot(fit, pars = "alpha", inc_warmup = TRUE)
traceplot(fit, pars = betapars, inc_warmup = FALSE)
traceplot(fit, pars = thetapars, inc_warmup = FALSE)

fit.array <- as.array(fit)

mcmc_acf(fit.array, 
         pars = paste0("alpha[",c(1, ceiling(I / 2), I), "]"), 
         lags = 20)

mcmc_acf(fit.array, 
         pars = betapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array, 
         pars = thetapars[c(1, 5, 9)], 
         lags = 20)

plot(IP[, K], sum.fit$alpha[, 1], pch = 20,
     xlab = "True alpha",
     ylab = "Estimated alpha",
     xlim = c(0, 2.5),
     ylim = c(0, 2.5),
     main = paste0("Discrimination; cor = ", round(cor(IP[, K], sum.fit$alpha[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = IP[, K], 
         y0 = sum.fit$alpha[, 4], 
         y1 = sum.fit$alpha[, 8],
         col = rgb(0, 0, 0, 0.25))

plot(c(t(IP[, 1:(K - 1)])), sum.fit$beta[, 1], pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(t(IP[, 1:(K - 1)])), sum.fit$beta[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = c(t(IP[, 1:(K - 1)])), 
         y0 = sum.fit$beta[, 4], 
         y1 = sum.fit$beta[, 8],
         col = rgb(0, 0, 0, 0.25))


plot(theta, sum.fit$theta[, 1], pch = 20,
     xlab = "True theta",
     ylab = "Estimated theta",
     xlim = c(-3, 3),
     ylim = c(-3, 3),
     main = paste0("Theta; cor = ", 
                   round(cor(theta, sum.fit$theta[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)

# Scatter plot of time and theta. Including the estimated theta in red (With 95% credibility interval).
plot(1:nT, theta, type = "l")
polygon(c(1:nT, rev(1:nT)),
        c(sum.fit$theta[, 4], rev(sum.fit$theta[, 8])),
        border = NA, 
        col = rgb(1, 0, 0, 0.25))
lines(1:nT, sum.fit$theta[, 1], col = "red", lwd = 1)

# Every 10th:
supp <- seq(0, nT, by = 10)[-1]
plot(supp, theta[supp], type = "l", main = "Every 10th time point")
polygon(c(supp, rev(supp)),
        c(sum.fit$theta[supp, 4], rev(sum.fit$theta[supp, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
lines(supp, sum.fit$theta[supp, 1], col = "red", lwd = 1)

dev.off()
rm(list = ls())
