# AR-IRT PCM ----
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)

source("R/IRT_models.R")

set.seed(123)

# This files simulates a simple AR-IRT model for a single individual.
# The Thetas over time are described by a lag 1 autoregressive model.
# The relation between the items and the person are modeled through a PCM model.

nT      <- 200 # Number of time points.  
lambda  <- 0.5 # Autoregressive effect.
#inn_var <- 0.7 # Variance of the innovation.


theta <- rep(NA, nT)

theta[1] <- rnorm(1)

for (i in 2:nT) {
  theta[i] <- lambda * theta[i - 1] + rnorm(1, 0, 1) 
}

# Next, we generate data based on the GRM and the thetas we just created.

I <- 10 # Number of items.
K <- 5  # Number of categories per items.
M <- K - 1

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

# Fitting an AR-IRT (PCM) model in stan. 

standata <- list(nT = nT,
                 I  = I,
                 K  = K,
                 N  = nT * I,
                 tt = rep(1:nT, I),
                 ii = rep(1:I, each = nT),
                 y  = c(responses))

fit <- stan(file = "Stan/ar_irt_pcm.stan",   # Stan model. 
            data = standata,                  # Data.
            iter = 2000,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 1000,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

sum.fit <- list()

sum.fit$beta    <- summary(fit, pars = "beta")$summary
sum.fit$theta   <- summary(fit, pars = "theta")$summary

betapars  <- paste0("beta[", rep(c(1, ceiling(I / 2), I), each = 3), 
                    ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", ceiling(unname(quantile(1:nT, 
                                                      probs = seq(0, 1, length.out = 9)))), 
                    "]")

pdf(file = paste0("AR_IRT_PCM_nT", nT, 
                  "_I", I, ".pdf"))
if (length(warnings()) != 0) {
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  mtext(names(warnings()), side = 3, line = (1:length(warnings())) * -5, adj = 0)
}

mcmc_rhat(rhat(fit))

traceplot(fit, pars = betapars, inc_warmup = FALSE)
traceplot(fit, pars = thetapars, inc_warmup = FALSE)

fit.array <- as.array(fit)

mcmc_acf(fit.array, 
         pars = betapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array, 
         pars = thetapars[c(1, 5, 9)], 
         lags = 20)

plot(c(t(thresholds[, 1:M])), sum.fit$beta[, 1], pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(t(thresholds[, 1:M])), sum.fit$beta[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = c(t(thresholds)), 
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
