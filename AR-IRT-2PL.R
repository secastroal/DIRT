# IRT-splines ----
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)

source("R/IRT_models.R")

set.seed(123)

# This files simulates a simple AR-IRT model for a single individual.
# The Thetas over time are described by a lag 1 autoregressive model.
# The relation between the items and the person are modeled through a 2PL model.

nT      <- 500 # Number of time points.  
lambda  <- 0.5 # Autoregressive effect.
#inn_var <- 0.7 # Variance of the innovation.


theta <- rep(NA, nT)

theta[1] <- rnorm(1)

for (i in 2:nT) {
  theta[i] <- lambda * theta[i - 1] + rnorm(1, 0, 1) 
}

# Next, we generate data based on the 2PL and the thetas we just created.

I <- 5 # Number of items.

# Create item parameters
alpha <- rlnorm(I, 0, 0.25)   # Discrimination parameters.
beta  <- sort(rnorm(I, 0, 1)) # Difficulty parameters. #They don´t need to be sorted.

responses <- gen.2PL(theta = theta, alpha = alpha, beta = beta)

# Fitting an AR-IRT model in stan. 

standata <- list(nT = nT,
                 I = I,
                 N = nT * I,
                 tt = rep(1:nT, I),
                 ii = rep(1:I, each = nT),
                 y = c(responses))

# Fit either the 2pl with ar_irt_2pl.stan or the 2 parameter
# normal ogive model with ar_irt_2pprobit.stan.
fit <- stan(file = "Stan/ar_irt_2pprobit.stan",    # Stan model. 
            data = standata,                  # Data.
            iter = 2000,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 1000,                   # Burn-in samples.
            control = list(adapt_delta=0.99)) # Other parameters to control sampling behavior.

sum.fit <- list()

sum.fit$alpha   <- summary(fit, pars = "alpha")$summary
sum.fit$beta    <- summary(fit, pars = "beta")$summary
sum.fit$theta   <- summary(fit, pars = "theta")$summary

thetapars <- paste0("theta[", ceiling(unname(quantile(1:nT, 
                                                      probs = seq(0, 1, length.out = 9)))), 
                    "]")

pdf(file = paste0("AR_IRT_2PL_nT", nT, 
                  "_I", I, ".pdf"))
if (length(warnings()) != 0) {
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  mtext(names(warnings()), side = 3, line = (1:length(warnings())) * -5, adj = 0)
}

mcmc_rhat(rhat(fit))

traceplot(fit, pars = "alpha", inc_warmup = TRUE)
traceplot(fit, pars = "beta", inc_warmup = TRUE)
traceplot(fit, pars = thetapars, inc_warmup = FALSE)

fit.array <- as.array(fit)

mcmc_acf(fit.array, 
         pars = paste0("alpha[",c(1, ceiling(I / 2), I), "]"), 
         lags = 20)

mcmc_acf(fit.array, 
         pars = paste0("beta[",c(1, ceiling(I / 2), I), "]"), 
         lags = 20)

mcmc_acf(fit.array, 
         pars = thetapars[c(1, 5, 9)], 
         lags = 20)

# If the probit model is used, the true item parameters have 
# to be divided by 1.7 to match the logit parameterization.
plot(alpha, sum.fit$alpha[, 1], pch = 20,
     xlab = "True alpha",
     ylab = "Estimated alpha",
     xlim = c(0, 1.2),
     ylim = c(0, 1.2),
     main = paste0("Discrimination; cor = ", round(cor(alpha, sum.fit$alpha[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = alpha, 
         y0 = sum.fit$alpha[, 4], 
         y1 = sum.fit$alpha[, 8],
         col = rgb(0, 0, 0, 0.25))

plot(beta, sum.fit$beta[, 1], pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(beta, sum.fit$beta[, 1]), 3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = beta, 
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
