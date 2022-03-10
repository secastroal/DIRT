# Fitting the tv-dpcm in jags

library(splines)
library(rjags)
library(runjags)

# Load required functions
source("R/IRT_models.R")
source("R/DBDA2E-utilities.R")

# Generate data ----
nT       <- 200  # Number of time points
I        <- 6    # Number of items
K        <- 5    # Number of categories per item
M        <- K - 1
in_var   <- 1    # Variance of the innovations
lambda   <- 0.25 # Autoregressive effect

seed   <- 1000

set.seed(seed)

# True states
time <- 1:nT

# Create time-varying intercept based on cosine. 
tv_int <- cos(3 * pi * time/(nT))

# Repeat lambda
tv_lambda <- rep(lambda, nT)

# Generate the latent scores theta
theta <- rep(NA, nT)

# First theta based on an stationary marginal distribution 
# (see Bringmann et al., 2017).
theta[1] <- rnorm(1, tv_int[1]/(1 - tv_lambda[1]), sqrt(in_var/(1 - tv_lambda[1] ^ 2)))

for (t in 2:nT) {
  theta[t] <- tv_int[t] + tv_lambda[t] * theta[t - 1] + rnorm(1, 0, sqrt(in_var))
}
rm(t)

attractor <- tv_int / (1 - tv_lambda)
p_var     <- in_var / (1 - lambda ^ 2)

# Create item parameters
thresholds <- t(apply(matrix(runif(I * M, .3, 1), I), 1, cumsum))
if (M == 1) {thresholds <- t(thresholds)}
thresholds <- thresholds - rowMeans(thresholds)
thresholds <- thresholds + rnorm(I)

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


n_knots  <- 8  # Define number of knots.
s_degree <- 3  # Define degree of splines.
n_basis  <- n_knots + s_degree - 1 # Compute number of basis.

jags_data <- tvdpcm2jags_data2(responses,
                               I = I,
                               K = K,
                               time     = time, 
                               n_basis  = n_basis,
                               s_degree = s_degree)

j.inits <- function() {
  list(lambda    = runif(1, -1, 1),
       alpha_raw = rnorm(n_basis, 0, 3),
       #beta      = cbind(NA, array(rnorm(I * (K - 1), 0, 3), dim = c(I, K - 1))),
       pr_alpha  = rgamma(1, .1, .1),
       pr_inno   = rgamma(1, .1, .1))
}

# run model in runjags-rjags

t0 <- proc.time()
fit <- run.jags(model = "jags/tv_dpcm.jags",
                monitor = c("beta", "theta", "attractor",
                            "lambda", "sigma2", "p_var"),
                data = jags_data,
                method = "rjags",
                #inits = j.inits,
                n.chains = 3,
                adapt = 1000,
                burnin = 1000,
                sample = 5000,
                thin = 1,
                summarise = FALSE,
                plots = FALSE)
t1 <- proc.time() - t0

fit <- as.mcmc.list(fit)

fit_sum <- summary(fit)

# Diagnostics
acor  <- autocorr.diag(fit[, dimnames(fit[[1]])[[2]][-(1:I)]], lags = 1)
rhats <- gelman.diag(fit[, dimnames(fit[[1]])[[2]][-(1:I)]], 
                     multivariate = FALSE, 
                     autoburnin = FALSE)
summary(as.vector(acor))
summary(rhats$psrf[, 1])

diagMCMC(fit, parName = "beta[1,2]")
diagMCMC(fit, parName = "beta[3,4]")
diagMCMC(fit, parName = "beta[6,5]")
diagMCMC(fit, parName = "theta[1]")
diagMCMC(fit, parName = "theta[100]")
diagMCMC(fit, parName = "theta[200]")
diagMCMC(fit, parName = "attractor[2]")
diagMCMC(fit, parName = "attractor[99]")
diagMCMC(fit, parName = "attractor[199]")
diagMCMC(fit, parName = "lambda")
diagMCMC(fit, parName = "sigma2")
diagMCMC(fit, parName = "p_var")

plotPost(fit[, "lambda"], cenTend = "median", xlab = "lambda",
         compVal = lambda)
plotPost(fit[, "sigma2"], cenTend = "median", xlab = "sigma2",
         compVal = in_var)
plotPost(fit[, "p_var"], cenTend = "median", xlab = "p_var",
         compVal = p_var)

plot(c(thresholds), 
     fit_sum$statistics[grep("beta", row.names(fit_sum$statistics)), 1][-(1:6)], 
     pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(thresholds), 
                             fit_sum$statistics[grep("beta", row.names(fit_sum$statistics)), 1][-(1:6)]), 
                         3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = c(thresholds), 
         y0 = fit_sum$quantiles[grep("beta", row.names(fit_sum$statistics)), 1][-(1:I)], 
         y1 = fit_sum$quantiles[grep("beta", row.names(fit_sum$statistics)), 5][-(1:I)],
         col = rgb(0, 0, 0, 0.25))

plot(theta, pch = 20, ylim =c(-4, 4))
polygon(c(time, rev(time)),
        c(fit_sum$quantiles[grep("theta", row.names(fit_sum$statistics)), 1], 
          rev(fit_sum$quantiles[grep("theta", row.names(fit_sum$statistics)), 5])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
polygon(c(time, rev(time)),
        c(fit_sum$quantiles[grep("attractor", row.names(fit_sum$statistics)), 1], 
          rev(fit_sum$quantiles[grep("attractor", row.names(fit_sum$statistics)), 5])),
        border = NA,
        col = rgb(0, 0, 1, 0.15))
lines(time, fit_sum$statistics[grep("theta", row.names(fit_sum$statistics)), 1], 
      col = "red", lwd = 2)
lines(time, attractor, col = gray(0.5), lwd = 2)
lines(time, fit_sum$statistics[grep("attractor", row.names(fit_sum$statistics)), 1], 
      col = "blue", lwd = 2)

# plot(fit[, "lambda"])
# traceplot(fit[, "lambda"])
# cumuplot(fit[, "lambda"])
# autocorr.plot(fit[, "lambda"], lag.max = 30)
# densplot(fit[, "lambda"])
# gelman.plot(fit[, "lambda"])

# Try the tv-dpcm with jagam ----

library(mgcv)
source("R/tvdpcm2jags.R")

jags_data <- tvdpcm2jags_data(resp = responses,
                              I    = I,
                              K    = K,
                              time = time,
                              n_basis = 30,
                              bs   = "tp")

t0 <- proc.time()
fit_jagam <- run.jags(model = "jags/tv_dpcm_jagam.jags",
                 monitor = c("beta", "theta", "attractor",
                             "lambda", "sigma2", "p_var"),
                 data = jags_data,
                 method = "rjags",
                 #inits = j.inits,
                 n.chains = 3,
                 adapt = 1000,
                 burnin = 1000,
                 sample = 5000,
                 thin = 1,
                 summarise = FALSE,
                 plots = FALSE)
t1 <- proc.time() - t0

fit <- as.mcmc.list(fit_jagam)

fit_sum <- summary(fit)

# Diagnostics
acor  <- autocorr.diag(fit[, dimnames(fit[[1]])[[2]][-(1:I)]], lags = 1)
rhats <- gelman.diag(fit[, dimnames(fit[[1]])[[2]][-(1:I)]], 
                     multivariate = FALSE, 
                     autoburnin = FALSE)
summary(as.vector(acor))
summary(rhats$psrf[, 1])

diagMCMC(fit, parName = "beta[1,2]")
diagMCMC(fit, parName = "beta[3,4]")
diagMCMC(fit, parName = "beta[6,5]")
diagMCMC(fit, parName = "theta[1]")
diagMCMC(fit, parName = "theta[100]")
diagMCMC(fit, parName = "theta[200]")
diagMCMC(fit, parName = "attractor[2]")
diagMCMC(fit, parName = "attractor[99]")
diagMCMC(fit, parName = "attractor[199]")
diagMCMC(fit, parName = "lambda")
diagMCMC(fit, parName = "sigma2")
diagMCMC(fit, parName = "p_var")

plotPost(fit[, "lambda"], cenTend = "median", xlab = "lambda",
         compVal = lambda)
plotPost(fit[, "sigma2"], cenTend = "median", xlab = "sigma2",
         compVal = in_var)
plotPost(fit[, "p_var"], cenTend = "median", xlab = "p_var",
         compVal = p_var)

plot(c(thresholds), 
     fit_sum$statistics[grep("beta", row.names(fit_sum$statistics)), 1][-(1:6)], 
     pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(thresholds), 
                             fit_sum$statistics[grep("beta", row.names(fit_sum$statistics)), 1][-(1:6)]), 
                         3)))
abline(0, 1, col = 2, lwd = 2)
segments(x0 = c(thresholds), 
         y0 = fit_sum$quantiles[grep("beta", row.names(fit_sum$statistics)), 1][-(1:I)], 
         y1 = fit_sum$quantiles[grep("beta", row.names(fit_sum$statistics)), 5][-(1:I)],
         col = rgb(0, 0, 0, 0.25))

plot(theta, pch = 20, ylim =c(-4, 4))
polygon(c(time, rev(time)),
        c(fit_sum$quantiles[grep("theta", row.names(fit_sum$statistics)), 1], 
          rev(fit_sum$quantiles[grep("theta", row.names(fit_sum$statistics)), 5])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
polygon(c(time, rev(time)),
        c(fit_sum$quantiles[grep("attractor", row.names(fit_sum$statistics)), 1], 
          rev(fit_sum$quantiles[grep("attractor", row.names(fit_sum$statistics)), 5])),
        border = NA,
        col = rgb(0, 0, 1, 0.15))
lines(time, fit_sum$statistics[grep("theta", row.names(fit_sum$statistics)), 1], 
      col = "red", lwd = 2)
lines(time, attractor, col = gray(0.5), lwd = 2)
lines(time, fit_sum$statistics[grep("attractor", row.names(fit_sum$statistics)), 1], 
      col = "blue", lwd = 2)

