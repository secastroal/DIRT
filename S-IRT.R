# IRT-splines ----
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)
library("splines")

source("R/IRT_models.R")

## Optional ##
# If this script is executed through a batch file. One can specify
# additional arguments to change the number of time points, the number
# of knots, and the number of items.

#arg <- commandArgs(trailingOnly = TRUE)
#arg <- as.numeric(arg)

set.seed(123)

# In this section, I try to make an IRT-splines model. The idea is to assume that
# we have the time series of a test of a single individual. The relation between the 
# items' responses and the individual are modeled through the GRM. Moreover, the relation 
# between theta and time are modeled through a b-splines regression.

# First, let's assume the person was measured on 200 different occasions. Hence, our 
# dependent variable X is actually the time at which the person was measured.

time <- 1:200
# time <- 1:arg[1] # When batched.

# If the person was measured 10 times a day during 20 days, it makes sense to define 
# knots to separate measurements between days. This means 21 knots.
# We keep using cube splines.
num_knots     <- 21 # Define number of knots.
# num_knots     <- arg[2] # Define number of knots. When batched.
spline_degree <- 3  # Define degree of splines.
num_basis     <- num_knots + spline_degree - 1 # Compute number of basis.

num_data  <- length(time)
knots <- unname(quantile(time, probs = seq(0, 1, length.out = num_knots)))

# Compute b-splines within the time interval.
# Here, we generate b-splines of order fourth, which are cubic polynomial functions. These 
# splines are built based on num_knots, num_basis, and spline_degree.
B_true <- t(bs(time, df = num_basis, degree = spline_degree, intercept = TRUE)) # creating the B-splines

# Define intercept and coefficients to generate theta based on time and the b-splines.
a0 <- 0                      # intercept
a  <- rnorm(num_basis, 0, 1) # coefficients of B-splines
#!# How do we interpret this intercept? is it the trait theta?

# Compute theta based on time and b-splines. 
theta <- as.vector(a0 * time + a %*% B_true) # generating theta
#theta <- theta + rnorm(length(X),0,.2) 
#!# Should I add some noise on the relation between time and theta?

# Next, we generate data based on the GRM and the thetas we just created.

I <- 5 # Number of items.
# I <- arg[3] # Number of items. When batched.
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

# Fitting an IRT-splines model in stan. This first try intends to combine the splines model 
# in example 2 with the GRM model. 

standata <- list(num_data      = num_data, 
                 num_knots     = num_knots,
                 knots         = knots,
                 spline_degree = spline_degree,
                 I = I,
                 K = K,
                 time = time, 
                 Y = responses)

fit <- stan(file = "Stan/splines_irt.stan",   # Stan model. 
            data = standata,                  # Data.
            iter = 1000,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 500,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

sum.fit <- list()

sum.fit$alpha   <- summary(fit, pars = "alpha")$summary
sum.fit$beta    <- summary(fit, pars = "beta")$summary
sum.fit$theta   <- summary(fit, pars = "theta")$summary

betapars  <- paste0("beta[", rep(c(1, ceiling(I / 2), I), each = 3), 
                    ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", ceiling(unname(quantile(1:num_data, 
                                                      probs = seq(0, 1, length.out = 9)))), 
                    "]")

pdf(file = paste0("Splines_IRT_nT", num_data, 
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
     main = paste0("Discrimination; cor = ", round(cor(IP[, 5], sum.fit$alpha[, 1]), 3)))
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
plot(time, theta, pch = 20)
abline(v = knots, col = gray(0.75))
polygon(c(time, rev(time)),
        c(sum.fit$theta[, 4], rev(sum.fit$theta[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
lines(time, sum.fit$theta[, 1], col = "red", lwd = 2)

# Scatter plot of time and raw scores (average).
plot(time, apply(responses, 1, mean) - mean(apply(responses, 1, mean)), 
     type = "l", ylab = "Centered Raw Scores", lwd = 1.5)
lines(time, sum.fit$theta[, 1], col = "red", lwd = 2)

dev.off()
rm(list = ls())

standata <- list(num_data      = num_data, 
                 num_knots     = num_knots,
                 knots         = knots,
                 spline_degree = spline_degree,
                 I = I,
                 K = K,
                 N  = num_data * I,
                 tt = rep(1:num_data, I),
                 ii = rep(1:I, each = num_data),
                 y  = c(responses),
                 time = time)

fit.long <- stan(file = "Stan/splines_irt_long.stan",   # Stan model. 
                 data = standata,                  # Data.
                 iter = 1000,                      # Number of iterations.
                 chains  = 3,                      # Number of chains.
                 warmup  = 500,                    # Burn-in samples.
                 control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

