# Dynamic Partial Credit Model with Time-varying Intercept ----

# In this file, we combine two models: The partial credit model and the time-varying
# autoregressive model (Bringmann et al., 2017).

# 0. Prepare environment ----
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)
library("splines")

source("R/IRT_models.R")
source("R/IRT_plots.R")

# 1. Data Simulation ----

# We assume the data represents the psychological time series of one individual.
# In this time series, the same items were repeatedly used to measure the subject.
# The measurement model between the observations and the latent construct is
# assumed to be a partial credit model. Moreover, at the latent level, we assume
# that the latent variable follows a time-varying autoregressive model with a 
# time-varying intercept. Here, the time-varying intercept is simulated and 
# modeled by means of a cubic splines model.
seed <- 123
set.seed(seed)

nT     <- 200   # Number of time points
I      <- 6     # Number of items
K      <- 5     # Number of categories per item
M      <- K - 1 # Number of thresholds per item
lambda <- 0.5   # Size of the autoregressive effect
in_var <- 1     # Variance of the innovations

# Simulate the time varying intercept with cubic splines
time     <- 1:nT
n_knots  <- 4  # Define number of knots.
s_degree <- 3  # Define degree of splines.
n_basis  <- n_knots + s_degree - 1 # Compute number of basis.

knots <- unname(quantile(time, probs = seq(0, 1, length.out = n_knots)))

# Create the b-splines within the time interval.
B_true <- t(bs(time, df = n_basis, degree = s_degree, intercept = TRUE))

# Define intercept and coefficients to generate the time-varying intercept based 
# on the b-splines.
a0 <- 0                    # intercept
a  <- rnorm(n_basis, 0, 1) # coefficients of B-splines

# Compute time-varying intercept based on time and b-splines. 
tv_int <- as.vector(a0 * time + a %*% B_true)

# Repeat lambda
lambda <- rep(lambda, nT)

# Generate the latent scores theta
theta <- rep(NA, nT)

# First theta based on an stationary marginal distribution (see Bringmann et al., 2017).
theta[1] <- rnorm(1, tv_int[1]/(1 - lambda[1]), sqrt(in_var/(1 - lambda[1] ^ 2)))

for (t in 2:nT) {
  theta[t] <- tv_int[t] + lambda[t] * theta[t - 1] + rnorm(1, 0, sqrt(in_var))
}
rm(t)

attractor <- tv_int / (1 - lambda)
p_var     <- in_var / (1 - lambda ^ 2)

# Next, we generate data based on the PCM and the thetas we just created.

# Create item parameters

# Thresholds
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

# 2. Fitting the true model in stan ----

# model <- stan_model(file = "Stan/ar_irt_pcm_na.stan", verbose = FALSE)
# model <- stan_model(file = "Stan/ar_irt_pcm_na_free.stan", verbose = FALSE)
# model <- stan_model(file = "Stan/splines_irt_pcm_long.stan", verbose = FALSE)
model <- stan_model(file = "Stan/tv_dpcm_int_v7.stan", verbose = FALSE)
# model <- stan_model(file = "Stan/tv_dpcm_slo.stan", verbose = FALSE)

# standata <- list(num_data      = nT, 
#                  num_knots     = n_knots,
#                  knots         = knots,
#                  spline_degree = s_degree,
#                  I = I,
#                  K = K,
#                  N  = nT * I,
#                  tt = rep(1:nT, I),
#                  ii = rep(1:I, each = nT),
#                  y  = c(responses),
#                  time = time)

standata <- list(nT = nT,
                 n_knots = n_knots,
                 knots = knots,
                 s_degree = s_degree,
                 I  = I,
                 K  = K,
                 N  = nT * I,
                 N_obs = sum(!is.na(c(responses))),
                 tt = rep(1:nT, I),
                 ii = rep(1:I, each = nT),
                 tt_obs = rep(1:nT, I)[!is.na(c(responses))],
                 ii_obs = rep(1:I, each = nT)[!is.na(c(responses))],
                 y_obs  = c(responses)[!is.na(c(responses))],
                 time = time)

arpcm_inits <- function() {
  list(lambda = runif(1, 0, 1),
       beta   = array(rnorm(I * (K - 1), 0, 3), dim = c(I, K - 1)),
       inno   = rnorm(nT, 0, 3))
}

# arpcm_inits <- function() {
#   list(a_raw  = rep(0.01, n_basis),
#        a0     = 0.01,
#        tau    = 0.01,
#        beta   = array(rnorm(I * (K - 1), 0, 3), dim = c(I, K - 1)),
#        inno   = rnorm(nT, 0, 3))
# }

begin.time <- proc.time()
fit <- sampling(model,                            # Stan model. 
                data = standata,                  # Data.
                iter = 2000,                      # Number of iterations.
                chains  = 3,                      # Number of chains.
                warmup  = 1000,                   # Burn-in samples.
                init    = arpcm_inits,
                seed    = 123,
                #pars = c("beta", "theta", "lambda"),
                control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.
run.time <- proc.time() - begin.time
rm(begin.time)

sum.fit <- list()

sum.fit$beta   <- summary(fit, pars = "beta")$summary
sum.fit$theta  <- summary(fit, pars = "theta")$summary
sum.fit$lambda <- summary(fit, pars = "lambda")$summary
sum.fit$sigma2 <- summary(fit, pars = "sigma2")$summary
sum.fit$p_var  <- summary(fit, pars = "p_var")$summary
sum.fit$attractor <- summary(fit, pars = "attractor")$summary

# 3. Checking convergence and parameter recovery ----

betapars  <- paste0("beta[", rep(c(1, ceiling(I / 2), I), each = 3), 
                    ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", trunc(seq(1, nT, length.out = 9)), "]")
attrapars <- paste0("attractor[", trunc(seq(1, nT, length.out = 9)), "]")

pdf(file = paste0("Figures/TV_DPCM_nT", nT, 
                  "_I", I, "_seed", seed, ".pdf"))
if (length(warnings()) != 0) {
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  mtext(names(warnings()), side = 3, line = (1:length(warnings())) * -5, adj = 0)
}

# Plot Rhats
mcmc_rhat(rhat(fit))

# Traceplots of selected parameters
traceplot(fit, pars = betapars, inc_warmup = FALSE)
traceplot(fit, pars = thetapars, inc_warmup = FALSE)
traceplot(fit, pars = attrapars, inc_warmup = FALSE)
traceplot(fit, pars = c("lambda", "sigma2", "p_var"), inc_warmup = FALSE)

# Autocorrelation plots selected parameters
fit.array <- as.array(fit)

mcmc_acf(fit.array, 
         pars = betapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array, 
         pars = thetapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array, 
         pars = attrapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array, 
         pars = c("lambda", "sigma2", "p_var"), 
         lags = 20)

# Recovery of parameters: True vs. estimated.
plot(c(t(thresholds)), sum.fit$beta[, 1], pch = 20,
     xlab = "True beta",
     ylab = "Estimated beta",
     xlim = c(-3.5, 3.5),
     ylim = c(-5, 5),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(t(thresholds)), sum.fit$beta[, 1]), 3)))
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

hist(fit.array[, , "lambda"], freq = FALSE)
abline(v = lambda[1], lwd = 5, col = 2)
abline(v = sum.fit$lambda[, 1], lwd = 5, col = 1)
lines(sum.fit$lambda[, c(4, 8)], c(0, 0), lwd = 5, col = 3)

hist(fit.array[, , "sigma2"], freq = FALSE)
abline(v = in_var, lwd = 5, col = 2)
abline(v = sum.fit$sigma2[, 1], lwd = 5, col = 1)
lines(sum.fit$sigma2[, c(4, 8)], c(0, 0), lwd = 5, col = 3)

hist(fit.array[, , "p_var"], freq = FALSE)
abline(v = p_var[1], lwd = 5, col = 2)
abline(v = sum.fit$p_var[, 1], lwd = 5, col = 1)
lines(sum.fit$p_var[, c(4, 8)], c(0, 0), lwd = 5, col = 3)

# Scatter plot of time and theta. 
# Including the estimated theta in red (With 95% credibility interval) and the 
# estimated attractor in blue.
plot(theta, pch = 20)
abline(v = knots, col = gray(0.75))
polygon(c(time, rev(time)),
        c(sum.fit$theta[, 4], rev(sum.fit$theta[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
polygon(c(time, rev(time)),
        c(sum.fit$attractor[, 4], rev(sum.fit$attractor[, 8])),
        border = NA,
        col = rgb(0, 0, 1, 0.15))
lines(time, sum.fit$theta[, 1], col = "red", lwd = 2)
lines(time, attractor, col = gray(0.5), lwd = 2)
lines(time, sum.fit$attractor[, 1], col = "blue", lwd = 2)

# Same plot but only 25 points
supp <- trunc(seq(1, nT, length.out = 25))
plot(supp, theta[supp], pch = 20, 
     main = "25 time points", ylim = c(min(theta) - 0.5, max(theta) + 0.5))
abline(v = knots, col = gray(0.75))
polygon(c(supp, rev(supp)),
        c(sum.fit$theta[supp, 4], rev(sum.fit$theta[supp, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
polygon(c(supp, rev(supp)),
        c(sum.fit$attractor[supp, 4], rev(sum.fit$attractor[supp, 8])),
        border = NA,
        col = rgb(0, 0, 1, 0.15))
lines(supp, sum.fit$theta[supp, 1], col = "red", lwd = 2)
lines(supp, attractor[supp], col = gray(0.5), lwd = 2)
lines(supp, sum.fit$attractor[supp, 1], col = "blue", lwd = 2)

# ICCs and IIFs
plot.ICC(fit, standata, quiet = TRUE)

plot.IIF(fit, standata, type = "IIF")

plot.IIF(fit, standata, type = "TIF")

dev.off()

assign("last.warning", NULL, envir = baseenv())
rm(list = setdiff(ls(), lsf.str()))

# 4. End ----

