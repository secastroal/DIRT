# In this document I try out some of the longitudinal IRT models.
# The next sections are included:
# 0. Prepare enviroment
# 1. Vandemeulebroecke et al (2017)


# 0. Prepare enviroment ----
rm(list = ls())
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)

#rfiles <- list.files("R/", full.names = TRUE)
#sapply(rfiles,source,.GlobalEnv)
#rm(rfiles)
source("R/IRT_models.R")

# 1. Vandemeulebroecke et al (2017) ----
# We will generate date based on the model proposed by Vandemeulebroecke et al (2017), 
# and fit the model to the data. 

# Simulating data
# For this, we need: 
#   The parameters of the model,
#   A dataset of covariates,
#   A dataset that includes responses of persons to items in long format.

set.seed(123)
# Let's define the conditions of the data
n  <- 100  # Number of subjects. 
p  <- 3    # Number of items. 
K  <- 5    # Number of categories per items.
nC <- 2    # Number of covariates.
nt <- 15   # Number of time points. #! In the paper not all persons were measured the same number of times.

# Now, we generate the thetas given the longitudinal model proposed by Vandemeulebroecke
# Priors are based on the supplementary Table 1 (Vandemeulebroecke et al., 2017).

gamma0    <- rnorm(n)                               # Random intercept, which varies over persons.
#J# Let's closely follow the priors (please check, Sebas!!):
#gammastar_mu <- rnorm(1, 0, 2.5)
#gammastar_pr <- rgamma(1, shape=.2, rate=.2)
#gammastar_sd <- 1 / sqrt(gammastar_pr)
gammastar    <- rnorm(n, 0, 0.01) # Random intercept of the random slope, which varies over persons.

#J# Probably it won't matter, but let's draw these from N(0,1), its prior:
#beta   <- rnorm(nC, 0, 0.005) # Coefficients of the covariates.
beta <- runif(nC, -0.05, 0.05)
#J# Z is not standardized. Better like this:
Z      <- replicate(nC, scale(rnorm(n))[1:n, ]) # Generate nC standardized covariates.
gamma1 <- gammastar + Z %*% beta                # Compute random slope.

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

#theta[, 3] <- rnorm(length(theta[, 3]), 0, 1)   

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

simdatalong <- reshape(simdata, 
                       idvar         = c("id", "time"),
                       varying       = list(3:(2 + p)),
                       new.row.names = 1:(n * nt * p),
                       direction     = "long") 

names(simdatalong) <- c("id", "item", "X")
time               <- rep(1:nt, n * p)
simdatalong        <- data.frame(simdatalong, time)
rm(time, simdata, responses)

simdatalong <- simdatalong[order(simdatalong[, 1], simdatalong[, 4]), c(1, 4, 2, 3)] # order data by id and time

nr <- nrow(simdatalong) # Number of observed responses. Might be fewer than p * nt * i given missing values.

obsxper <- tapply(simdatalong[,2], simdatalong[,1], function(x) length(unique(x)))
Tobs <- sum(obsxper)
indoxp <- rep(1:n, times = obsxper)


# Run Vandemeulebroecke's model

# In the paper, Vandemeulebroecke run 10000 iterations in stan (half burn-in) with 10 chains and no thining. 
# Nine out of the ten chains converge to the same results, one chain got stuck.The running time was 24 hours
# per chain. The reported effective size was 3700 on average using the function effectiveSize of the coda 
# package.

# Prepare data for stan

standata <- list(nr           = nr,                     # Number of rows.                
                 n            = n,                      # Number of persons.
                 p            = p,                      # Number of items.
                 K            = K,                      # Number of categories per item.
                 nC           = nC,                     # Number of standardized covariates.
                 Y            = simdatalong[, "X"] + 1, # Vector of responses.
                 TP           = simdatalong[, "time"],  # Vector of time indicators.
                 X1           = simdatalong[, "id"],    # Vector of ID.
                 item         = simdatalong[, "item"],  # Vector of item indicator.
                 Z            = Z,                      # Matrix of standardized covariates.
                 Tobs         = Tobs,                   # Total number ob observations.
                 indoxp       = indoxp,                 # Indicator of person per observation.
                 m_mu_gamma1  = 0,                      # Mean of the Hyperprior gamma 1.
                 sd_mu_gamma1 = 2.5,                    # SD of the Hyperprior gamma 1.
                 m_alpha      = 1,                      # Mean of the Hyperprior alpha.
                 sd_alpha     = 2.5,                    # SD of the Hyperprior alpha.
                 m_kappa      = 0,                      # Mean of the Hyperprior kappa.
                 sd_kappa     = 2.5,                    # SD of the Hyperprior kappa.
                 a_pr_gamma1  = 0.2,                    # shape of the Hyperprior gamma 1
                 b_pr_gamma1  = 0.2                     # rate of the Hyperprior gamma 1
)
t0 <- proc.time()
fit.vande <- stan(file   = "Stan/vandemeulebroecke_mlindex.stan", 
                  data   = standata,
                  pars   = c("alpha", "kappa", "theta"),
                  iter   = 5250,
                  warmup = 250, 
                  chains = 10, 
                  thin   = 5, 
                  cores  = 10)
time.vande <- proc.time() - t0
rm(t0)

sum.vande <- list()

sum.vande$alpha   <- summary(fit.vande, pars = "alpha")$summary
sum.vande$kappa   <- summary(fit.vande, pars = "kappa")$summary
#sum.vande$gamma0  <- summary(fit.vande, pars = "gamma0")$summary
#sum.vande$gamma1s <- summary(fit.vande, pars = "gamma1s")$summary
#sum.vande$gamma1  <- summary(fit.vande, pars = "gamma1")$summary
sum.vande$theta   <- summary(fit.vande, pars = "theta")$summary

kappapars <- paste0("kappa[", rep(c(1, ceiling(p / 2), p), each = 3), 
                    ",", rep(c(1, ceiling(K / 2), K - 1), times = 3), "]")
thetapars <- paste0("theta[", rep(c(0, ceiling(n / 2) - 1, n - 1), each = 3) * nt + 
                      rep(c(1, ceiling(nt / 2), nt), times = 3), "]")

pdf(file = paste0("Vandeplots_mlindex_N", n, "_nT", nt, "_I", p, ".pdf"))
if (length(warnings()) != 0) {
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  mtext(names(warnings()), side = 3, line = (1:length(warnings())) * -5, adj = 0)
}

mcmc_rhat(rhat(fit.vande))

traceplot(fit.vande, pars = "alpha", inc_warmup = TRUE)
traceplot(fit.vande, pars = kappapars, inc_warmup = TRUE)
traceplot(fit.vande, pars = thetapars, inc_warmup = TRUE)

fit.array <- as.array(fit.vande)

mcmc_acf(fit.array[,c(1, 5, 10),], 
         regex_pars = "alpha", 
         lags = 20)

mcmc_acf(fit.array[,c(1, 5, 10),], 
         pars = kappapars[c(1, 5, 9)], 
         lags = 20)

mcmc_acf(fit.array[,c(1, 5, 10),], 
         pars = thetapars[c(1, 5, 9)], 
         lags = 20)

plot(IP[,5], sum.vande$alpha[,1], pch = 20,
     xlab = "True alpha",
     ylab = "Estimated alpha",
     xlim = c(0.5, 2.5),
     ylim = c(0.5, 2.5),
     main = paste0("Discrimination; cor = ", round(cor(IP[,5], sum.vande$alpha[,1]), 3)))
abline(0, 1, col = 2, lwd = 2)

plot(c(t(IP[,1:4])), sum.vande$kappa[,1] / rep(sum.vande$alpha[,1], each = K - 1), pch = 20,
     xlab = "True kappa",
     ylab = "Estimated kappa",
     xlim = c(-3, 3),
     ylim = c(-3, 3),
     main = paste0("Thresholds; cor = ", 
                   round(cor(c(t(IP[,1:4])), sum.vande$kappa[,1] / rep(sum.vande$alpha[,1], each = K - 1)), 3)))
abline(0, 1, col = 2, lwd = 2)

#plot(gamma0, sum.vande$gamma0[,1], pch = 20,
 #    xlab = "True gamma 0",
  #   ylab = "Estimated gamma 0",
   #  main = paste0("Gamma 0; cor = ", 
    #               round(cor(gamma0, sum.vande$gamma0[,1]), 3)))
#abline(0, 1, col = 2, lwd = 2)

#plot(gammastar, sum.vande$gamma1s[,1], pch = 20,
 #    xlab = "True gamma star",
  #   ylab = "Estimated gamma star",
   #  main = paste0("Gamma star; cor = ", 
    #               round(cor(gammastar, sum.vande$gamma1s[,1]), 3)))
#abline(0, 1, col = 2, lwd = 2)

#plot(gamma1, sum.vande$gamma1[,1], pch = 20,
 #    xlab = "True gamma 1",
  #   ylab = "Estimated gamma 1",
   #  main = paste0("Gamma 1; cor = ", 
    #               round(cor(gamma1, sum.vande$gamma1[,1]), 3)))
#abline(0, 1, col = 2, lwd = 2)

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