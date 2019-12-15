# Jorge: Trying out a new Stan model (vandemeulebroecke2.stan).
# 
# In this document I try out some of the longitudinal IRT models.
# The next sections are included:
# 0. Prepare enviroment
# 1. Vandemeulebroecke et al (2017)



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
beta <- rnorm(nC)

Z         <- replicate(nC, rnorm(n)) # Generate nC standardized covariates. 
#J# Z is not standardized. Better like this:
Z <- replicate(nC, scale(rnorm(n))[1:n, ])

gamma1    <- gammastar + Z %*% beta  # Compute random slope.

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
fit.vande <- stan(file   = "Stan/vandemeulebroecke2.stan", 
                  data   = standata,
                  iter   = 2000, 
                  warmup = 200, # it seems enough from previous runs
                  chains = 3, 
                  thin   = 1, 
                  cores  = 3, 
                  pars   = c("alpha", "kappa",  "beta", "theta"))
time.vande <- proc.time() - t0
rm(t0)

# Warning message:
#   Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 

print(object.size(fit.vande), units="Mb") # 470Mb (n=2000, p=10, nt=10, nC=2, K=5, 
                                          #        iter=1000, chains=3, thin=1, cores=3, 
                                          #        pars = c("alpha", "kappa", "theta"))
                                          # 1855Mb (n=2000, p=10, nt=20, nC=2, K=5, 
                                          #        iter=2000, chains=3, thin=1, cores=3, 
                                          #        pars = c("alpha", "kappa", "beta", "theta")), 18.5hr

saveRDS(fit.vande, "fit_vande.rds")

traceplot(fit.vande, pars = "alpha", inc_warmup = TRUE)
traceplot(fit.vande, pars = "kappa", inc_warmup = TRUE)
traceplot(fit.vande, pars = "beta",  inc_warmup = TRUE)
traceplot(fit.vande, pars = "theta", inc_warmup = TRUE)

summary(fit.vande, pars = "alpha")$summary
plot(IP[, 5], summary(fit.vande, pars = "alpha")$summary[, "mean"])
abline(0, 1, col = 2, lwd = 2)

summary(fit.vande, pars = "kappa")$summary
plot(c(t(IP[, 1:4])), summary(fit.vande, pars = "kappa")$summary[, "mean"])
abline(0, 1, col = 2, lwd = 2)

beta
summary(fit.vande, pars = "beta")$summary

summary(fit.vande, pars = "theta")$summary
plot(theta[, 3], summary(fit.vande, pars = "theta")$summary[, "mean"])
abline(0, 1, col = 2, lwd = 2)

save(fit.vande, file = "Fits/fit.vande.RData")




#J#: NOT edited below.
# Fit vandemeleubroecke model in jags

jagsdata <- list(nr           = nr,                     # Number of rows.                
                 n            = n,                      # Number of persons.
                 p            = p,                      # Number of items.
                 K            = K,                      # Number of categories per item.
                 nC           = nC,                     # Number of standardized covariates.
                 Y            = simdatalong[, "X"] + 1, # Vector of responses.
                 TP           = simdatalong[, "time"],  # Vector of time indicators.
                 X1           = simdatalong[, "id"],    # Vector of ID.
                 item         = simdatalong[, "item"],  # Vector of item indicator.
                 Z            = Z,                      # Matrix of standardized covariates.
                 m.mu.gamma1  = 0,                      # Mean of the Hyperprior gamma 1.
                 s.mu.gamma1  = 2.5,                    # SD of the Hyperprior gamma 1.
                 m.alpha      = 1,                      # Mean of the Hyperprior alpha.
                 s.alpha      = 2.5,                    # SD of the Hyperprior alpha.
                 m.kappa      = 0,                      # Mean of the Hyperprior kappa.
                 s.kappa      = 2.5,                    # SD of the Hyperprior kappa.
                 a.pr.gamma1  = 0.2,                    # shape of the Hyperprior gamma 1
                 b.pr.gamma1  = 0.2                     # rate of the Hyperprior gamma 1
)

# Compile the model:
jagscompiled <- jags.model("jags/vandemeulebroecke.txt", 
                           data     = jagsdata, 
                           n.chains = 3, 
                           n.adapt  = 1000)

# Warm up:
update(jagscompiled, 1000)

 # Draw samples:
vande.fit.jags <- coda.samples(jagscompiled, 
                               data           = jagsdata, 
                               variable.names = c("alpha", "kappa"), 
                               n.iter         = 1000, 
                               thin           = 1)

save(vande.fit.jags, file = "Fits/fit.vande.jags.RData")

load(file = "Fits/fit.vande.jags.RData")

