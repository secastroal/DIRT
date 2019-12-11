# In this document I try out some of the longitudinal IRT models.
# The next sections are included:
# 0. Prepare enviroment
# 1. Vandemeulebroecke et al (2017)
# 2. Ram et al (2005)


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
n  <- 100  # Number of subjects. 
p  <- 5    # Number of items. 
K  <- 5    # Number of categories per items.
nC <- 2    # Number of covariates.
nt <- 10   # Number of time points. #! In the paper not all persons were measured the same number of times.

# Now, we generate the thetas given the longitudinal model proposed by Vandemeulebroecke
# Priors are based on the supplementary Table 1 (Vandemeulebroecke et al., 2017).

gamma0    <- rnorm(n)                               # Random intercept, which varies over persons.
#J# Let's closely follow the priors (please check, Sebas!!):
#S# These are in fact the priors used in the original paper, which would result in generating gammastar
# as follows:
#gammastar_mu <- rnorm(1, 0, 2.5)
#gammastar_pr <- rgamma(1, shape=.2, rate=.2)
#gammastar_sd <- 1 / sqrt(gammastar_pr)
#gammastar    <- rnorm(n, gammastar_mu, gammastar_sd) 

# Given the stan code, gammastar should be generated as follows:
#gammastar <- gammastar_mu + gammastar_sd * rnorm(n)

# However, this way of generating gammastar can be problematic if gammastar_sd is too big. Therefore,
# here we generate gammastar from a normal standard distribution. Still, the priors of gammastar_mu and
# gammastar_sd are used in the stan model.
gammastar <- rnorm(n) # Random intercept of the random slope, which varies over persons.

#J# Probably it won't matter, but let's draw these from N(0,1), its prior:
beta   <- rnorm(nC) # Coefficients of the covariates.

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
fit.vande <- stan(file   = "Stan/vandemeulebroecke.stan", 
                  data   = standata,
                  iter   = 500, 
                  chains = 3, 
                  thin   = 1, 
                  cores  = 3, 
                  pars   = c("alpha", "kappa"))
time.vande <- proc.time() - t0
rm(t0)

traceplot(fit.vande, pars = "alpha", inc_warmup = TRUE)

summary(fit.vande, pars = "alpha")$summary

saveRDS(fit.vande, file = "Fits/fit.vande.rds")

detach("package:rstan")

# Fit vandemeleubroecke model in jags

library(rjags)

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

saveRDS(vande.fit.jags, file = "Fits/fit.vande.jags.rds")

loadRDS(file = "Fits/fit.vande.jags.rds")

# 2. Ram et al (2005) ----
# LIRT model proposed by Ram et al. (2005) is based on the rating scale model. 
# The longitudinal aspect is modeled through spectral analysis.

# Simulating data
# We will generate data given the spectral analysis model of Ram and the GRM.
# For this, we need: 
#   The parameters of the spectral model,
#   The parameters of the GRM,

set.seed(123)
# Let's define the conditions of the data
n  <- 200 # Number of subjects.
p  <- 8   # Number of items.
K  <- 5   # Number of categories per items.
nt <- 35  # Number of time points. #! not necessarily the same for all the persons

# Now, we generate the thetas given spectral analysis proposed by Ram
# We need the next random parameters:

mu    <- rnorm(n, -0.5, 1)  # Random intercept.
beta  <- rnorm(n, 0, .1)   # Random linear slope, which multiplies time_it. Detrend data.
R     <- rnorm(n, 0.2, 0.2) # Random amplitude of the cycle.
omega <- (2 * pi) / 7       # Period of oscillattion assumed to be 7 days.
phi   <- rnorm(n, 0.2, 0.4) # Random phase shift.

# Create matrix to store thetas over time
theta <- matrix(NA, n, nt)
for (i in 1:nt) {
  theta[, i] <- mu + R * (cos((omega * i) + phi)) + rnorm(1, 0, 0.5) #+ beta * (i - 1)
}
rm(i)

theta <- reshape(as.data.frame(theta), 
                 varying       = list(1:nt), 
                 new.row.names = 1:(n * nt),
                 direction     = "long") # Turn theta matrix into long format.
theta <- theta[order(theta[, 3]), c(3, 1, 2)]                                    # Order theta matrix by id.    

# Generate RSM responses given the generated thetas.

# Create item parameters
alpha <- rep(1, p)            # Discrimination parameters.
lambda <- sort(rnorm(p))
taus <- runif(K - 2, -1, 1) 
taus <- sort(c(taus, -sum(taus)))

probs.array       <- array(NA, dim = c(length(theta[, 3]), p, K))
for (y in 0:(K - 1)) 
{
  probs.array[, , y + 1] <- P.GPCM(y, alpha = alpha, delta = lambda, taus = taus, theta = theta[, 3], M = K - 1)
}
responses               <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
rm(probs.array)

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

nr <- nrow(simdatalong) # Number of observed responses.

# Run Ram's model

# There is no information on the paper about the conditions in which the analyses were run. They only
# mentioned that it took over 30 hours.


standata <- list(nr           = nr,
                 n            = n,
                 p            = p,
                 K            = K,
                 Y            = simdatalong[, "X"] + 1,
                 TP           = simdatalong[, "time"],
                 X1           = simdatalong[, "id"],
                 item         = simdatalong[, "item"]
                 )

t0 <- proc.time()
fit.ram <- stan(file = "Stan/ram.stan", 
                data = standata,
                pars = c("lambda", 
                         "delta",
                         "mu_gamma0",
                         "v_gamma0",
                         #"mu_gamma1",
                         #"v_gamma1",
                         "mu_gamma2",
                         "v_gamma2",
                         "mu_gamma3",
                         "v_gamma3",
                         "sigma_2"),
                iter   = 10000, 
                chains = 3, 
                thin   = 1, 
                cores  = 3)
time.ram <- proc.time() - t0
rm(t0)

traceplot(fit.ram, pars = "lambda", inc_warmup = TRUE)

summary(fit.ram, pars = "delta")$summary

saveRDS(fit.ram, file = "Fits/fit.ram.rds")

