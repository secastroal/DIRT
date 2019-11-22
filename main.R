# In this document I try out some of the longitudinal IRT models.
# The next sections are included:
# 0. Prepare enviroment
# 1. Vandemeulebroecke et al (2017)
# 2. Ram et al (2005)


# 0. Prepare enviroment ----
rm(list = ls())
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

rfiles <- list.files("R/", full.names = TRUE)
sapply(rfiles,source,.GlobalEnv)
rm(rfiles)

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
n  <- 100 # Number of subjects.
p  <- 5   # Number of items.
K  <- 5   # Number of categories per items.
nC <- 2   # Number of covariates.
nt <- 5  # Number of time points. #! not necessarily the same for all the persons

# Now, we generate the thetas given the longitudinal model proposed by Vandemeulebroecke

gamma0    <- rnorm(n)                             # Random intercept, which varies over persons.
gammastar <- rnorm(n, 0, 0.3)                     # Random intercept of the random slope, which varies over persons.
beta      <- sample((-3:3)/50, nC, replace = TRUE) # Coefficients of the covariates. 
Z         <- replicate(nC, rnorm(n))              # Generate nC standardized covariates. 

gamma1    <- gammastar + Z %*% beta               # Compute random slope.

# Create matrix to store thetas over time
theta <- matrix(NA, n, nt)
for (i in 1:nt) {
  theta[, i] <- gamma0 + gamma1 * i 
}
rm(i)

theta <- reshape(as.data.frame(theta), 
                 varying       = list(1:nt), 
                 new.row.names = 1:(n * nt),
                 direction     = "long") # Turn theta matrix into long format.
theta <- theta[order(theta[, 3]), c(3, 1, 2)]                                    # Order theta matrix by id.    

# Generate GRM responses given the generated thetas.

# Create item parameters
alpha      <- runif(p, 0.5 , 2)    # Discrimination parameters.
delta      <- matrix(NA, p, K - 1) # Matrix to store difficulty parameters.
delta[, 1] <- rnorm(p)
for (i in 2:(K - 1)){
  delta[, i] <- delta[, i - 1] + runif(1, 0.4, 0.9)
}
rm(i)

IP          <- cbind(delta, alpha) # Matrix of item parameters to input in P.GRM

probs.array <- P.GRM(K - 1, IP, theta[, 3])
responses   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
rm(probs.array, delta, alpha)

# Restructure data to add id, time, and item factors.
simdata <- data.frame(theta[, 1:2], responses)

simdatalong <- reshape(simdata, 
                       idvar         = c("id", "time"),
                       varying       = list(3:7),
                       new.row.names = 1:(n * nt * p),
                       direction     = "long") 

names(simdatalong) <- c("id", "item", "X")
time               <- rep(1:nt, n * p)
simdatalong        <- data.frame(simdatalong, time)
rm(time, simdata, responses)

simdatalong <- simdatalong[order(simdatalong[, 1], simdatalong[, 4]), c(1, 4, 2, 3)] # order data by id and time

nr <- nrow(simdatalong) # Number of observed responses. Might be fewer than p * nt * i given missing values.

# Run Vandemeulebroecke's model

# Prepare data for stan

standata <- list(nr           = nr,
                 n            = n,
                 p            = p,
                 K            = K,
                 nC           = nC,
                 Y            = simdatalong[, "X"] + 1,
                 TP           = simdatalong[, "time"],
                 X1           = simdatalong[, "id"],
                 item         = simdatalong[, "item"],
                 Z            = Z,
                 m_mu_gamma1  = 0,
                 sd_mu_gamma1 = 1,
                 m_alpha      = 0,
                 sd_alpha     = 1,
                 m_kappa      = 0,
                 sd_kappa     = 1,
                 a_pr_gamma1  = 1,
                 b_pr_gamma1  = 1
                 )

fit.vande <- stan(file = "Stan/vandemeulebroecke.stan", data = standata, iter = 2000, chains = 3, 
                  thin = 10, cores = 3)

traceplot(fit.vande, pars = "alpha[1]", inc_warmup = TRUE)

summary(fit.vande, pars = "alpha")$summary

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
n  <- 100 # Number of subjects.
p  <- 5   # Number of items.
K  <- 5   # Number of categories per items.
nt <- 10  # Number of time points. #! not necessarily the same for all the persons

# Now, we generate the thetas given spectral analysis proposed by Ram
# We need the next random parameters:

mu    <- rnorm(n, -0.5, 1)  # Random intercept.
beta  <- rnorm(n, 0, .01)   # Random linear slope, which multiplies time_it. Detrend data.
R     <- rnorm(n, 0.2, 0.2) # Random amplitude of the cycle.
omega <- (2 * pi) / 5       # Period of oscillattion assumed to be 7 days.
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
delta <- sort(rnorm(p))
taus <- runif(p - 2, -1, 1)
taus <- sort(c(taus, -sum(taus)))

probs.array       <- array(NA, dim = c(length(theta[, 3]), p, p))
for (y in 0:(p - 1)) 
{
  probs.array[, , y + 1] <- P.GPCM(y, alpha = alpha, delta = delta, taus = taus, theta = theta[, 3], M = p - 1)
}
responses               <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
rm(probs.array)

# Restructure data to add id, time, and item factors.
simdata <- data.frame(theta[, 1:2], responses)

simdatalong <- reshape(simdata, 
                       idvar         = c("id", "time"),
                       varying       = list(3:7),
                       new.row.names = 1:(n * nt * p),
                       direction     = "long") 

names(simdatalong) <- c("id", "item", "X")
time               <- rep(1:nt, n * p)
simdatalong        <- data.frame(simdatalong, time)
rm(time, simdata, responses)

simdatalong <- simdatalong[order(simdatalong[, 1], simdatalong[, 4]), c(1, 4, 2, 3)] # order data by id and time

nr <- nrow(simdatalong) # Number of observed responses.

# Run Ram's model

standata <- list(nr           = nr,
                 n            = n,
                 p            = p,
                 K            = K,
                 Y            = simdatalong[, "X"] + 1,
                 TP           = simdatalong[, "time"],
                 X1           = simdatalong[, "id"],
                 item         = simdatalong[, "item"]
                 )

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
                iter   = 2000, 
                chains = 3, 
                thin   = 10, 
                cores  = 3)

traceplot(fit.ram, pars = "lambda", inc_warmup = TRUE)

summary(fit.ram, pars = "delta")$summary


