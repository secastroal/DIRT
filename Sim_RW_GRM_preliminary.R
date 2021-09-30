# Simulation testing Dynamic IRT models
# Testing RW-GRM

# CONTENTS

# 0.0 Prepare Environment
# 1.0 Set up Conditions
# 2.0 Run simulation
# 3.0 Export output

# 0.0 Prepare Environment----
# Load required packages
library(doParallel)
library(foreach)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(bayesplot)

# Load required functions
source("R/IRT_models.R")

# Create folder to save output
if (!dir.exists(paste(getwd(), "Simulation", sep = "/"))) {
  dir.create(paste(getwd(), "Simulation", sep = "/"), recursive = TRUE)
}

# Create function needed to combine foreach output ----
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# 1.0 Set up Conditions----

# Fix conditions

K      <- 5   # Number of categories per item

# Manipulated conditions

N.timepoints <- c(60, 120, 200, 350, 500)
N.items      <- c(3, 6, 9, 12)

Cond        <- expand.grid(N.timepoints, N.items)
names(Cond) <- c("nT", "I")

rm(N.timepoints, N.items)

# Compile the model
model <- stan_model(file = "Stan/rw_irt_fixvar.stan", verbose = TRUE)

# 2.0 Run simulation ----

# Setup parallel backend to use parallel tasks
cl <- makeCluster(24)
registerDoParallel(cl, cores = 24)

# Get conditions and replications from the batch file
args <- commandArgs(trailingOnly = TRUE)

# Nested foreach loop (across conditions and across replications)

outcome.simulation <- foreach(cond = args[1]:args[2], .combine = 'list', .multicombine = TRUE) %:%
  foreach(r = args[3]:args[4], .combine = 'comb', .multicombine = TRUE, 
          .packages = c("rstan", "bayesplot")) %dopar% {
            
            # Define manipulated factors and seed
            nT <- Cond[cond, 1]
            I  <- Cond[cond, 2]
            seed <- 1000 * cond + r
            
            # Simulate data ----
            set.seed(seed)
            
            # True states
            theta    <- rep(NA, nT)
            theta[1] <- rnorm(1)
            
            for (i in 2:nT) {
              step <- rnorm(1, 0, 1)
              if (abs(theta[i - 1] + step) > 3) {
                theta[i] <- theta[i - 1]
              }else{
                theta[i] <- theta[i - 1] + step 
              } 
            }
            rm(i)
            
            # Create item parameters
            alpha      <- rlnorm(I, 0, 0.25)    # Discrimination parameters.
            delta      <- matrix(NA, I, K - 1)  # Matrix to store difficulty parameters.
            delta[, 1] <- rnorm(I, -2)
            for (k in 2:(K - 1)){
              delta[, k] <- delta[, k - 1] + runif(I, 0.4, 1.2)  #!#!# edited!
            }
            delta <- delta - rowMeans(delta)                     #!#!# edited!
            delta <- delta + rnorm(I)                            #!#!# edited!
            rm(k)
            
            IP <- cbind(delta, alpha) # Matrix of item parameters to input in P.GRM
            
            probs.array <- P.GRM(K - 1, IP, theta)
            probs.array[probs.array < 0] <- 0
            responses   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
            responses   <- responses + 1 # To fit the GRM in stan, items should be coded starting from 1. 
            rm(probs.array, delta, alpha)
            
            # Fit the RW-GRM model with stan ----
            standata <- list(nT = nT,
                             I  = I,
                             K  = K,
                             N  = nT * I,
                             tt = rep(1:nT, I),
                             ii = rep(1:I, each = nT),
                             y  = c(responses))
            begin.time <- proc.time()
            fit <- sampling(model,                            # Stan model. 
                            data = standata,                  # Data.
                            iter = 2000,                      # Number of iterations.
                            chains  = 3,                      # Number of chains.
                            warmup  = 1000,                   # Burn-in samples.
                            pars = c("alpha", "beta", "theta"),
                            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.
            run.time <- proc.time() - begin.time
            
            # Extract stan results ----
            sum.fit <- list()
            
            sum.fit$alpha    <- summary(fit, pars = "alpha")$summary
            sum.fit$beta     <- summary(fit, pars = "beta")$summary
            sum.fit$theta    <- summary(fit, pars = "theta")$summary
            
            #max Rhat
            maxRhat <- round(max(rhat(fit, pars = c("alpha", "beta", "theta"))), 4)
            nRhat   <- sum(rhat(fit, pars = c("alpha", "beta", "theta")) > 1.05)
            
            # correlations theta alpha and beta
            alpha.cor <- round(cor(IP[, K], sum.fit$alpha[, 1]), 4)
            beta.cor  <- round(cor(c(t(IP[, 1:(K - 1)])), sum.fit$beta[, 1]), 4)
            theta.cor <- round(cor(theta, sum.fit$theta[, 1]), 4)
            
            # bias theta alpha beta lambda
            alpha.bias <- round(sum(sum.fit$alpha[, 1] - IP[, K]) / I, 4)
            beta.bias  <- round(sum(sum.fit$beta[, 1] - c(t(IP[, 1:(K - 1)]))) / (I * (K - 1)), 4)
            theta.bias <- round(sum(sum.fit$theta[, 1] - theta) / nT, 4)
            
            # max CI theta alpha beta lambda
            alpha.CI  <-  round(max(sum.fit$alpha[, 8] - sum.fit$alpha[, 4]), 4)
            beta.CI   <-  round(max(sum.fit$beta[, 8] - sum.fit$beta[, 4]), 4)
            theta.CI  <-  round(max(sum.fit$theta[, 8] - sum.fit$theta[, 4]), 4)
            
            # Conditional to check whether the results are reliable (convergence and reasonable estimates)
            if (nRhat != 0   | alpha.cor < 0 | beta.cor < 0 | theta.cor < 0 | 
                alpha.CI > 2 | beta.CI > 4   | theta.CI > 4) {
              corrupt <- 1
            } else {
              corrupt <- 0
            }
            
            # Export output to outer foreach ----
            result <- list(c(cond,
                             r,
                             round(run.time[3]),
                             maxRhat,
                             nRhat,
                             alpha.cor,
                             alpha.bias,
                             alpha.CI,
                             beta.cor,
                             beta.bias,
                             beta.CI,
                             theta.cor,
                             theta.bias,
                             theta.CI,
                             corrupt))
            rm(standata, fit, sum.fit, IP, theta, responses,
               maxRhat, nRhat, alpha.cor, alpha.bias, alpha.CI,
               beta.cor, beta.bias, beta.CI, theta.cor, theta.bias,
               theta.CI, run.time, begin.time,
               corrupt)
            result
          }

# 3.0 Export output ----

for (i in args[1]:args[2]) {
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[1]]
  colnames(tmp) <- c("cond",
                     "r",
                     "run.time",
                     "maxRhat",
                     "nRhat",
                     "alpha.cor",
                     "alpha.bias",
                     "alpha.CI",
                     "beta.cor",
                     "beta.bias",
                     "beta.CI",
                     "theta.cor",
                     "theta.bias",
                     "theta.CI",
                     "corrupt")
  write.table(tmp, file = paste0(getwd(), "/Simulation/Sim_RW_GRM_cond_", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
}

# Stop parallel cluster
stopCluster(cl)
rm(cl)

# END