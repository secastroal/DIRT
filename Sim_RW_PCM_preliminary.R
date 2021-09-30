# Simulation testing Dynamic IRT models
# Testing RW-PCM

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
M      <- K - 1

# Manipulated conditions

N.timepoints <- c(60, 120, 200, 350, 500)
N.items      <- c(3, 6, 9, 12)

Cond        <- expand.grid(N.timepoints, N.items)
names(Cond) <- c("nT", "I")

rm(N.timepoints, N.items)

# Compile the model
model <- stan_model(file = "Stan/rw_irt_pcm_fixvar.stan", verbose = TRUE)

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
            rm(i, step)
            
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
            
            # Fit the RW-PCM model with stan ----
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
                            pars = c("beta", "theta"),
                            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.
            run.time <- proc.time() - begin.time
            
            # Extract stan results ----
            sum.fit <- list()
            
            sum.fit$beta     <- summary(fit, pars = "beta")$summary
            sum.fit$theta    <- summary(fit, pars = "theta")$summary
            
            #max Rhat
            maxRhat <- round(max(rhat(fit, pars = c("beta", "theta"))), 4)
            nRhat   <- sum(rhat(fit, pars = c("beta", "theta")) > 1.05)
            
            # correlations theta alpha and beta
            beta.cor  <- round(cor(c(t(thresholds[, 1:M])), sum.fit$beta[, 1]), 4)
            theta.cor <- round(cor(theta, sum.fit$theta[, 1]), 4)
            
            # bias theta alpha beta lambda
            beta.bias  <- round(sum(sum.fit$beta[, 1] - c(t(thresholds[, 1:M]))) / (I * M), 4)
            theta.bias <- round(sum(sum.fit$theta[, 1] - theta) / nT, 4)
            
            # max CI theta alpha beta lambda
            beta.CI   <-  round(max(sum.fit$beta[, 8] - sum.fit$beta[, 4]), 4)
            theta.CI  <-  round(max(sum.fit$theta[, 8] - sum.fit$theta[, 4]), 4)
            
            # Conditional to check whether the results are reliable (convergence and reasonable estimates)
            if (nRhat != 0   | beta.cor < 0 | theta.cor < 0 | 
                beta.CI > 4  | theta.CI > 4) {
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
                             beta.cor,
                             beta.bias,
                             beta.CI,
                             theta.cor,
                             theta.bias,
                             theta.CI,
                             corrupt))
            rm(standata, fit, sum.fit, thresholds, delta, taus, theta, responses,
               maxRhat, nRhat, beta.cor, beta.bias, beta.CI, theta.cor, theta.bias,
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
                     "beta.cor",
                     "beta.bias",
                     "beta.CI",
                     "theta.cor",
                     "theta.bias",
                     "theta.CI",
                     "corrupt")
  write.table(tmp, file = paste0(getwd(), "/Simulation/Sim_RW_PCM_cond_", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
}

# Stop parallel cluster
stopCluster(cl)
rm(cl)

# END