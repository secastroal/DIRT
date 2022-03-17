# Simulation testing Dynamic IRT models
# Testing AR-PCM

# CONTENTS

# 0.0 Prepare Environment
# 1.0 Set up Conditions
# 2.0 Run simulation
# 3.0 Export output

# 0.0 Prepare Environment----
# Load required packages
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
suppressPackageStartupMessages(library(bayesplot))

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

N.timepoints <- c(100, 200, 350, 500) # Number of timepoints
N.items      <- c(3, 6)               # Number of items 
S.lambda     <- c(0, 0.25, 0.5)       # Size of the autoregressive effect
M.prop       <- c(0, 0.15, 0.25)      # Proportion of missing values

Cond        <- expand.grid(N.timepoints, N.items, S.lambda, M.prop)
names(Cond) <- c("nT", "I", "lambda", "NAprop")

rm(N.timepoints, N.items, S.lambda, M.prop)

# Compile the model
# To run the simulation we commented out the generated quantities block to 
# reduce memory usage.
model <- stan_model(file = "Stan/ar_irt_pcm_na.stan", verbose = TRUE)

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
            nT     <- Cond[cond, 1]
            I      <- Cond[cond, 2]
            lambda <- Cond[cond, 3]
            naprop <- Cond[cond, 4]
            seed   <- 1000 * cond + r
            
            # Simulate data ----
            set.seed(seed)
            
            # True states
            theta    <- rep(NA, nT)
            theta[1] <- rnorm(1)
            
            for (i in 2:nT) {
              theta[i] <- lambda * theta[i - 1] + rnorm(1, 0, 1) 
            }
            rm(i)
            
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
            
            # Introduce missing values per "beep" (random rows with complete NAs)
            
            if (naprop != 0) {
              naind <- rbinom(dim(responses)[1], 1, naprop)
              responses[naind == 1, ] <- NA
              rm(naind)
            }
                        
            # Fit the AR-PCM model with stan ----
            standata <- list(nT = nT,
                             I  = I,
                             K  = K,
                             N  = nT * I,
                             N_obs = sum(!is.na(c(responses))),
                             tt = rep(1:nT, I),
                             ii = rep(1:I, each = nT),
                             tt_obs = rep(1:nT, I)[!is.na(c(responses))],
                             ii_obs = rep(1:I, each = nT)[!is.na(c(responses))],
                             y_obs  = c(responses)[!is.na(c(responses))])
            
            begin.time <- proc.time()
            fit <- sampling(model,                            # Stan model. 
                            data = standata,                  # Data.
                            iter = 10000,                      # Number of iterations.
                            chains  = 3,                      # Number of chains.
                            warmup  = 1000,                   # Burn-in samples.
                            pars = c("beta", "theta", "lambda"),
                            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.
            run.time <- proc.time() - begin.time
            rm(begin.time)
            
            # Extract stan results ----
            sum.fit <- list()
            
            sum.fit$beta     <- summary(fit, pars = "beta")$summary
            sum.fit$theta    <- summary(fit, pars = "theta")$summary
            sum.fit$lambda   <- summary(fit, pars = "lambda")$summary
            
            #stan checks
            stan.diag <- monitor(extract(fit, permuted = FALSE, inc_warmup = FALSE), warmup = 0)
            
            ndiv  <- get_num_divergent(fit)     # number of divergent transitions
            nbfmi <- get_low_bfmi_chains(fit)   # number of chains with a low bayesian fraction missing information 
            ntree <- get_num_max_treedepth(fit) # number of transitions that exceeded the maximum treedepth
            nbulk <- sum(stan.diag$Bulk_ESS < 100 * dim(fit)[2]) # number of parameters with low bulk ESS
            ntail <- sum(stan.diag$Tail_ESS < 100 * dim(fit)[2]) # number of parameters with low tail ESS
          
            #max Rhat
            maxRhat <- round(max(rhat(fit, pars = c("beta", "theta", "lambda"))), 4)
            nRhat   <- sum(rhat(fit, pars = c("beta", "theta", "lambda")) > 1.05)
            
            
            # correlations theta alpha and beta
            beta.cor  <- round(cor(c(t(thresholds[, 1:M])), sum.fit$beta[, 1]), 4)
            theta.cor <- round(cor(theta, sum.fit$theta[, 1]), 4)
            
            # bias theta alpha beta lambda
            beta.bias   <- round(sum(sum.fit$beta[, 1] - c(t(thresholds[, 1:M]))) / (I * M), 4)
            theta.bias  <- round(sum(sum.fit$theta[, 1] - theta) / nT, 4)
            lambda.bias <- round(sum.fit$lambda[, 1] - lambda, 4)
            
            # abbias theta alpha beta lambda
            beta.abbias   <- round(sum(abs(sum.fit$beta[, 1] - c(t(thresholds[, 1:M])))) / (I * M), 4)
            theta.abbias  <- round(sum(abs(sum.fit$theta[, 1] - theta)) / nT, 4)
            lambda.abbias <- round(abs(sum.fit$lambda[, 1] - lambda), 4)
            
            # rmse theta alpha beta lambda
            beta.rmse  <- round(sqrt(sum((sum.fit$beta[, 1] - c(t(thresholds[, 1:M])))  ^ 2) / (I * M)), 4)
            theta.rmse <- round(sqrt(sum((sum.fit$theta[, 1] - theta) ^ 2) / nT), 4)
            
            # max CI theta alpha beta lambda
            beta.CI   <-  round(max(sum.fit$beta[, 8] - sum.fit$beta[, 4]), 4)
            theta.CI  <-  round(max(sum.fit$theta[, 8] - sum.fit$theta[, 4]), 4)
            lambda.CI <-  round(max(sum.fit$lambda[, 8] - sum.fit$lambda[, 4]), 4)
            
            # Conditional to check whether the mcmc algorithm is reliable (convergence and reasonable estimates)
            if (nRhat != 0   | ndiv != 0 | length(nbfmi) != 0) {
              corrupt <- 1
            } else {
              corrupt <- 0
            }
            
            # Conditional to check is more iterations or increasing the treedepth are needed.
            if (ntree != 0 | nbulk != 0 | ntail != 0) {
              efficiency <- 1
            } else {
              efficiency <- 0
            }
            
            # Export output to outer foreach ----
            result <- list(c(cond,
                             r,
                             round(run.time[3]),
                             ndiv,
                             length(nbfmi),
                             ntree,
                             nbulk,
                             ntail,
                             maxRhat,
                             nRhat,
                             beta.cor,
                             beta.bias,
                             beta.abbias,
                             beta.rmse,
                             beta.CI,
                             theta.cor,
                             theta.bias,
                             theta.abbias,
                             theta.rmse,
                             theta.CI,
                             lambda.bias,
                             lambda.abbias,
                             lambda.CI, 
                             corrupt,
                             efficiency))
            rm(standata, fit, sum.fit, thresholds, delta, taus, theta, responses,
               ndiv, nbfmi, ntree, nbulk, ntail, stan.diag,
               maxRhat, nRhat, beta.cor, beta.bias, beta.abbias,
               beta.rmse, beta.CI, theta.cor, theta.bias, theta.abbias,
               theta.rmse, theta.CI, lambda.bias, lambda.abbias,
               lambda.CI, corrupt, efficiency)
            result
          }

# 3.0 Export output ----

for (i in args[1]:args[2]) {
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[1]]
  colnames(tmp) <- c("cond",
                     "r",
                     "run.time",
                     "ndiv",
                     "nbfmi",
                     "ntree",
                     "nbulk",
                     "ntail",
                     "maxRhat",
                     "nRhat",
                     "beta.cor",
                     "beta.bias",
                     "beta.abbias",
                     "beta.rmse",
                     "beta.CI",
                     "theta.cor",
                     "theta.bias",
                     "theta.abbias",
                     "theta.rmse",
                     "theta.CI",
                     "lambda.bias",
                     "lambda.abbias",
                     "lambda.CI",
                     "corrupt",
                     "efficiency")
  write.table(tmp, file = paste0(getwd(), "/Simulation/Sim_AR_PCM_cond_", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
}

# Stop parallel cluster
stopCluster(cl)
rm(cl)

# END
