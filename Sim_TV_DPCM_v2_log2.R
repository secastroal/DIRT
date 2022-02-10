# Simulation testing Dynamic IRT models
# Testing TV-DPCM with time-varying intercept

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
library(splines)

# Load required functions
source("R/IRT_models.R")
#source("R/IRT_plots.R")

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

K        <- 5   # Number of categories per item
M        <- K - 1
in_var   <- 1   # Variance of the innovations
n_knots  <- 4  # Define number of knots.
s_degree <- 3  # Define degree of splines.
n_basis  <- n_knots + s_degree - 1 # Compute number of basis.

# Manipulated conditions

N.timepoints <- c(200, 300)           # Number of timepoints
N.items      <- c(6)                  # Number of items 
S.lambda     <- c(0.25, 0.5)          # Size of the autoregressive effect
M.prop       <- c(0)                  # Proportion of missing values

Cond        <- expand.grid(N.timepoints, N.items, S.lambda, M.prop)
names(Cond) <- c("nT", "I", "lambda", "NAprop")

rm(N.timepoints, N.items, S.lambda, M.prop)

# Compile the model
# To run the simulation we commented out the generated quantities block to 
# reduce memory usage.
model <- stan_model(file = "Stan/tv_dpcm_int_v2.stan", verbose = FALSE)

# 2.0 Run simulation ----

# Setup parallel backend to use parallel tasks
cl <- makeCluster(24)
registerDoParallel(cl, cores = 24)

# Get conditions and replications from the batch file
args <- commandArgs(trailingOnly = TRUE)

# Nested foreach loop (across conditions and across replications)

outcome.simulation <- foreach(cond = args[1]:args[2], .combine = 'list', .multicombine = TRUE) %:%
  foreach(r = args[3]:args[4], .combine = 'comb', .multicombine = TRUE, 
          .packages = c("rstan", "bayesplot", "splines")) %dopar% {
            
            # Define manipulated factors and seed
            nT     <- Cond[cond, 1]
            I      <- Cond[cond, 2]
            lambda <- Cond[cond, 3]
            naprop <- Cond[cond, 4]
            seed   <- 1000 * cond + r
            
            # Simulate data ----
            set.seed(seed)
            
            # True states
            time <- 1:nT
            
            # Create time varying intercept with s-shape growth
            tv_int <- 4 * (exp(0.05 * (x - nT/2))/( 1 + exp(0.05 * (x - nT/2)))) - 2
            
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
            
            # Introduce missing values per "beep" (random rows with complete NAs)
            if (naprop != 0) {
              naind <- rbinom(dim(responses)[1], 1, naprop)
              responses[naind == 1, ] <- NA
              rm(naind)
            }
                        
            # Define twice the knots for fitting
            knots <- trunc(seq(1, nT, length.out = n_knots * 2))
            
            # Fit the TV-DPCM model with stan ----
            standata <- list(nT = nT,
                             n_knots  = n_knots * 2,
                             knots    = knots,
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
            
            tvdpcm_inits <- function() {
              list(#a_raw  = rnorm(n_basis),
                   #a0     = rnorm(1),
                   #tau    = rnorm(1),
                   lambda = runif(1, -1, 1),
                   beta   = array(rnorm(I * (K - 1), 0, 3), dim = c(I, K - 1)),
                   inno   = rnorm(nT, 0, 3),
                   sigma2 = rlnorm(1, 1))
            }

            begin.time <- proc.time()
            fit <- sampling(model,                            # Stan model. 
                            data = standata,                  # Data.
                            iter = 2000,                      # Number of iterations.
                            chains  = 3,                      # Number of chains.
                            warmup  = 1000,                   # Burn-in samples.
                            init    = tvdpcm_inits,
                            seed = seed,
                            pars = c("beta", "theta", "lambda",
                                     "sigma2", "p_var", "attractor"),
                            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.
            run.time <- proc.time() - begin.time
            rm(begin.time)
            
            # Extract stan results ----
            sum.fit <- list()
            
            sum.fit$beta     <- summary(fit, pars = "beta")$summary
            sum.fit$theta    <- summary(fit, pars = "theta")$summary
            sum.fit$lambda   <- summary(fit, pars = "lambda")$summary
            sum.fit$sigma2 <- summary(fit, pars = "sigma2")$summary
            sum.fit$p_var  <- summary(fit, pars = "p_var")$summary
            sum.fit$attractor <- summary(fit, pars = "attractor")$summary
            
            #stan checks
            stan.diag <- monitor(extract(fit, permuted = FALSE, inc_warmup = FALSE), warmup = 0)
            
            ndiv  <- get_num_divergent(fit)     # number of divergent transitions
            nbfmi <- get_low_bfmi_chains(fit)   # number of chains with a low bayesian fraction missing information 
            ntree <- get_num_max_treedepth(fit) # number of transitions that exceeded the maximum treedepth
            nbulk <- sum(stan.diag$Bulk_ESS < 100 * dim(fit)[2]) # number of parameters with low bulk ESS
            ntail <- sum(stan.diag$Tail_ESS < 100 * dim(fit)[2]) # number of parameters with low tail ESS
          
            #max Rhat
            maxRhat <- round(max(rhat(fit, pars = c("beta", "theta", "lambda", 
                                                    "sigma2", "p_var", "attractor"))), 4)
            nRhat   <- sum(rhat(fit, pars = c("beta", "theta", "lambda",
                                              "sigma2", "p_var", "attractor")) > 1.05)
            
            # correlations
            beta.cor  <- round(cor(c(t(thresholds)), sum.fit$beta[, 1]), 4)
            theta.cor <- round(cor(theta, sum.fit$theta[, 1]), 4)
            attra.cor <- round(cor(attractor, sum.fit$attractor[, 1]), 4)
            
            # bias
            beta.bias   <- round(sum(sum.fit$beta[, 1] - c(t(thresholds))) / (I * M), 4)
            theta.bias  <- round(sum(sum.fit$theta[, 1] - theta) / nT, 4)
            attra.bias  <- round(sum(sum.fit$attractor[, 1] - attractor) / nT, 4)
            pvar.bias   <- round(sum.fit$p_var[, 1] - p_var, 4)
            lambda.bias <- round(sum.fit$lambda[, 1] - lambda, 4)
            sigma2.bias <- round(sum.fit$sigma2[, 1] - in_var, 4)
            
            # abbias
            beta.abbias   <- round(sum(abs(sum.fit$beta[, 1] - c(t(thresholds)))) / (I * M), 4)
            theta.abbias  <- round(sum(abs(sum.fit$theta[, 1] - theta)) / nT, 4)
            attra.abbias  <- round(sum(abs(sum.fit$attractor[, 1] - attractor)) / nT, 4)
            pvar.abbias   <- round(abs(sum.fit$p_var[, 1] - p_var), 4)
            lambda.abbias <- round(abs(sum.fit$lambda[, 1] - lambda), 4)
            sigma2.abbias <- round(abs(sum.fit$sigma2[, 1] - in_var), 4)
            
            # rmse
            beta.rmse  <- round(sqrt(sum((sum.fit$beta[, 1] - c(t(thresholds)))  ^ 2) / (I * M)), 4)
            theta.rmse <- round(sqrt(sum((sum.fit$theta[, 1] - theta) ^ 2) / nT), 4)
            attra.rmse <- round(sqrt(sum((sum.fit$attractor[, 1] - attractor) ^ 2) / nT), 4)
            
            # relative bias
            beta.rbias   <- round(sum((sum.fit$beta[, 1] - c(t(thresholds))) / c(t(thresholds))) / (I * M), 4)
            theta.rbias  <- round(sum((sum.fit$theta[, 1] - theta) / theta) / nT, 4)
            attra.rbias  <- round(sum((sum.fit$attractor[, 1] - attractor) / attractor) / nT, 4)
            pvar.rbias   <- round((sum.fit$p_var[, 1] - p_var) / p_var, 4)
            lambda.rbias <- round((sum.fit$lambda[, 1] - lambda) / lambda, 4)
            sigma2.rbias <- round((sum.fit$sigma2[, 1] - in_var) / in_var, 4)
            
            # relative bias (Schultzberg & Muthen, 2018)
            beta.rbias2   <- round(sum(sum.fit$beta[, 1] / c(t(thresholds))) / (I * M), 4)
            theta.rbias2  <- round(sum(sum.fit$theta[, 1] / theta) / nT, 4)
            attra.rbias2  <- round(sum(sum.fit$attractor[, 1] / attractor) / nT, 4)
            pvar.rbias2   <- round(sum.fit$p_var[, 1] / p_var, 4)
            lambda.rbias2 <- round(sum.fit$lambda[, 1] / lambda, 4)
            sigma2.rbias2 <- round(sum.fit$sigma2[, 1] / in_var, 4)
            
            # coverage
            beta.cover   <-  round(sum(sum.fit$beta[, 4] <= c(t(thresholds)) & 
                                   c(t(thresholds)) <= sum.fit$beta[, 8]) / (I * M), 4)
            theta.cover  <-  round(sum(sum.fit$theta[, 4] <= theta &
                                   theta <= sum.fit$theta[, 8]) / nT, 4)
            attra.cover  <-  round(sum(sum.fit$attractor[, 4] <= attractor &
                                   attractor <= sum.fit$attractor[, 8]) / nT, 4)
            pvar.cover   <-  round(sum(sum.fit$p_var[, 4] <= p_var &
                                   p_var <= sum.fit$p_var[, 8]), 4)
            lambda.cover <-  round(sum(sum.fit$lambda[, 4] <= lambda & 
                                   lambda <= sum.fit$lambda[, 8]), 4)
            sigma2.cover <-  round(sum(sum.fit$sigma2[, 4] <= in_var & 
                                   in_var <= sum.fit$sigma2[, 8]), 4)
            
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
                             beta.rbias,
                             beta.rbias2,
                             beta.rmse,
                             beta.cover,
                             theta.cor,
                             theta.bias,
                             theta.abbias,
                             theta.rbias,
                             theta.rbias2,
                             theta.rmse,
                             theta.cover,
                             attra.cor,
                             attra.bias,
                             attra.abbias,
                             attra.rbias,
                             attra.rbias2,
                             attra.rmse,
                             attra.cover,
                             lambda.bias,
                             lambda.abbias,
                             lambda.rbias,
                             lambda.rbias2,
                             lambda.cover, 
                             pvar.bias,
                             pvar.abbias,
                             pvar.rbias,
                             pvar.rbias2,
                             pvar.cover, 
                             sigma2.bias,
                             sigma2.abbias,
                             sigma2.rbias,
                             sigma2.rbias2,
                             sigma2.cover, 
                             corrupt,
                             efficiency))
            rm(standata, fit, sum.fit, thresholds, delta, taus, theta, responses,
               B_true, a, a0, attractor, knots, tv_lambda, tv_int, time,
               run.time, p_var, nT, naprop, lambda, I, cond, r, seed,
               ndiv, nbfmi, ntree, nbulk, ntail, stan.diag, maxRhat, nRhat, 
               beta.cor, beta.bias, beta.abbias, beta.rmse, beta.rbias, beta.rbias2,
               beta.cover, theta.cor, theta.bias, theta.abbias, theta.rbias,
               theta.rmse, theta.rbias2, theta.cover, attra.abbias, attra.bias,
               attra.cor, attra.cover, attra.rbias, attra.rbias2, attra.rmse,  
               lambda.bias, lambda.abbias, lambda.cover, lambda.rbias, lambda.rbias2,
               sigma2.abbias, sigma2.bias, sigma2.cover, sigma2.rbias, sigma2.rbias2,
               pvar.abbias, pvar.bias, pvar.cover, pvar.rbias, pvar.rbias2,
               corrupt, efficiency)
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
                     "beta.rbias",
                     "beta.rbias2",
                     "beta.rmse",
                     "beta.cover",
                     "theta.cor",
                     "theta.bias",
                     "theta.abbias",
                     "theta.rbias",
                     "theta.rbias2",
                     "theta.rmse",
                     "theta.cover",
                     "attra.cor",
                     "attra.bias",
                     "attra.abbias",
                     "attra.rbias",
                     "attra.rbias2",
                     "attra.rmse",
                     "attra.cover",
                     "lambda.bias",
                     "lambda.abbias",
                     "lambda.rbias",
                     "lambda.rbias2",
                     "lambda.cover", 
                     "pvar.bias",
                     "pvar.abbias",
                     "pvar.rbias",
                     "pvar.rbias2",
                     "pvar.cover", 
                     "sigma2.bias",
                     "sigma2.abbias",
                     "sigma2.rbias",
                     "sigma2.rbias2",
                     "sigma2.cover",
                     "corrupt",
                     "efficiency")
  write.table(tmp, file = paste0(getwd(), "/Simulation/Sim_TV_DPCM_v2_log2_cond_", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
}

# Stop parallel cluster
stopCluster(cl)
rm(cl)

# END
