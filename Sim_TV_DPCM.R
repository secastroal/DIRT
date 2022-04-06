# Simulation testing Dynamic IRT models
# Testing TV-DPCM with time-varying intercept

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
source("R/tvdpcm2stan.R")
source("R/genTVDPCM.R")

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
n_knots  <- 8   # Define number of knots.
s_degree <- 3   # Define degree of splines.

# Manipulated conditions

N.timepoints <- c(100, 200, 300, 500) # Number of timepoints
N.items      <- c(3, 6)               # Number of items 
S.lambda     <- c(0, 0.25, 0.5)       # Size of the autoregressive effect
M.prop       <- c(0, 0.3)       # Proportion of missing values
T.trend      <- c("sinusoidal") # Shape of the true trend

Cond        <- expand.grid(N.timepoints, N.items, S.lambda, M.prop, T.trend)
names(Cond) <- c("nT", "I", "lambda", "NAprop", "Trend")

rm(N.timepoints, N.items, S.lambda, M.prop, T.trend)

# Compile the model
# To run the simulation we commented out the generated quantities block to 
# reduce memory usage.
model <- stan_model(file = "Stan/tv_dpcm_int_v5.1.stan", verbose = FALSE)

# 2.0 Run simulation ----

# Setup parallel backend to use parallel tasks
cl <- makeCluster(24, outfile = "")
registerDoParallel(cl, cores = 24)

# Get conditions and replications from the batch file
args <- commandArgs(trailingOnly = TRUE)

# Nested foreach loop (across conditions and across replications)

outcome.simulation <- foreach(cond = args[1]:args[2], .combine = 'list', .multicombine = TRUE) %:%
  foreach(r = args[3]:args[4], .combine = 'comb', .multicombine = TRUE, 
          .packages = c("rstan", "bayesplot"), 
          .export = c("linear", "sinusoidal", "logarithmic")) %dopar% {
            
            # Define manipulated factors and seed
            nT     <- Cond[cond, 1]
            I      <- Cond[cond, 2]
            lambda <- Cond[cond, 3]
            naprop <- Cond[cond, 4]
            trend  <- Cond[cond, 5]
            seed   <- 1000 * cond + r
            
            cat(sprintf("This is the %dth replication of condition %d\n", r, cond))
            
            # Simulate data ----
            
            gen.Data <- gen.TVDPCM(nT = nT,
                                   I  = I,
                                   K  = K,
                                   pop.param = list(lambda = lambda),
                                   seed = seed,
                                   FUN  = as.character(trend),
                                   maxAbsValue = 1)
            
            # Introduce missing values per "beep" (random rows with complete NAs)
            if (naprop != 0) {
              naind    <- rbinom(nT, 1, naprop)
              naind[1] <- 0 # The first time point is always observed.
              gen.Data$data[naind == 1, ] <- NA
              rm(naind)
            }
                        
            # Fit the TV-DPCM model with stan ----
            standata <- tvdpcm2stan_data(resp = gen.Data$data,
                                         I    = I,
                                         K    = K,
                                         nT   = nT,
                                         n_knots  = n_knots,
                                         s_degree = s_degree)
            
            tvdpcm_inits <- function() {
              list(lambda = runif(1, -1, 1),
                   beta   = array(rnorm(I * (K - 1), 0, 3), dim = c(I, K - 1)),
                   inno   = rnorm(nT, 0, 3),
                   sigma  = rlnorm(1, 1))
            }

            begin.time <- proc.time()
            fit <- sampling(model,                  # Stan model. 
                            data = standata,        # Data.
                            iter = 2000,            # Number of iterations.
                            chains  = 3,            # Number of chains.
                            warmup  = 500,          # Burn-in samples.
                            init    = tvdpcm_inits, # Initial values
                            seed = seed,            # Seed
                            pars = c("beta", "theta", "lambda",
                                     "sigma2", "pvar", "attractor"),
                            control = list(adapt_delta   = 0.99,
                                           max_treedepth = 15), # Other parameters to control sampling behavior.
                            show_messages = FALSE,
                            refresh = 0) 
            run.time <- proc.time() - begin.time
            rm(begin.time)
            
            # Extract stan results ----
            # Create list to store the summaries of the parameters of interest.
            sum.fit <- list()
            
            sum.fit$beta   <- summary(fit, pars = "beta")$summary
            sum.fit$theta  <- summary(fit, pars = "theta")$summary
            sum.fit$lambda <- summary(fit, pars = "lambda")$summary
            sum.fit$sigma2 <- summary(fit, pars = "sigma2")$summary
            sum.fit$pvar   <- summary(fit, pars = "pvar")$summary
            sum.fit$attractor <- summary(fit, pars = "attractor")$summary
            
            # Get stan convergence checks. 
            # All these values are expected to be 0.
            stan.diag <- monitor(extract(fit, permuted = FALSE, inc_warmup = FALSE), 
                                 warmup = 0, print = FALSE)
            
            ndiv  <- get_num_divergent(fit)     # number of divergent transitions
            nbfmi <- get_low_bfmi_chains(fit)   # number of chains with a low bayesian fraction missing information 
            ntree <- get_num_max_treedepth(fit) # number of transitions that exceeded the maximum treedepth
            nbulk <- sum(stan.diag$Bulk_ESS < 100 * dim(fit)[2]) # number of parameters with low bulk ESS
            ntail <- sum(stan.diag$Tail_ESS < 100 * dim(fit)[2]) # number of parameters with low tail ESS
          
            #max Rhat
            maxRhat <- round(max(rhat(fit, pars = c("beta", "theta", "lambda", 
                                                    "sigma2", "pvar", "attractor"))), 4)
            nRhat   <- sum(rhat(fit, pars = c("beta", "theta", "lambda",
                                              "sigma2", "pvar", "attractor")) > 1.05)
            
            # Compare the true and the estimated parameters. 
            # correlations 
            beta.cor  <- round(cor(c(t(gen.Data$thresholds.gen)), sum.fit$beta[, 1]), 4)
            theta.cor <- round(cor(gen.Data$theta.gen, sum.fit$theta[, 1]), 4)
            attra.cor <- round(cor(gen.Data$attractor.gen, sum.fit$attractor[, 1]), 4)
            
            # bias
            beta.bias   <- round(sum(sum.fit$beta[, 1] - c(t(gen.Data$thresholds.gen))) / (I * M), 4)
            theta.bias  <- round(sum(sum.fit$theta[, 1] - gen.Data$theta.gen) / nT, 4)
            attra.bias  <- round(sum(sum.fit$attractor[, 1] - gen.Data$attractor.gen) / nT, 4)
            pvar.bias   <- round(sum.fit$pvar[, 1] - gen.Data$pvar.gen, 4)
            lambda.bias <- round(sum.fit$lambda[, 1] - gen.Data$lambda.gen, 4)
            sigma2.bias <- round(sum.fit$sigma2[, 1] - gen.Data$sigma2.gen, 4)
            
            # abbias
            beta.abbias   <- round(sum(abs(sum.fit$beta[, 1] - c(t(gen.Data$thresholds.gen)))) / (I * M), 4)
            theta.abbias  <- round(sum(abs(sum.fit$theta[, 1] - gen.Data$theta.gen)) / nT, 4)
            attra.abbias  <- round(sum(abs(sum.fit$attractor[, 1] - gen.Data$attractor.gen)) / nT, 4)
            pvar.abbias   <- round(abs(sum.fit$pvar[, 1] - gen.Data$pvar.gen), 4)
            lambda.abbias <- round(abs(sum.fit$lambda[, 1] - gen.Data$lambda.gen), 4)
            sigma2.abbias <- round(abs(sum.fit$sigma2[, 1] - gen.Data$sigma2.gen), 4)
            
            # rmse
            beta.rmse  <- round(sqrt(sum((sum.fit$beta[, 1] - c(t(gen.Data$thresholds.gen)))  ^ 2) / (I * M)), 4)
            theta.rmse <- round(sqrt(sum((sum.fit$theta[, 1] - gen.Data$theta.gen) ^ 2) / nT), 4)
            attra.rmse <- round(sqrt(sum((sum.fit$attractor[, 1] - gen.Data$attractor.gen) ^ 2) / nT), 4)
            
            # relative bias (does not work when the true parameter tends to 0)
            # only computed for the variance parameters and the autoregressive 
            # effect when it is different from 0.
            pvar.rbias   <- round((sum.fit$pvar[, 1] - gen.Data$pvar.gen) / gen.Data$pvar.gen, 4)
            lambda.rbias <- ifelse(gen.Data$lambda.gen == 0,
                                   NA,
                                   round((sum.fit$lambda[, 1] - gen.Data$lambda.gen) / gen.Data$lambda.gen, 4))
            sigma2.rbias <- round((sum.fit$sigma2[, 1] - gen.Data$sigma2.gen) / gen.Data$sigma2.gen, 4)
            
            # relative bias (Schultzberg & Muthen, 2018)
            # only for variances and autoregressive effect
            pvar.rbias2   <- round(sum.fit$pvar[, 1] / gen.Data$pvar.gen, 4)
            lambda.rbias2 <- ifelse(gen.Data$lambda.gen == 0,
                                    NA,
                                    round(sum.fit$lambda[, 1] / gen.Data$lambda.gen, 4))
            sigma2.rbias2 <- round(sum.fit$sigma2[, 1] / gen.Data$sigma2.gen, 4)
            
            # coverage of the C.I.
            beta.cover   <-  round(sum(sum.fit$beta[, 4] <= c(t(gen.Data$thresholds.gen)) & 
                                   c(t(gen.Data$thresholds.gen)) <= sum.fit$beta[, 8]) / (I * M), 4)
            theta.cover  <-  round(sum(sum.fit$theta[, 4] <= gen.Data$theta.gen &
                                   gen.Data$theta.gen <= sum.fit$theta[, 8]) / nT, 4)
            attra.cover  <-  round(sum(sum.fit$attractor[, 4] <= gen.Data$attractor.gen &
                                   gen.Data$attractor.gen <= sum.fit$attractor[, 8]) / nT, 4)
            pvar.cover   <-  round(sum(sum.fit$pvar[, 4] <= gen.Data$pvar.gen &
                                   gen.Data$pvar.gen <= sum.fit$pvar[, 8]), 4)
            lambda.cover <-  round(sum(sum.fit$lambda[, 4] <= gen.Data$lambda.gen & 
                                   gen.Data$lambda.gen <= sum.fit$lambda[, 8]), 4)
            sigma2.cover <-  round(sum(sum.fit$sigma2[, 4] <= gen.Data$sigma2.gen & 
                                   gen.Data$sigma2.gen <= sum.fit$sigma2[, 8]), 4)
            
            # Average size of the C.I.
            beta.CI   <-  round(sum(sum.fit$beta[, 8] - sum.fit$beta[, 4]) / (I * M), 4)
            theta.CI  <-  round(sum(sum.fit$theta[, 8] - sum.fit$theta[, 4]) / nT, 4)
            attra.CI  <-  round(sum(sum.fit$attractor[, 8] - sum.fit$attractor[, 4]) / nT, 4)
            pvar.CI   <-  round(sum(sum.fit$pvar[, 8] - sum.fit$pvar[, 4]), 4)
            lambda.CI <-  round(sum(sum.fit$lambda[, 8] - sum.fit$lambda[, 4]), 4)
            sigma2.CI <-  round(sum(sum.fit$sigma2[, 8] - sum.fit$sigma2[, 4]), 4)
            
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
                             beta.cover,
                             beta.CI,
                             theta.cor,
                             theta.bias,
                             theta.abbias,
                             theta.rmse,
                             theta.cover,
                             theta.CI,
                             attra.cor,
                             attra.bias,
                             attra.abbias,
                             attra.rmse,
                             attra.cover,
                             attra.CI,
                             lambda.bias,
                             lambda.abbias,
                             lambda.rbias,
                             lambda.rbias2,
                             lambda.cover,
                             lambda.CI,
                             pvar.bias,
                             pvar.abbias,
                             pvar.rbias,
                             pvar.rbias2,
                             pvar.cover,
                             pvar.CI,
                             sigma2.bias,
                             sigma2.abbias,
                             sigma2.rbias,
                             sigma2.rbias2,
                             sigma2.cover,
                             sigma2.CI,
                             corrupt,
                             efficiency))
            rm(standata, fit, sum.fit, gen.Data, trend,
               run.time, nT, naprop, lambda, I, cond, r, seed,
               ndiv, nbfmi, ntree, nbulk, ntail, stan.diag, maxRhat, nRhat, 
               beta.cor, beta.bias, beta.abbias, beta.rmse, beta.cover,
               beta.CI, theta.cor, theta.bias, theta.abbias, theta.rmse,
               theta.cover, theta.CI, attra.cor, attra.bias, attra.abbias,
               attra.rmse, attra.cover, attra.CI, lambda.bias, lambda.abbias,
               lambda.rbias, lambda.rbias2, lambda.cover, lambda.CI,
               pvar.bias, pvar.abbias, pvar.rbias, pvar.rbias2, pvar.cover,
               pvar.CI, sigma2.bias, sigma2.abbias, sigma2.rbias, sigma2.rbias2,
               sigma2.cover, sigma2.CI, corrupt, efficiency)
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
                     "beta.cover",
                     "beta.CI",
                     "theta.cor",
                     "theta.bias",
                     "theta.abbias",
                     "theta.rmse",
                     "theta.cover",
                     "theta.CI",
                     "attra.cor",
                     "attra.bias",
                     "attra.abbias",
                     "attra.rmse",
                     "attra.cover",
                     "attra.CI",
                     "lambda.bias",
                     "lambda.abbias",
                     "lambda.rbias",
                     "lambda.rbias2",
                     "lambda.cover",
                     "lambda.CI",
                     "pvar.bias",
                     "pvar.abbias",
                     "pvar.rbias",
                     "pvar.rbias2",
                     "pvar.cover",
                     "pvar.CI",
                     "sigma2.bias",
                     "sigma2.abbias",
                     "sigma2.rbias",
                     "sigma2.rbias2",
                     "sigma2.cover",
                     "sigma2.CI",
                     "corrupt",
                     "efficiency")
  write.table(tmp, file = paste0(getwd(), "/Simulation/Sim_TV_DPCM_cond_", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
}

# Stop parallel cluster
stopCluster(cl)
rm(cl)

# END
