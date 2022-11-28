# Simulation study to test the PPMC methods for the TV-DPCM

# Data generating models

# TV-DPCM (control)
# TV-MDPCM (maybe low vs high correlation)
# TV-GPCM (1/3 or 2/3 different)
# TV-DPCM-IPD (1/3 or 2/3 different)
# TV-DPCM-Meaning (1/3 or 2/3 different)


# Conditions
# Time series length nT <- 300
# Either 3 or 6 items (preferably 6)

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
options(mc.cores = 1)
suppressPackageStartupMessages(library(bayesplot))

# Load required functions
source("R/IRT_models.R")
source("R/tvdpcm2stan.R")
source("R/genTVDPCM.R")
source("R/PPMC.R")

# Create folder to save output
if (!dir.exists(paste(getwd(), "Simulation", sep = "/"))) {
  dir.create(paste(getwd(), "Simulation", sep = "/"), recursive = TRUE)
}

# Create function needed to combine foreach output ----
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# 1.0 Set up Conditions ----

# Number of parallel analyses.
ncores <- 24

# Fix conditions

nT     <- 300   # Number of time points
I      <- 6     # Number of items
K      <- 5     # Number of categories per item
M      <- K - 1 # Number of thresholds per item
lambda <- 0.5   # Size of the autoregressive effect
in_var <- 1     # Variance of the innovations

# Manipulated conditions
# Generating model.
#   'TV-DPCM': No violation to assumptions of TV-DPCM
#   'BiDim': Unidimensionality is violated.
#   'GPCM': Discrimination is not equal for all items.
#   'Drift': Longitudinal measurement invariance does not hold. Parameter Drift.
#   'Meaning': The meaning of the items changes. 
#              This is also measurement non-invariance
gen.model <- c("TV-DPCM", "BiDim", "BiDim", "GPCM", "GPCM", 
               "DRIFT", "DRIFT", "Meaning")
# Parameters of the generating model. Correlation between latent variables for 
# 'BiDim' and proportion of items with discrimination different from 1 or with
# item parameter drift for 'GPCM' and 'DRIFT' respectively.
par.model <- c(NA, 0.3, 0.6, 1/3, 2/3, 1/3, 2/3, NA)

# Matrix of conditions to loop through
Cond <- data.frame(gen.model, par.model)
names(Cond) <- c("genModel", "parModel")
rm(gen.model, par.model)

# Compile the model
# The generated quantities block is required to compute the PPMC methods.
model <- stan_model(file = "Stan/tv_dpcm_int_v5.1.stan", verbose = FALSE)

# 2.0 Run Simulation ----

# Setup parallel backend to use parallel tasks
cl <- makeCluster(ncores, outfile = "")
registerDoParallel(cl, cores = ncores)

# Get conditions and replications from the batch file
args <- commandArgs(trailingOnly = TRUE)

# 3.0 Export output ----

outcome.simulation <- foreach(cond = args[1]:args[2], .combine = 'list', .multicombine = TRUE) %:%
  foreach(r = args[3]:args[4], .combine = 'comb', .multicombine = TRUE, 
          .packages = c("rstan", "bayesplot", "scatterpie", "abind"), 
          .export = c("logarithmic")) %dopar% {
            
            # Define manipulated factors and seed
            gen.model <- Cond[cond, 1]
            par.model <- Cond[cond, 2]
            seed      <- 1000 + r
            
            cat(sprintf("This is the %dth replication of condition %d\n", r, cond))
            
            # Simulate data ----
            
            gen.Data <- gen.TVDPCM(nT = nT,
                                   I  = I,
                                   K  = K,
                                   pop.param = list(lambda = lambda),
                                   seed = seed,
                                   FUN  = "logarithmic",
                                   maxAbsValue = 0.5)
            
            # Modify data given gen.model
            
            if (gen.model == "TV-DPCM") {
              responses <- gen.Data$data
            }
            
            if (gen.model == "BiDim") {
              
            }
            
            if (gen.model == "GPCM") {
              
            }
            
            if (gen.model == "DRIFT") {
              
            }
            
            if (gen.model == "Meaning") {
              
            }
            
            # Fit the TV-DPCM model with stan ----
            standata <- tvdpcm2stan_data(resp = as.matrix(responses),
                                         I    = I,
                                         K    = K,
                                         nT   = nT,
                                         n_knots  = 8,
                                         s_degree = 3)
            
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
                                     "sigma2", "pvar", "attractor", "rep_y"),
                            control = list(adapt_delta   = 0.99,
                                           max_treedepth = 15), # Other parameters to control sampling behavior.
                            show_messages = FALSE,
                            refresh = 0) 
            run.time <- proc.time() - begin.time
            rm(begin.time)
            
            # Stan diagnostic checks ----
            
            # All these values are expected to be 0.
            stan.diag <- monitor(
              extract(fit, 
                      pars = c("beta", "theta", "lambda",
                               "sigma2", "pvar", "attractor"), 
                      permuted = FALSE, inc_warmup = FALSE),
              warmup = 0, 
              print = FALSE)
            
            ndiv  <- get_num_divergent(fit)     # number of divergent transitions
            nbfmi <- get_low_bfmi_chains(fit)   # number of chains with a low bayesian fraction missing information 
            ntree <- get_num_max_treedepth(fit) # number of transitions that exceeded the maximum treedepth
            nbulk <- sum(stan.diag$Bulk_ESS < 100 * dim(fit)[2]) # number of parameters with low bulk ESS
            ntail <- sum(stan.diag$Tail_ESS < 100 * dim(fit)[2]) # number of parameters with low tail ESS
            
            #max Rhat
            maxRhat <- round(max(rhat(fit,
                                      pars = c("beta", "theta", 
                                               "lambda", "sigma2", 
                                               "pvar", "attractor"))), 4)
            nRhat   <- sum(rhat(fit, 
                                pars = c("beta", "theta", 
                                         "lambda", "sigma2", 
                                         "pvar", "attractor")) > 1.05)
            # nbfmi as number
            nbfmi <- length(nbfmi)
            
            # Conditional to check whether the mcmc algorithm is reliable 
            # (convergence and reasonable estimates)
            if (nRhat != 0   | ndiv != 0 | nbfmi != 0) {
              corrupt <- 1
            } else {
              corrupt <- 0
            }
            
            # Conditional to check is more iterations or increasing the 
            # treedepth are needed.
            if (ntree != 0 | nbulk != 0 | ntail != 0) {
              efficiency <- 1
            } else {
              efficiency <- 0
            }
            
            # Compute PPMC methods ----
            
            ppmc01 <- ppmc.acf(object = fit, data = standata, lag.max = 3)
            ppmc02 <- ppmc.racf(object = fit, data = standata)
            ppmc03 <- ppmc.lpacf(object = fit, data = standata, quiet = TRUE, 
                                 sumscores = TRUE)
            ppmc04 <- ppmc.mssd(object = fit, data = standata)
            ppmc05 <- ppmc.itcor(object = fit, data = standata, quiet = TRUE, 
                                 method = "pearson")
            ppmc06 <- ppmc.Q1(object = fit, data = standata, quiet = TRUE)
            ppmc07 <- ppmc.Q1.alt(object = fit, data = standata, quiet = TRUE)
            ppmc08 <- ppmc.lpacf(object = fit, data = standata, quiet = TRUE)
            ppmc09 <- ppmc.Q3(object = fit, data = standata)$ppp
            ppmc10 <- ppmc.OR(object = fit, data = standata)$ppp
            ppmc11 <- ppmc.ORDiff(object = fit, data = standata)$ppp
            ppmc12 <- ppmc.cov.resid(object = fit, data = standata)$ppp
            ppmc13 <- ppmc.cov.rediff(object = fit, data = standata)$ppp
            
            # Export output to outer foreach ----
            result <- list(c(cond,
                             r,
                             gen.model,
                             par.model,
                             round(run.time[3]),
                             ndiv,
                             length(nbfmi),
                             ntree,
                             nbulk,
                             ntail,
                             maxRhat,
                             nRhat,
                             corrupt,
                             efficiency,
                             ppmc01,
                             ppmc02,
                             ppmc03,
                             ppmc04,
                             ppmc05,
                             ppmc06,
                             ppmc07,
                             ppmc08,
                             ppmc09,
                             ppmc10,
                             ppmc11,
                             ppmc12,
                             ppmc13))
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

# END ----