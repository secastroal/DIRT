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
clus <- makeCluster(ncores, outfile = "")
registerDoParallel(clus, cores = ncores)

# Get conditions and replications from the batch file
args <- commandArgs(trailingOnly = TRUE)

outcome.simulation <- foreach(cond = args[1]:args[2], .combine = 'list', .multicombine = TRUE) %:%
  foreach(r = args[3]:args[4], .combine = 'comb', .multicombine = TRUE, 
          .packages = c("rstan", "bayesplot", "scatterpie", "abind", 
                        "parallel", "doParallel"), 
          .export = c("acomb")
          ) %dopar% {
            
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
              # Generate two correlated factors
              thetabi <- matrix(NA, nrow = nT, ncol = 2)
              mubi    <- replicate(2, logarithmic(nT, maxAbsValue = 0.5))
              Sigma   <- matrix(c(1, par.model, par.model, 1), 2) 
              
              thetabi[1, ] <- MASS::mvrnorm(1, 
                                            mu = mubi[1, ]/(1 - gen.Data$lambda.gen), 
                                            Sigma = Sigma)
              
              for (i in 2:nT) {
                thetabi[i, ] <- mubi[i, ] + lambda * thetabi[i - 1, ] +
                  MASS::mvrnorm(1, mu = rep(0, 2), Sigma = Sigma)
              }
              rm(i, Sigma)
              
              
              gen.factor1 <- gen.TVDPCM(nT = nT,
                                        I  = ceiling(I / 2),
                                        K  = K,
                                        pop.param = list(
                                          lambda     = lambda,
                                          thresholds = gen.Data$thresholds.gen[1:ceiling(I / 2), ],
                                          tv_int     = mubi[, 1],
                                          theta      = thetabi[, 1]),
                                        seed = seed,
                                        FUN  = "logarithmic",
                                        maxAbsValue = 0.5)
              
              gen.factor2 <- gen.TVDPCM(nT = nT,
                                        I  = I - ceiling(I / 2),
                                        K  = K,
                                        pop.param = list(
                                          lambda     = lambda,
                                          thresholds = gen.Data$thresholds.gen[(ceiling(I / 2) + 1):I, ],
                                          tv_int     = mubi[, 2],
                                          theta      = thetabi[, 2]),
                                        seed = seed,
                                        FUN  = "logarithmic",
                                        maxAbsValue = 0.5)
              
              responses <- cbind(gen.factor1$data, gen.factor2$data)
              rm(mubi, gen.factor1, gen.factor2, thetabi)
            }
            
            if (gen.model == "GPCM") {
              # Discrimination parameters different from 1
              alpha.gpcm <- c(rep(1, I - I * par.model),
                              rep(c(.5, 1.5), each = I * par.model / 2))
              
              gen.tvgpcm <- gen.TVDPCM(nT = nT,
                                       I  = I,
                                       K  = K,
                                       pop.param = list(
                                         lambda     = lambda,
                                         thresholds = gen.Data$thresholds.gen,
                                         alpha      = alpha.gpcm,
                                         theta      = gen.Data$theta.gen),
                                       seed = seed,
                                       FUN  = "logarithmic",
                                       maxAbsValue = 0.5)
              responses <- gen.tvgpcm$data
              rm(alpha.gpcm, gen.tvgpcm)
            }
            
            if (gen.model == "DRIFT") {
              # Item Parameter Drift
              drift <- c(rep(0, I - I * par.model),
                         rep(c(-2, 2), each = I * par.model / 2))
              
              gen.drift1 <- gen.TVDPCM(
                nT = ceiling(2 * nT/3),
                I  = I,
                K  = K,
                pop.param = list(
                  lambda     = lambda,
                  thresholds = gen.Data$thresholds.gen,
                  theta      = gen.Data$theta.gen[1:ceiling(2 * nT/3)],
                  tv_int     = logarithmic(
                    nT, 
                    maxAbsValue = 0.5)[1:ceiling(2 * nT/3)]),
                seed = seed,
                FUN  = "logarithmic",
                maxAbsValue = 0.5)
              
              gen.drift2 <- gen.TVDPCM(
                nT = nT - ceiling(2 * nT/3),
                I  = I,
                K  = K,
                pop.param = list(
                  lambda     = lambda,
                  thresholds = gen.Data$thresholds.gen + drift,
                  theta      = gen.Data$theta.gen[(ceiling(2 * nT/3) + 1):nT],
                  tv_int     = logarithmic(
                    nT, 
                    maxAbsValue = 0.5)[(ceiling(2 * nT/3) + 1):nT]),
                seed = seed,
                FUN  = "logarithmic",
                maxAbsValue = 0.5)
              
              responses <- rbind(gen.drift1$data, gen.drift2$data)
              rm(drift, gen.drift1, gen.drift2)
            }
            
            if (gen.model == "Meaning") {
              responses <- gen.Data$data
              responses[ceiling(2*nT/3):nT, I] <- 6 - responses[ceiling(2*nT/3):nT, I]
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
            
            tot.time <- proc.time() - begin.time
            rm(begin.time)
            
            # Export output to outer foreach ----
            result <- list(c(cond,
                             r,
                             gen.model,
                             par.model,
                             round(run.time[3]),
                             round(tot.time[3]),
                             ndiv,
                             nbfmi,
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
            rm(standata, fit, gen.Data, responses, tot.time,
               run.time, cond, r, seed, gen.model, par.model,
               ndiv, nbfmi, ntree, nbulk, ntail, stan.diag, maxRhat, nRhat, 
               corrupt, efficiency, ppmc01, ppmc02, ppmc03, ppmc04, ppmc05,
               ppmc06, ppmc07, ppmc08, ppmc09, ppmc10, ppmc11, ppmc12,
               ppmc13)
            result
          }

# 3.0 Export output ----

for (i in args[1]:args[2]) {
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[1]]
  colnames(tmp) <- c("cond",
                     "r",
                     "genModel",
                     "parModel",
                     "run.time",
                     "tot.time",
                     "ndiv",
                     "nbfmi",
                     "ntree",
                     "nbulk",
                     "ntail",
                     "maxRhat",
                     "nRhat",
                     "corrupt",
                     "efficiency",
                     paste0("acf", 1:3),
                     "racf",
                     "lpacf",
                     "mssd",
                     paste0("itcor", 1:I),
                     paste0("q1_", 1:I),
                     paste0("q1alt_", 1:I),
                     paste0("ilpacf_", 1:I),
                     apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, 
                           function(x) paste0("q3_", paste(x, collapse = "_"))),
                     apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, 
                           function(x) paste0("or_", paste(x, collapse = "_"))),
                     apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, 
                           function(x) paste0("ordiff_", paste(x, collapse = "_"))),
                     apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, 
                           function(x) paste0("resid_", paste(x, collapse = "_"))),
                     apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, 
                           function(x) paste0("rediff_", paste(x, collapse = "_")))
                     )
  write.table(tmp, file = paste0(getwd(), "/Simulation/Sim_PPMC_cond_", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
}

# Stop parallel cluster
stopCluster(clus)
rm(clus)

# END ----