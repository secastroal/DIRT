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
options(mc.cores = parallel::detectCores())
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


# 2.0 Run Simulation ----

# 3.0 Export output ----



# END ----