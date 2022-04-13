# Utility functions to simulate tv-dpcm data and create datalist to export to stan.

# linear ----
# Create decreasing linear trend from the maxAbsValue to -maxAbsValue.
linear <- function(nT, maxAbsValue) {
  out <- seq(maxAbsValue, -maxAbsValue, length.out = nT)
  return(out)
}

# sinusoidal ----
# Create sinusoidal trend of one and a half cycle that starts at maxAbsValue and 
# ends at -maxAbsValue.
sinusoidal <- function(nT, maxAbsValue) {
  out <- maxAbsValue * cos(3 * pi * (1:nT)/(nT))
  return(out)
}

# logarithmic ----
# Create logistic trend based on the item characteristic function of dichotomous
# IRT models, which starts at -maxAbsValue and ends at maxAbsValue
logarithmic <-  function(nT, maxAbsValue) {
  out <- 2 * maxAbsValue * 
    (exp(0.05 * ((1:nT) - nT/2))/( 1 + exp(0.05 * ((1:nT) - nT/2)))) - 
    maxAbsValue
  return(out)
}


# select_trend ----
# Function to select the trend for the latent time series
select_trend <- function(nT, FUN, ...) {
  FUN <- match.fun(FUN)
  out <- forceAndCall(1, FUN, nT, ...)
  return(out)
}

# gen.TVDPCM
# Function to simulate data for the TV-DPCM
# The required arguments are the number of time points nT; the number of items I;
# the number of response categories K; the population parameters pop.param,
# which must be given as a list with the true autoregressive effect (lambda), 
# the variance of the innovations (sigma2), the thresholds parameters 
# (thresholds), the discrimination parameters (alpha), and the theta parameters 
# (theta); the seed used to simulate the data; a function FUN which can be 
# predefined or customized as in apply; and additional arguments for FUN.
# The list of true parameters can be ommited or partially given. If missing,
# the autoregressive effect is randomly sampled from a uniform distribution,
# the variance of the innovations is set to 1, the threshold parameters 
# are randomly generated, the discrimination parameters are set to 1, and the 
# theta parameters are generated based on the TV-DPCM and the given function for 
# the trend of the intercept. The function must create a vector of length nT. 
# Even when theta is provided in pop.param, the trend of the intercept must be 
# specified via FUN. To have the intercept constant over time define 
# FUN = function(x) {rep(0, nT)}.
gen.TVDPCM <-  function (nT, I, K, pop.param = NULL, seed, FUN, ...) {
  # Set seed to allow replicability
  set.seed(seed)
  
  # Get the time indicator and the number of thresholds per item
  time <- 1:nT
  M    <- K - 1
  
  # Check population parameters and create them if missing.
  if (is.null(pop.param)) {pop.param <- list()}
  
  if (is.null(pop.param$lambda)) {
    pop.param$lambda <- runif(1, -1, 1)
  }
  
  if (is.null(pop.param$sigma2)) {
    pop.param$sigma2 <- 1
  }
  
  if (is.null(pop.param$thresholds)) {
    # Create threshold parameters
    thresholds <- t(apply(matrix(runif(I * M, .3, 1), I), 1, cumsum))
    if (M == 1) {thresholds <- t(thresholds)}
    thresholds <- thresholds - rowMeans(thresholds)
    thresholds <- thresholds + rnorm(I)
    
    pop.param$thresholds <- thresholds
  }
  
  if (is.null(pop.param$alpha)) {
    pop.param$alpha <- rep(1, I)
  }
  
  # Generate the latent dynamic process if theta is missing ----
  
  # Create time varying intercept
  tv_int <- select_trend(nT, FUN, ...)
  
  if (is.null(pop.param$theta)) {
    # Create latent state disposition scores
    # First theta based on an stationary marginal distribution 
    # (see Bringmann et al., 2017).
    theta <- rep(NA, nT)
    
    theta[1] <- rnorm(1, tv_int[1]/(1 - pop.param$lambda), 
                      sqrt(pop.param$sigma2/(1 - pop.param$lambda ^ 2)))
    
    for (t in 2:nT) {
      theta[t] <- tv_int[t] + pop.param$lambda * theta[t - 1] + 
        rnorm(1, 0, sqrt(pop.param$sigma2))
    }
    rm(t)
    
    pop.param$theta <- theta
  }
  
  # Compute the attractor and the variance of the dynamic process
  attractor <- tv_int / (1 - pop.param$lambda)
  pvar      <- pop.param$sigma2 / (1 - pop.param$lambda ^ 2)
  
  # Generate item responses given the PCM ----
  # The function used to generate the responses uses another parameterization 
  # of the PCM, which splits the thresholds into location and step parameters. 
  
  # Location
  delta <- rowMeans(pop.param$thresholds)
  
  # Steps
  taus  <- pop.param$thresholds - delta
  
  # Generate responses
  probs.array <- array(NA, dim = c(length(pop.param$theta), I, K))
  
  for (y in 0:M) {
    probs.array[, , y + 1] <- P.GPCM(y     = y, 
                                     alpha = pop.param$alpha, 
                                     delta = delta, 
                                     taus  = taus, 
                                     theta = pop.param$theta, 
                                     M     = M)
  }
  # The reponses are coded from 1 to I
  responses   <- apply(probs.array, 1:2, function(vec) {
    which( rmultinom(1, 1, vec) == 1)
    })
  rm(probs.array, y)
  
  # Return list with the simulated data and the true population parameters
  out <- list(thresholds.gen = pop.param$thresholds,
              alpha.gen      = pop.param$alpha,
              theta.gen      = pop.param$theta,
              attractor.gen  = attractor,
              lambda.gen     = pop.param$lambda,
              sigma2.gen     = pop.param$sigma2,
              pvar.gen       = pvar,
              data           = responses)
  
  return(out)
}



