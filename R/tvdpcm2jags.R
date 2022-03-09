# function to create the data list for the tv-dpcm jags models with different 
# smooth functions

# The required arguments are: A matrix of responses resp, the number of items,
# the number of response categories, the time variable as numeric (integers), 
# the number of basis functions, and the type of smoother.
tvdpcm2jags_data <- function(resp, I, K, time, n_basis, bs = c("bs", "ps", "tp")){
  # Get the number of time points
  nT <- length(time)
  
  # We use some functions of the mgcv package to create the basis and the S 
  # matrices for the penalization.
  
  # create temporal data
  theta <- rowSums(resp)
  
  tmp <- data.frame(theta, time)
  
  # Write formula with smooth term
  gp <- interpret.gam(formula(theta ~ s(time, bs = bs, k = n_basis)))
  
  # replicate code form jagam to be able to use gam.setup
  # using default values of jagam and gam.setup
  #cl <- match.call() # call needed in gam object for update to work
  #mf <- match.call(expand.dots=FALSE)
  #mf <- match.call(model.frame, model.frame(formula = gp$fake.formula, data = tmp))
  mf <- call("model.frame")
  mf$formula <- gp$fake.formula 
  mf$family <- mf$knots <- mf$sp <- mf$file <- mf$control <- 
    mf$centred <- mf$sp.prior <- mf$diagonalize <- NULL
  mf$drop.unused.levels <- TRUE
  #mf[[1]] <- quote(stats::model.frame) ##as.name("model.frame")
  pmf <- mf
  
  pmf$formula <- gp$pf
  pmf <- eval(pmf, tmp) # pmf contains all data for parametric part
  pterms <- attr(pmf,"terms") ## pmf only used for this
  rm(pmf, mf)
  
  # Use gam.setup to get the basis functions and the S matrices
  G <- mgcv:::gam.setup(gp,
                        pterms=pterms,
                        data=tmp,
                        sparse.cons=FALSE,
                        select=TRUE)
  
  rm(tmp)
  
  # Make list to store data
  out <- list()
  
  out$I       <- I
  out$K       <- K
  out$nT      <- nT
  out$N       <- length(c(responses))
  out$n_basis <- n_basis
  out$ii      <- rep(1:I, each = nT)
  out$tt      <- rep(1:nT, times = I)
  out$X       <- G$X
  out$S       <- cbind(G$smooth[[1]]$S[[1]], G$smooth[[1]]$S[[2]])
  out$zero    <- rep(0, n_basis)
  out$y       <- c(responses)
  
  return(out)
}
