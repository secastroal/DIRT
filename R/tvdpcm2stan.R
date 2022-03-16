# tvdpcm2stan_data ----

# Function to create data list to run the tv-dpcm model in stan.

# The function requires the following arguments: The matrix resp with the 
# responses, the number of items I, the number of response categories K,
# the number of timepoints nT, the number of knots n_knots, and the degree
# of the splines s_degree.
tvdpcm2stan_data <- function(resp, I, K, nT, n_knots, s_degree){
  # Get time indicator variable
  time <- 1:nT
  
  # Compute knots
  knots <- seq(min(time), max(time), length.out = n_knots)
  
  out <- list()
  
  out$I        <- I
  out$K        <- K
  out$nT       <- nT
  out$N        <- length(c(resp))
  out$N_obs    <- sum(!is.na(c(resp)))
  out$s_degree <- s_degree
  out$n_knots  <- n_knots
  out$knots    <- knots
  out$ii       <- rep(1:I, each = nT)
  out$tt       <- rep(1:nT, times = I)
  out$ii_obs   <- rep(1:I, each = nT)[!is.na(c(resp))]
  out$tt_obs   <- rep(1:nT, times = I)[!is.na(c(resp))]
  out$y_obs    <- c(resp)[!is.na(c(resp))]
  out$time     <- time
  
  return(out)
}