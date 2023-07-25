# This files contains functions to compute and plot the posterior predictive
# model checking methods for the time-varying dynamic partial credit model. 
# All the functions have as input the stanfit object, the data list as
# it is required for stan functions, and some extra arguments to control
# other parameters.
require(parallel)
require(foreach)
require(doParallel)
require(abind)
require(scatterpie)

# Load saved fit for testing.
# fit    <- readRDS("Fits/ppmc.fit.Rdata")
# object <- fit$fit
# data   <- fit$standata

# Function to combine foreach output to array.
# see stackoverflow.com/questions/17570806/parallel-for-loop-with-an-array-as-output
acomb <- function(...) {abind(..., along = 3)}

# Sumscores time series ----
ppmc.sumscore.ts <- function(object, data, mc.cores = getOption("mc.cores", 2L)) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- mclapply(as.data.frame(t(repy)), function(x) {
    tapply(x, data$tt_obs, sum)
  }, mc.cores = cores)
  
  sumscoresrepy <- as.matrix(as.data.frame(sumscoresrepy))

  plot(sort(unique(data$tt_obs)), sumscores, type = "p", pch = 19, las = 1,
       ylim = c(I - 0.5, I * K + 0.5), xlab = "Time",
       ylab = "Sum Scores")
  polygon(c(sort(unique(data$tt_obs)), 
            rev(sort(unique(data$tt_obs)))),
          c(apply(sumscoresrepy, 1, quantile, probs = 0.05), 
            rev(apply(sumscoresrepy, 1, quantile, probs = 0.95))),
          border = NA,
          col = rgb(1, 0, 0, 0.15))
  polygon(c(sort(unique(data$tt_obs)), 
            rev(sort(unique(data$tt_obs)))),
          c(apply(sumscoresrepy, 1, quantile, probs = 0.25), 
            rev(apply(sumscoresrepy, 1, quantile, probs = 0.75))),
          border = NA,
          col = rgb(1, 0, 0, 0.25))
  lines(sort(unique(data$tt_obs)), 
        apply(sumscoresrepy, 1, quantile, probs = 0.5), 
        col = "red")
  
  # Count observations out of the 95 percent of the replications
  out95 <- apply(sumscoresrepy, 1, quantile, probs = 0.05) <= sumscores &
    sumscores <= apply(sumscoresrepy, 1, quantile, probs = 0.95)
  prop_out <- sum(!out95) / length(out95)
  
  mtext(paste0("Prop Out = ", round(prop_out, 3)), side = 3, line = -1.5)
  
  return(prop_out)
}

# Partial Autocorrelation of the residuals ----
ppmc.racf <- function(object, data, col.ppp = "black", xlab = NULL,
                      col.y = "black", col.yrep = "lightgray", 
                      mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- mclapply(as.data.frame(t(repy)), function(x) {
    tapply(x, data$tt_obs, sum)
  }, mc.cores = cores)
  
  res  <- rstandard(lm(sumscores[-1] ~ head(sumscores, -1))) 
  acor <- acf(res, lag.max = 1, plot = FALSE)$acf[2, ,]
  
  acorrep <- mclapply(sumscoresrepy, function(x) {
    res <- rstandard(lm(x[-1] ~ head(x, -1)))
    out <- acf(res, lag.max = 1, plot = FALSE)$acf[2, ,]
    return(out)
  }, mc.cores = cores)
  
  acorrep <- unlist(acorrep)
  
  # Compute posterior predictive p-values
  out <- sum(acorrep >= acor) / length(acorrep)
  
  if (is.null(xlab)) {xlab = expression(PACF(y^rep))}
  
  hist(acorrep, main = "",
       xlab = xlab,
       xlim = c(min(acorrep, acor) - IQR(acorrep)/3,
                max(acorrep, acor) + IQR(acorrep)/3),
       las = 1, freq = FALSE, col = col.yrep, ...)
  abline(v = acor, lwd = 3, col = col.y)
  mtext(paste0("  PPP = ", round(out, 3)), line = -1.5, col = col.ppp, 
        cex = 0.8, adj = 0)
  
  return(round(out, 3))
}

# Autocorrelation of the sumscores ----
ppmc.acf <- function(object, data, lag.max = 5, col.ppp = "black", 
                     col.y = "black", col.yrep = "lightgray", 
                     xlab = NULL, default.xlab = is.null(xlab),  
                     mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- mclapply(as.data.frame(t(repy)), function(x) {
    tapply(x, data$tt_obs, sum)
  }, mc.cores = cores)
  
  acor <- acf(sumscores, lag.max = lag.max, plot = FALSE)$acf[2:(lag.max + 1), ,]
  acorrep <- mclapply(sumscoresrepy, function(x) {
    acf(x, lag.max = lag.max, plot = FALSE)$acf[2:(lag.max + 1), ,]
  }, mc.cores = cores)
  
  acorrep <- as.matrix(as.data.frame(acorrep))
  
  # Compute posterior predictive p-values
  out <- apply(acorrep >= acor, 1, function(x) sum(x)/dim(acorrep)[2])
  
  for (i in 1:lag.max) {
    if (default.xlab) {xlab = bquote(ACF[.(i)](y^rep))}
    
    hist(acorrep[i, ], main = "",
         xlab = xlab,
         xlim = c(min(acorrep[i, ], acor) - IQR(acorrep[i, ])/3,
                  max(acorrep[i, ], acor) + IQR(acorrep[i, ])/3),
         las = 1, freq = FALSE, col = col.yrep, ...)
    abline(v = acor[i], lwd = 3, col = col.y)
    mtext(paste0("  PPP = ", round(out[i], 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
  }
  
  return(round(out, 3))
}

# Mean Square Successive Difference ----

ppmc.mssd <- function(object, data, col.ppp = "black", xlab = NULL, 
                      col.y = "black", col.yrep = "lightgray", 
                      mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- mclapply(as.data.frame(t(repy)), function(x) {
    tapply(x, data$tt_obs, sum)
  }, mc.cores = cores)
  
  # Compute mssd
  mssd <- sum(diff(sumscores) ^ 2)/(nT - 1)
  
  mssdrep <- mclapply(sumscoresrepy, function(x) {
    out <- sum(diff(x) ^ 2)/(nT - 1)
    return(out)
  }, mc.cores = cores)
  
  mssdrep <- unlist(mssdrep)
  
  # Compute posterior predictive p-values
  out <- sum(mssdrep >= mssd) / length(mssdrep)
  
  if (is.null(xlab)) {xlab = expression(paste("MSSD(", y^rep, ")"))}
  
  hist(mssdrep, main = "",
       xlab = xlab,
       xlim = c(min(mssdrep, mssd) - IQR(mssdrep)/3,
                max(mssdrep, mssd) + IQR(mssdrep)/3),
       las = 1, freq = FALSE, col = col.yrep, ...)
  abline(v = mssd, lwd = 3, col = col.y)
  mtext(paste0("  PPP = ", round(out, 3)), line = -1.5, col = col.ppp, 
        cex = 0.8, adj = 0)
  
  return(round(out, 3))
}

# Item score time series ----
ppmc.item.ts <- function(object, data, quiet = FALSE, items = NULL) {
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  prop_out <- rep(NA, length(items))
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    plot(data$tt_obs[data$ii_obs == items[i]], y[data$ii_obs == items[i]], 
         type = "p", pch = 19, las = 1,
         ylim = c(min(y) - 0.5, max(y) + 0.5), xlab = "Time",
         ylab = paste0("Scores of Item ", i))
    polygon(c(data$tt_obs[data$ii_obs == items[i]], 
              rev(data$tt_obs[data$ii_obs == items[i]])),
            c(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.05), 
              rev(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.95))),
            border = NA,
            col = rgb(1, 0, 0, 0.15))
    polygon(c(data$tt_obs[data$ii_obs == items[i]], 
              rev(data$tt_obs[data$ii_obs == items[i]])),
            c(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.25), 
              rev(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.75))),
            border = NA,
            col = rgb(1, 0, 0, 0.25))
    lines(data$tt_obs[data$ii_obs == items[i]], 
          apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.5), 
          col = "red")
    
    # Count observations out of the 95 percent of the replications
    out95 <- apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.05) <= y[data$ii_obs == items[i]] &
      y[data$ii_obs == items[i]] <= apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.95)
    prop_out[items[i]] <- sum(!out95) / length(out95)
    
    mtext(paste0("Prop Out = ", round(prop_out[items[i]], 3)), side = 3, line = -1.5)
  }
  
  return(round(prop_out, 3))
}

# Item total correlation (Not modified) ----
# This is the typical item-total correlation in which the item scores of item j
# are correlated with the sumscores of the scale minus de items scores of 
# such item j. 
# The correlation between item scores and rescores can be either polyserial
# or pearson correlation. 
# There might be warning messages with the polyserial correlation because the 
# estimated correlation is larger than 1. For example, "initial correlation 
# inadmissible, 1.01866242085942, set to 0.9999"
ppmc.itcor <- function(object, data, method = c("polyserial", "pearson"), 
                       items = NULL, col.ppp = "black", col.y = "black", 
                       col.yrep = "lightgray", quiet = FALSE, 
                       xlab = NULL, default.xlab = is.null(xlab),  
                       mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  if (method == "polyserial") {
    polcor <- rep(NA, length(items))
    for(i in 1:length(items)) {
      polcor[items[i]] <- polycor::polyserial(tapply(y[data$ii_obs != items[i]], 
                                                     data$tt_obs[data$ii_obs != items[i]], 
                                                     sum), 
                                              y[data$ii_obs == items[i]])
    }
    
    repdat <- mclapply(as.data.frame(t(repy)), function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        out[items[i]] <- polycor::polyserial(tapply(x[data$ii_obs != items[i]], 
                                                    data$tt_obs[data$ii_obs != items[i]], 
                                                    sum), 
                                             x[data$ii_obs == items[i]])
      }
      return(out)
    }, mc.cores = cores)
  }
  
  if (method == "pearson") {
    polcor <- rep(NA, length(items))
    for(i in 1:length(items)) {
      polcor[items[i]] <- cor(tapply(y[data$ii_obs != items[i]],
                                     data$tt_obs[data$ii_obs != items[i]], 
                                     sum),
                              y[data$ii_obs == items[i]])
    }
    
    repdat <- mclapply(as.data.frame(t(repy)), function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        out[items[i]] <- cor(tapply(x[data$ii_obs != items[i]],
                                    data$tt_obs[data$ii_obs != items[i]], 
                                    sum),
                             x[data$ii_obs == items[i]])
        if (is.na(out[items[i]])) {out[items[i]] <- polcor[items[i]]}
      }
      return(out)
    }, mc.cores = cores)
  }
  
  repdat <- as.matrix(as.data.frame(repdat))
  
  # Compute posterior predictive p-values
  out <- apply(repdat >= polcor, 1, function(x) sum(x, na.rm = TRUE)/dim(repdat)[2])
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    if (default.xlab) {xlab = substitute(paste("Item-Total Correlation Item ", 
                                               ii, " (", y^rep, ")"), 
                                         list(ii = items[i]))}
    
    hist(repdat[items[i], ], 
         main = "",
         xlab = xlab,
         xlim = c(min(repdat[items[i], ], polcor[items[i]], na.rm = TRUE)
                  - IQR(repdat[items[i], ]/3),
                  max(repdat[items[i], ], polcor[items[i]], na.rm = TRUE)
                  + IQR(repdat[items[i], ]/3)),
         las = 1, freq = FALSE, col = col.yrep, ...)
    abline(v = polcor[items[i]], lwd = 3, col = col.y)
    mtext(paste0("  PPP = ", round(out[items[i]], 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
  }
  
  # return posterior predictive p-values
  return(round(out, 3))
}

# Item-total Correlation Version 2 ----
# In order to account for the time dependency of the observations, we modified the 
# item-total correlation so instead of correlating the item scores with the
# rescores, we correlate the item scores with the residuals of an autoregressive 
# model computed on rescores. 
# Again, one can use either polyserial or pearson correlation.
# Apparently with this methods, we do not get the warnings about correlations 
# larger than 1, when using polyserial correlation.
ppmc.itcor2 <- function(object, data, method = c("polyserial", "pearson"), 
                        items = NULL, quiet = FALSE, col.ppp = "black", 
                        col.y = "black", col.yrep = "lightgray",
                        xlab = NULL, default.xlab = is.null(xlab),  
                        mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  if (method == "polyserial") {
    polcor <- rep(NA, length(items))
    for(i in 1:length(items)) {
      sums <- tapply(y[data$ii_obs != items[i]], 
                     data$tt_obs[data$ii_obs != items[i]], 
                     sum)
      res  <- rstandard(lm(sums[-1] ~ head(sums, -1)))
      polcor[items[i]] <- polycor::polyserial(res, 
                                              y[data$ii_obs == items[i]][-1])
    }
    
    repdat <- mclapply(as.data.frame(t(repy)), function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        sums <- tapply(x[data$ii_obs != items[i]], 
                       data$tt_obs[data$ii_obs != items[i]], 
                       sum)
        res  <- rstandard(lm(sums[-1] ~ head(sums, -1)))
        out[items[i]] <- polycor::polyserial(res, 
                                             x[data$ii_obs == items[i]][-1])
      }
      return(out)
    }, mc.cores = cores)
  }

  if (method == "pearson") {
    polcor <- rep(NA, length(items))
    for(i in 1:length(items)) {
      sums <- tapply(y[data$ii_obs != items[i]], 
                     data$tt_obs[data$ii_obs != items[i]], 
                     sum)
      res  <- rstandard(lm(sums[-1] ~ head(sums, -1)))
      polcor[items[i]] <- cor(res,
                              y[data$ii_obs == items[i]][-1])
    }
    
    repdat <- mclapply(as.data.frame(t(repy)), function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        sums <- tapply(x[data$ii_obs != items[i]], 
                       data$tt_obs[data$ii_obs != items[i]], 
                       sum)
        res  <- rstandard(lm(sums[-1] ~ head(sums, -1)))
        out[items[i]] <- cor(res,
                             x[data$ii_obs == items[i]][-1])
      }
      return(out)
    }, mc.cores = cores)
  }
  
  repdat <- as.matrix(as.data.frame(repdat))
  
  # Compute posterior predictive p-values
  out <- apply(repdat >= polcor, 1, function(x) sum(x, na.rm = TRUE)/dim(repdat)[2])
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    if (default.xlab) {xlab = substitute(paste("Item-Total Correlation (v2) Item ", 
                                               ii, " (", y^rep, ")"), 
                                         list(ii = items[i]))}
    
    hist(repdat[items[i], ], 
         main = "",
         xlim = c(min(repdat[items[i], ], polcor[items[i]], na.rm = TRUE) -
                    IQR(repdat[items[i], ])/3,
                  max(repdat[items[i], ], polcor[items[i]], na.rm = TRUE) +
                    IQR(repdat[items[i], ])/3),
         xlab = xlab,
         las = 1, freq = FALSE, col = col.yrep, ...)
    abline(v = polcor[items[i]], lwd = 3, col = col.y)
    mtext(paste0("  PPP = ", round(out[items[i]], 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
  }
  
  # return posterior predictive p-values
  return(round(out, 3))
}


# Item-total Correlation Version 3 ----
# In this version, we fit a autoregressive model to both the item scores and 
# the rescores. Then, the residuals from these model are correlated.
ppmc.itcor3 <- function(object, data, method = "pearson", 
                        items = NULL, quiet = FALSE, col.ppp = "black", 
                        col.y = "black", col.yrep = "lightgray", 
                        xlab = NULL, default.xlab = is.null(xlab),  
                        mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  if (method == "pearson") {
    polcor <- rep(NA, length(items))
    for(i in 1:length(items)) {
      sums <- tapply(y[data$ii_obs != items[i]], 
                     data$tt_obs[data$ii_obs != items[i]], 
                     sum)
      res  <- rstandard(lm(sums[-1] ~ head(sums, -1)))
      ires <- rstandard(lm(y[data$ii_obs == items[i]][-1] ~ 
                             head(y[data$ii_obs == items[i]], -1)))
      polcor[items[i]] <- cor(res, ires)
    }
    
    repdat <- mclapply(as.data.frame(t(repy)), function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        sums <- tapply(x[data$ii_obs != items[i]], 
                       data$tt_obs[data$ii_obs != items[i]], 
                       sum)
        res  <- rstandard(lm(sums[-1] ~ head(sums, -1)))
        ires <- rstandard(lm(x[data$ii_obs == items[i]][-1] ~ 
                               head(x[data$ii_obs == items[i]], -1)))
        out[items[i]] <- cor(res, ires)
      }
      return(out)
    }, mc.cores = cores)
  }
  
  repdat <- as.matrix(as.data.frame(repdat))
  
  # Compute posterior predictive p-values
  out <- apply(repdat >= polcor, 1, function(x) sum(x, na.rm = TRUE)/dim(repdat)[2])
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    if (default.xlab) {xlab = substitute(paste("Item-Total Correlation (v3) Item ", 
                                               ii, " (", y^rep, ")"), 
                                         list(ii = items[i]))}
    
    hist(repdat[items[i], ], 
         main = "",
         xlim = c(min(repdat[items[i], ], polcor[items[i]], na.rm = TRUE) -
                    IQR(repdat[items[i], ])/3,
                  max(repdat[items[i], ], polcor[items[i]], na.rm = TRUE) +
                    IQR(repdat[items[i], ])/3),
         xlab = xlab,
         las = 1, freq = FALSE, col = col.yrep, ...)
    abline(v = polcor[items[i]], lwd = 3, col = col.y)
    mtext(paste0("  PPP = ", round(out[items[i]], 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
  }
  
  # return posterior predictive p-values
  return(round(out, 3))
}

# Yen's Q1 unmodified ----

# Function to compute the Yen's Q1 for polytomous models
# The arguments are the estimated theta, threshold, and discrimination 
# parameters, the number of items, the number of response categories,
# the grouping variable for each time point, and a grouping and item index for
# the data in long format.
gpcm.Q1 <- function(y, theta, thresholds, alpha, I, K, group, group_index, i_index) {
  
  M <- K - 1
  
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (yy in 0:M) {
    probs.array[, , yy + 1] <- P.GPCM(y     = yy, 
                                      alpha = alpha, 
                                      delta = delta, 
                                      taus  = taus, 
                                      theta = theta, 
                                      M     = M)
  }
  
  E <- apply(probs.array, c(3, 2), function(x) {
    tapply(x, group, mean)
    })
  
  O <-  tapply(y, list(group_index, i_index), function(x) {
    x    <- factor(x, levels = 1:K)
    xtab <- prop.table(table(x, exclude = FALSE))
    return(as.vector(xtab))
    })
  
  q1 <- rep(0, I)
  g_size <- table(group)
  
  for (i in 1:I) {
    for (g in 1:10) {
      q1[i] <- q1[i] + sum(g_size[g] * ((O[g, i][[1]] - E[g, , i])^2 / E[g, , i]))
    }
  }
  
  return(q1)
}

ppmc.Q1 <- function(object, data, items = NULL, quiet = FALSE, 
                    col.ppp = "black", col.abline = "lightgray", lty.abline = 2,
                    xlab = NULL, default.xlab = is.null(xlab),  
                    ylab = NULL, default.ylab = is.null(ylab),  
                    mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  # Get the observed timepoints. Relevant if there were missing data.
  times_obs <- intersect(data$tt, data$tt_obs)
  
  # Parallel foreach loop
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl, cores = cores)
  
  discrepancy <- foreach(r = 1:nrow(repy), .combine = "acomb",
                         .export = c("gpcm.Q1", "P.GPCM")) %dopar% {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ][times_obs]
    
    # Create groups for Yen's Q1
    tmpx <- order(theta)
    tmpx <- cbind(tmpx, cut(seq_along(tmpx), 10, labels = FALSE))
    tmpx <- tmpx[order(tmpx[, 1]), ]
    
    group       <- tmpx[, 2]
    group_index <- rep(tmpx[, 2], I)
    rm(tmpx)
    
    result <- matrix(NA, I, 2)
    
    # Yen's Q1 for the observed scores
    result[, 1] <- gpcm.Q1(y           = y, 
                           theta       = theta, 
                           thresholds  = thresholds, 
                           alpha       = rep(1, I), 
                           I           = I, 
                           K           = K, 
                           group       = group,
                           group_index = group_index,
                           i_index     = data$ii_obs)
    # Yen's Q1 for the i-th replicated scores
    result[, 2] <- gpcm.Q1(y           = repy[r, ], 
                           theta       = theta, 
                           thresholds  = thresholds, 
                           alpha       = rep(1, I), 
                           I           = I, 
                           K           = K, 
                           group       = group,
                           group_index = group_index,
                           i_index     = data$ii_obs)
    
    result
  }
  
  stopCluster(cl = cl)
  
  discrepancy <- aperm(discrepancy, c(3, 1, 2))
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] >= discrepancy[, , 1], 2, mean)
  names(out) <- paste0("Item_", 1:I)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    if (default.ylab) {ylab = expression(paste("Yen's ", Q[1], "(", y^rep, ";", omega, ")"))}
    if (default.xlab) {xlab = expression(paste("Yen's ", Q[1], "(y;", omega, ")"))}
    
    plot(discrepancy[, items[i], 1], discrepancy[, items[i], 2], 
         las = 1, main = "",
         ylab = ylab,
         xlab = xlab,
         ...)
    abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
    mtext(paste0("  PPP = ", round(out[items[i]], 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
    mtext(paste("Item", items[i]), line = 0.5, cex = 1.2, font = 2, adj = 1)
  }
  
  return(round(out, 3))
  }

# Yen's Q1 modified ----
# This is an alternative version of Yen's Q1 in which the groups are not made by
# ordering the estimated latent values. Instead, the groups are simply defined 
# by the time in which the observation were made. The purpose of this is to take 
# the time into account.

ppmc.Q1.alt <- function(object, data, items = NULL, quiet = FALSE,
                        col.ppp = "black", col.abline = "lightgray", lty.abline = 2,
                        xlab = NULL, default.xlab = is.null(xlab),  
                        ylab = NULL, default.ylab = is.null(ylab),  
                        mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  # Get the observed timepoints. Relevant if there were missing data.
  times_obs <- intersect(data$tt, data$tt_obs)
  
  # Parallel foreach loop
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl, cores = cores)
  
  discrepancy <- foreach(r = 1:nrow(repy), .combine = "acomb",
                         .export = c("gpcm.Q1", "P.GPCM")) %dopar% {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ][times_obs]
    
    # Create groups for Yen's Q1 given the time order
    group       <- cut(seq_along(theta), 10, labels = FALSE)
    group_index <- rep(group, I)
    
    result <- matrix(NA, I, 2)
    
    # Yen's Q1 for the observed scores
    result[, 1] <- gpcm.Q1(y           = y, 
                           theta       = theta, 
                           thresholds  = thresholds, 
                           alpha       = rep(1, I), 
                           I           = I, 
                           K           = K, 
                           group       = group,
                           group_index = group_index,
                           i_index     = data$ii_obs)
    # Yen's Q1 for the i-th replicated scores
    result[, 2] <- gpcm.Q1(y           = repy[r, ], 
                           theta       = theta, 
                           thresholds  = thresholds, 
                           alpha       = rep(1, I), 
                           I           = I, 
                           K           = K, 
                           group       = group,
                           group_index = group_index,
                           i_index     = data$ii_obs)
    result
  }
  
  stopCluster(cl = cl)
  
  discrepancy <- aperm(discrepancy, c(3, 1, 2))
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] >= discrepancy[, , 1], 2, mean)
  names(out) <- paste0("Item_", 1:I)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    if (default.ylab) {ylab = expression(paste("Alt. Yen's ", Q[1], "(", y^rep, ";", omega, ")"))}
    if (default.xlab) {xlab = expression(paste("Alt. Yen's ", Q[1], "(y;", omega, ")"))}
    
    plot(discrepancy[, items[i], 1], discrepancy[, items[i], 2], 
         las = 1, main = "",
         ylab = ylab,
         xlab = xlab, ...)
    abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
    mtext(paste0("  PPP = ", round(out[items[i]], 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
    mtext(paste("Item", items[i]), line = 0.5, cex = 1.2, font = 2, adj = 1)
  }
  
  return(round(out, 3))
}

# Yen's Q3 ----

# First we create a function to compute Yen's Q3. The function requires the 
# vector of responses y, the estimated theta, thresholds, and discrimination
# parameters, the number of time points, the number of items, the number of 
# response categories, and two index variables t_index for time and i_index 
# for items given that the data is in long format.

gpcm.Q3 <- function(y, theta, thresholds, alpha, I, K, t_index) {
  
  M <- K - 1
  
  # Restructure responses in a matrix
  Y <- matrix(y, length(unique(t_index)), I)
  
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  theta <- theta[unique(t_index)]
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (yy in 0:M) {
    probs.array[, , yy + 1] <- P.GPCM(y    = yy, 
                                      alpha = alpha, 
                                      delta = delta, 
                                      taus  = taus, 
                                      theta = theta, 
                                      M     = M)
  }
  
  E <- apply(probs.array, c(1, 2), function(x) sum(x * 1:K))
  
  D <- Y - E
  
  CorD <- cor(D, use = "pairwise.complete.obs")
  q3 <- CorD[lower.tri(CorD)]
  names(q3) <- apply(which(lower.tri(CorD), arr.ind = TRUE), 1, function(x)
    paste0("(", paste(x, collapse = ","), ")"))
  
  return(q3)
}

ppmc.Q3 <- function(object, data, scatterplots = FALSE,
                    col.ppp = "black", col.abline = "lightgray", lty.abline = 2,
                    xlab = NULL, default.xlab = is.null(xlab),  
                    ylab = NULL, default.ylab = is.null(ylab),  
                    mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  # Parallel foreach loop
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl, cores = cores)
  
  discrepancy <- foreach(r = 1:nrow(repy), .combine = "acomb",
                         .export = c("gpcm.Q3", "P.GPCM")) %dopar% {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ]
    
    result <- matrix(NA, (I * (I - 1))/2, 2)
    
    # Yen's Q3 for the observed scores
    result[, 1] <- gpcm.Q3(y          = y,
                           theta      = theta,
                           thresholds = thresholds,
                           alpha      = rep(1, I),
                           I          = I,
                           K          = K,
                           t_index    = data$tt_obs)
    
    # Yen's Q3 for the i-th replicated scores
    result[, 2] <- gpcm.Q3(y          = repy[r, ],
                           theta      = theta,
                           thresholds = thresholds,
                           alpha      = rep(1, I),
                           I          = I,
                           K          = K,
                           t_index    = data$tt_obs)
    result
  }
  
  stopCluster(cl = cl)
  
  discrepancy <- aperm(discrepancy, c(3, 1, 2))
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] >= discrepancy[, , 1], 2, mean)
  names(out) <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
    paste0("ppp-q3(", paste(x, collapse = ","), ")"))
  
  if (scatterplots) {
    for (i in 1:length(out)) {
      if (default.ylab) {ylab = expression(paste("Yen's ", Q[3], "(", y^rep, ";", omega, ")"))}
      if (default.xlab) {xlab = expression(paste("Yen's ", Q[3], "(y;", omega, ")"))}
      
      plot(discrepancy[, i, 1], discrepancy[, i, 2], 
           las = 1, main = "",
           ylab = ylab,
           xlab = xlab, ...)
      abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
      mtext(paste0("  PPP = ", round(out[i], 3)), line = -1.5, col = col.ppp, 
            cex = 0.8, adj = 0)
      mtext(paste("Items", substring(names(out)[i], 7)), line = 0.5, cex = 1.2, 
            font = 2, adj = 1)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out <= 0.05 | tmp$out >= 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                             cols=c("out","out2"), color = NA,
                             show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    scale_x_continuous(breaks = 2:I) +
    scale_y_continuous(breaks = 1:(I - 1)) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(list(ppp = round(out, 3), sp = sp))
}

# Odds Ratio ----
# To compute the odds ratio for polytomous data, the items need to be 
# dichotomized. After this, the odds ratio can be compute as normally.
# Therefore, we first create to auxiliary functions, one used to dichotomize
# categorical responses and a second one to compute the odds ratio. These 
# functions assume that the responses are in long format.

# To dichotomize, this function requires the responses in long format, 
# an index variable indicating the item, the number of items, the number of
# response categories, and the cut value used to dichotomize the responses.
# The cutoff value might be different per item. The default cutoff is 
# ceiling(K/2) + 1
dicho <- function(y, I, K, i_index, cutoff = ceiling(K/2) + 1) {
  if (length(cutoff) == 1) {cutoff <- rep(cutoff, I)}
  
  out <- ifelse(y < cutoff[i_index], 0, 1)
  
  return(out)
}

# This function computes the odd ratio of a matrix. If any if the counts of the 
# denominator is equal to 0, the Haldane-Anscombe correction correction is 
# used (i.e., add 0.5 to each cell of the 2X2 matrix).
odds.ratio <- function(y, I, K, t_index) {
  
  # Restructure responses in a matrix
  Y <- matrix(y, length(unique(t_index)), I)
  
  # Vector to store Odds ratio
  or <- rep(NA, (I * (I - 1))/2)
  
  count <- 1
  for (i in 1:(I - 1)) {
    for (j in (i + 1):I) {
      n    <- rep(NA, 4)
      n[1] <- sum(Y[, i] == 1 & Y[, j] == 1, na.rm = TRUE)
      n[2] <- sum(Y[, i] == 0 & Y[, j] == 0, na.rm = TRUE)
      n[3] <- sum(Y[, i] == 1 & Y[, j] == 0, na.rm = TRUE)
      n[4] <- sum(Y[, i] == 0 & Y[, j] == 1, na.rm = TRUE)
      if (any(n[3:4] == 0)) {n <- n + 0.5}
      or[count] <- (n[1] * n[2])/(n[3] * n[4])
      names(or)[count] <- paste0('OR(', j, ',', i, ")")
      count <- count + 1
    }
  }
  
  return(or)
}

ppmc.OR <- function(object, data, cutoff = NULL, histograms = FALSE,
                    col.ppp = "black", col.y = "black", col.yrep = "lightgray",
                    xlab = NULL, default.xlab = is.null(xlab),  
                    mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(cutoff)) {cutoff <- rep(ceiling(K/2) + 1, I)}
  if (length(cutoff) == 1) {cutoff <- rep(cutoff, I)}
  
  ydich <- dicho(y = y, 
                 I = I, 
                 K = K, 
                 i_index = data$ii_obs, 
                 cutoff  = cutoff)
  or    <- odds.ratio(y  = ydich, 
                      I  = I, 
                      K  = K, 
                      t_index = data$tt_obs)
  
  
  orrep <- mclapply(as.data.frame(t(repy)), function(x) {
    repyd <- dicho(y = x, 
                   I = I, 
                   K = K, 
                   i_index = data$ii_obs, 
                   cutoff  = cutoff)
    out   <- odds.ratio(y  = repyd, 
                        I  = I, 
                        K  = K, 
                        t_index = data$tt_obs)
    return(out)
  }, mc.cores = cores)
  
  orrep <- as.matrix(as.data.frame(orrep))
  
  # Compute posterior predictive p-values
  out <- apply(orrep >= or, 1, function(x) sum(x)/dim(orrep)[2])
  
  if (histograms) {
    for (i in 1:length(or)) {
      if (default.xlab) {xlab = substitute(paste("OR", ii,"(", y^rep, ")"),
                                           list(ii = substring(names(out)[i], 3)))}
      
      hist(orrep[i, ], main = "",
           xlab = xlab,
           xlim = c(0, max(orrep[i, ], or[i]) + IQR(orrep[i, ])/3),
           las = 1, freq = FALSE, col = col.yrep, ...)
      abline(v = or[i], lwd = 3, col = col.y)
      mtext(paste0("  PPP = ", round(out[i], 3)), line = -1.5, col = col.ppp, 
            cex = 0.8, adj = 0)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out <= 0.05 | tmp$out >= 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                                   cols=c("out","out2"), color = NA,
                                   show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    scale_x_continuous(breaks = 2:I) +
    scale_y_continuous(breaks = 1:(I - 1)) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(list(ppp = round(out, 3), sp = sp))
}

# Odds Ratio Difference ----

ppmc.ORDiff <- function(object, data, cutoff = NULL, histograms = FALSE,
                        col.ppp = "black", col.y = "black", col.yrep = "lightgray",
                        xlab = NULL, default.xlab = is.null(xlab),  
                        mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  if (cores < 1L)
    stop("'mc.cores' must be >= 1")
  
  if (Sys.info()["sysname"] == "Windows") {
    if (cores > 1L)
      stop("'mc.cores' > 1 is not supported on Windows")
  }
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(cutoff)) {cutoff <- rep(ceiling(K/2) + 1, I)}
  if (length(cutoff) == 1) {cutoff <- rep(cutoff, I)}
  
  timesplit <- ifelse(data$tt_obs <= ceiling(nT/2), 0, 1)
  
  ydich  <- dicho(y = y, 
                  I = I, 
                  K = K, 
                  i_index = data$ii_obs, 
                  cutoff  = cutoff)
  
  orsplit <- matrix(NA, nrow = (I * (I - 1))/2, ncol = 2)
  
  for (s in 0:1) {
    tmp.nT <- ifelse(nT/2 - ceiling(nT/2) == 0, nT/2, ceiling(nT/2) - s)
    
    orsplit[, s + 1] <- odds.ratio(y  = ydich[timesplit == s],
                                   I  = I,
                                   K  = K,
                                   t_index = data$tt_obs[timesplit == s] - 
                                     ceiling(nT/2) * s)
  }
  
  ordiff <- c(diff(t(orsplit)))
  
  ordiffrep <- mclapply(as.data.frame(t(repy)), function(x) {
    repyd  <- dicho(y = x, 
                    I = I, 
                    K = K, 
                    i_index = data$ii_obs, 
                    cutoff  = cutoff)
    
    orsplitrep <- matrix(NA, nrow = (I * (I - 1))/2, ncol = 2)
    
    for (s in 0:1) {
      tmp.nT <- ifelse(nT/2 - ceiling(nT/2) == 0, nT/2, ceiling(nT/2) - s)
      
      orsplitrep[, s + 1] <- odds.ratio(y  = repyd[timesplit == s],
                                        I  = I,
                                        K  = K,
                                        t_index = data$tt_obs[timesplit == s] - 
                                          ceiling(nT/2) * s)
    }
    
    out <- as.vector(diff(t(orsplitrep)))
    return(out)
  }, mc.cores = cores)
  
  ordiffrep <- as.matrix(as.data.frame(ordiffrep))
  
  # Compute posterior predictive p-values
  out <- apply(ordiffrep >= ordiff, 1, function(x) sum(x)/dim(ordiffrep)[2])
  names(out) <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
    paste0("ppp-ordiff(", paste(x, collapse = ","), ")"))
  
  if (histograms) {
    for (i in 1:length(ordiff)) {
      if (default.xlab) {xlab = substitute(paste("ORDiff", ii,"(", y^rep, ")"), 
                                           list(ii = substring(names(out)[i], 11)))}
      
      hist(ordiffrep[i, ], main = "",
           xlab = xlab,
           xlim = c(min(ordiffrep[i, ], ordiff[i]) - IQR(ordiffrep[i, ])/3, 
                    max(ordiffrep[i, ], ordiff[i]) + IQR(ordiffrep[i, ])/3),
           las = 1, freq = FALSE, col = col.yrep, ...)
      abline(v = ordiff[i], lwd = 3, col = col.y)
      mtext(paste0("  PPP = ", round(out[i], 3)), line = -1.5, col = col.ppp, 
            cex = 0.8, adj = 0)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out <= 0.05 | tmp$out >= 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                                   cols=c("out","out2"), color = NA,
                                   show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    scale_x_continuous(breaks = 2:I) +
    scale_y_continuous(breaks = 1:(I - 1)) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(list(ppp = round(out, 3), sp = sp))
}

# Absolute Item Covariance Residual - RESID ----
# This function computes the absolute item covariance residual. It requires
# the vector of observed responses in long format, the estimated theta, 
# threshold, and discrimination parameters, the number of time points, 
# the number of items, the number of responses categories, and index variables 
# for the time points and the items.
cov.resid <- function(y, theta, thresholds, alpha, I, K, t_index) {
  M <- K - 1
  
  # Restructure responses in a matrix
  Y <- matrix(y, length(unique(t_index)), I)
  
  theta <- theta[unique(t_index)]
  
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (yy in 0:M) {
    probs.array[, , yy + 1] <- P.GPCM(y    = yy, 
                                      alpha = alpha, 
                                      delta = delta, 
                                      taus  = taus, 
                                      theta = theta, 
                                      M     = M)
  }
  
  E <- apply(probs.array, c(1, 2), function(x) sum(x * 1:K))
  
  out <- abs(cov(Y, use = "pairwise.complete.obs") - cov(E))
  out <- out[lower.tri(out)]
  names(out) <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
    paste0("(", paste(x, collapse = ","), ")"))
  
  return(out)
}

ppmc.cov.resid <- function(object, data, scatterplots = FALSE, 
                           col.ppp = "black", col.abline = "lightgray", lty.abline = 2,
                           xlab = NULL, default.xlab = is.null(xlab),  
                           ylab = NULL, default.ylab = is.null(ylab),  
                           mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  # Parallel foreach loop
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl, cores = cores)
  
  discrepancy <- foreach(r = 1:nrow(repy), .combine = "acomb",
                         .export = c("cov.resid", "P.GPCM")) %dopar% {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ]
    
    result <- matrix(NA, (I * (I - 1))/2, 2)
    
    # RESID for the observed scores
    result[, 1] <- cov.resid(y          = y,
                             theta      = theta,
                             thresholds = thresholds,
                             alpha      = rep(1, I),
                             I          = I,
                             K          = K,
                             t_index    = data$tt_obs)
    
    # RESID for the i-th replicated scores
    result[, 2] <- cov.resid(y          = repy[r, ],
                             theta      = theta,
                             thresholds = thresholds,
                             alpha      = rep(1, I),
                             I          = I,
                             K          = K,
                             t_index    = data$tt_obs)
    result
  }
  
  stopCluster(cl = cl)
  
  discrepancy <- aperm(discrepancy, c(3, 1, 2))
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] >= discrepancy[, , 1], 2, mean)
  names(out) <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
    paste0("ppp-resid(", paste(x, collapse = ","), ")"))
  
  if (scatterplots) {
    for (i in 1:length(out)) {
      if (default.ylab) {ylab = expression(paste("Cov RESID (", y^rep, ";", omega, ")"))}
      if (default.xlab) {xlab = expression(paste("Cov RESID (y;", omega, ")"))}
      
      plot(discrepancy[, i, 1], discrepancy[, i, 2], 
           las = 1, main = "",
           ylab = ylab,
           xlab = xlab, ...)
      abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
      mtext(paste0("  PPP = ", round(out[i], 3)), line = -1.5, col = col.ppp, 
            cex = 0.8, adj = 0)
      mtext(paste("Items", substring(names(out)[i], 10)), line = 0.5, cex = 1.2, 
            font = 2, adj = 1)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out <= 0.05 | tmp$out >= 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                                   cols=c("out","out2"), color = NA,
                                   show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    scale_x_continuous(breaks = 2:I) +
    scale_y_continuous(breaks = 1:(I - 1)) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(list(ppp = round(out, 3), sp = sp))
}

# Absolute Item Covariance Residual Difference ----

ppmc.cov.rediff <- function(object, data, scatterplots = FALSE,
                            col.ppp = "black", col.abline = "lightgray", lty.abline = 2,
                            xlab = NULL, default.xlab = is.null(xlab),  
                            ylab = NULL, default.ylab = is.null(ylab),  
                            mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  # Create dummy variable to split observations in two halves
  timesplit <- ifelse(data$tt_obs <= ceiling(nT/2), 0, 1)
  
  # Parallel foreach loop
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl, cores = cores)
  
  discrepancy <- foreach(r = 1:nrow(repy), .combine = "acomb",
                         .export = c("cov.resid", "P.GPCM")) %dopar% {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ]
    
    tmp <- array(NA, dim = c((I * (I - 1))/2, 2, 2))
    
    # Loop through the two halves.
    for (s in 0:1) {
      tmp.nT <- ifelse(nT/2 - ceiling(nT/2) == 0, nT/2, ceiling(nT/2) - s)
      
      # RESID difference of the observed scores
      tmp[, 1, s + 1] <- cov.resid(y          = y[timesplit == s],
                                   theta      = theta,
                                   thresholds = thresholds,
                                   alpha      = rep(1, I),
                                   I          = I,
                                   K          = K,
                                   t_index    = data$tt_obs[timesplit == s] - 
                                     ceiling(nT/2) * s)
      
      # RESID difference of the replicated scores
      tmp[, 2, s + 1] <- cov.resid(y          = repy[r, ][timesplit == s],
                                   theta      = theta,
                                   thresholds = thresholds,
                                   alpha      = rep(1, I),
                                   I          = I,
                                   K          = K,
                                   t_index    = data$tt_obs[timesplit == s] - 
                                     ceiling(nT/2) * s)
    }
    
    result <- matrix(NA, (I * (I - 1))/2, 2)
    
    result[, 1] <- c(diff(t(tmp[, 1, ])))
    result[, 2] <- c(diff(t(tmp[, 2, ])))
    
    result
  }
    
  stopCluster(cl = cl)
  
  discrepancy <- aperm(discrepancy, c(3, 1, 2))
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] >= discrepancy[, , 1], 2, mean)
  names(out) <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
    paste0("ppp-rediff(", paste(x, collapse = ","), ")"))
  
  if (scatterplots) {
    for (i in 1:length(out)) {
      if (default.ylab) {ylab = expression(paste("Cov REDIFF (", y^rep, ";", omega, ")"))}
      if (default.xlab) {xlab = expression(paste("Cov REDIFF (y;", omega, ")"))}
      
      plot(discrepancy[, i, 1], discrepancy[, i, 2], 
           las = 1, main = "",
           ylab = ylab,
           xlab = xlab, ...)
      abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
      mtext(paste0("  PPP = ", round(out[i], 3)), line = -1.5, col = col.ppp, 
            cex = 0.8, adj = 0)
      mtext(paste("Items", substring(names(out)[i], 11)), line = 0.5, cex = 1.2, 
            font = 2, adj = 1)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out <= 0.05 | tmp$out >= 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                                   cols=c("out","out2"), color = NA,
                                   show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    scale_x_continuous(breaks = 2:I) +
    scale_y_continuous(breaks = 1:(I - 1)) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(list(ppp = round(out, 3), sp = sp))
}

# Partial Autocorrelation of the Latent Scores Residuals ----

gpcm.lpacf <- function(y, theta, thresholds, alpha, I, K, 
                       t_index, sumscores = FALSE) {
  
  M <- K - 1
  
  # Restructure responses in a matrix
  Y <- matrix(y, length(unique(t_index)), I)
  
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (yy in 0:M) {
    probs.array[, , yy + 1] <- P.GPCM(y    = yy, 
                                      alpha = alpha, 
                                      delta = delta, 
                                      taus  = taus, 
                                      theta = theta, 
                                      M     = M)
  }
  
  E <- apply(probs.array, c(1, 2), function(x) sum(x * 1:K))
  
  if (sumscores) {
    D <- rowMeans(na.omit(Y)) - rowMeans(E)
    
    out <- acf(D, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2, ,]
  } else {
    D <- na.omit(Y) - E
    
    out <- apply(D, 2, function(x) {
      acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2, ,]
    })
  }
 
  return(out)
}

ppmc.lpacf <- function(object, data, items = NULL, quiet = FALSE, sumscores = FALSE, 
                       col.ppp = "black", col.abline = "lightgray", lty.abline = 2,
                       xlab = NULL, default.xlab = is.null(xlab),  
                       ylab = NULL, default.ylab = is.null(ylab), 
                       subtitle = TRUE, col = "darkgray", split = FALSE, 
                       mc.cores = getOption("mc.cores", 2L), ...) {
  
  cores <- as.integer(mc.cores)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  # Get the observed timepoints. Relevant if there were missing data.
  times_obs <- intersect(data$tt, data$tt_obs)
  
  # Parallel foreach loop
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl, cores = cores)
  
  discrepancy <- foreach(r = 1:nrow(repy), .combine = "acomb",
                         .export = c("gpcm.lpacf", "P.GPCM")) %dopar% {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ][times_obs]
    
    if (sumscores) {
      result <- matrix(NA, 1, 2)
    } else {
      result <- matrix(NA, I, 2)
    }
    
    # pacf for the observed scores
    result[, 1] <- gpcm.lpacf(y          = y,
                              theta      = theta,
                              thresholds = thresholds,
                              alpha      = rep(1, I),
                              I          = I,
                              K          = K,
                              t_index    = data$tt_obs,
                              sumscores  = sumscores)
    
    # pacf for the i-th replicated scores
    result[, 2] <- gpcm.lpacf(y          = repy[r, ],
                              theta      = theta,
                              thresholds = thresholds,
                              alpha      = rep(1, I),
                              I          = I,
                              K          = K,
                              t_index    = data$tt_obs,
                              sumscores  = sumscores)
    result
  }
  
  stopCluster(cl = cl)
  
  discrepancy <- aperm(discrepancy, c(3, 1, 2))
  
  # Compute posterior predictive p-values
  if (sumscores) {
    out  <- mean(discrepancy[, 1, 2] >= discrepancy[, 1, 1])
    names(out) <- "ppp"
    
    if (default.ylab) {ylab = expression(paste("LPACF", "(", y^rep, ";", omega, ")"))}
    if (default.xlab) {xlab = expression(paste("LPACF", "(y;", omega, ")"))}
    
    if (split) {
      col <- ifelse(discrepancy[, 1, 1] < discrepancy[, 1, 2], col[1], col[2])
    }
    
    plot(discrepancy[, 1, 1], discrepancy[, 1, 2], 
         las = 1, main = "",
         ylab = ylab,
         xlab = xlab, col = col,  ...)
    abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
    mtext(paste0("  PPP = ", round(out, 3)), line = -1.5, col = col.ppp, 
          cex = 0.8, adj = 0)
    if (subtitle) {
      mtext("Sumscores", line = 0.5, cex = 1.2, 
            font = 2, adj = 1)
    }
  } else {
    out  <- apply(discrepancy[, , 2] >= discrepancy[, , 1], 2, mean)
    names(out) <- paste0("Item_", 1:I)
    
    if (default.ylab) {ylab = expression(paste("LPACF", "(", y^rep, ";", omega, ")"))}
    if (default.xlab) {xlab = expression(paste("LPACF", "(y;", omega, ")"))}
    
    for (i in 1:length(items)) {
      if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
      if (split) {
        col <- ifelse(discrepancy[, items[i], 1] < discrepancy[, items[i], 2], col[1], col[2])
      }
      plot(discrepancy[, items[i], 1], discrepancy[, items[i], 2], 
           las = 1, main = "",
           ylab = ylab,
           xlab = xlab, col = col, ...)
      abline(a= 0, b = 1, lwd = 3, col = col.abline, lty = lty.abline)
      mtext(paste0("  PPP = ", round(out[items[i]], 3)), line = -1.5, col = col.ppp, 
            cex = 0.8, adj = 0)
      if (subtitle) {
        mtext(paste("Item ", items[i]), line = 0.5, cex = 1.2, 
              font = 2, adj = 1)
      }
    }
  }
  
  return(round(out, 3))
}

# rm(list = setdiff(ls(), c(lsf.str(), "object", "data", "fit")))

