# This files contains functions to compute and plot the posterior predictive
# model checking methods for the time-varying dynamic partial credit model. 
# All the functions have as input the stanfit object, the data list as
# it is required for stan functions, and some extra arguments to control
# other parameters.

# Sumscores time series ----
ppc.sumscore.ts <- function(object, data) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- apply(repy, 1, function(x) {
    tapply(x, data$tt_obs, sum)
  })
  
  plot(sort(unique(data$tt_obs)), sumscores, type = "p", pch = 19, las = 1,
       ylim = c(I - 0.5, I * K + 0.5), xlab = "Time",
       ylab = "Sum Scores")
  polygon(c(sort(unique(data$tt_obs)), 
            rev(sort(unique(data$tt_obs)))),
          c(apply(sumscoresrepy, 1, quantile, probs = 0.025), 
            rev(apply(sumscoresrepy, 1, quantile, probs = 0.975))),
          border = NA,
          col = rgb(1, 0, 0, 0.15))
  polygon(c(1:nT, rev(1:nT)),
          c(apply(sumscoresrepy, 1, quantile, probs = 0.25), 
            rev(apply(sumscoresrepy, 1, quantile, probs = 0.75))),
          border = NA,
          col = rgb(1, 0, 0, 0.25))
  lines(apply(sumscoresrepy, 1, quantile, probs = 0.5), 
        col = "red")
  
  # Count observations out of the 95 percent of the replications
  out95 <- apply(sumscoresrepy, 1, quantile, probs = 0.025) <= sumscores &
    sumscores <= apply(sumscoresrepy, 1, quantile, probs = 0.975)
  prop_out <- sum(!out95) / length(out95)
  
  mtext(paste0("Prop Out = ", prop_out), side = 3, line = -1.5)
  
  return(prop_out)
}

# Autocorrelation of the residuals ----
ppc.racf <- function(object, data) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- apply(repy, 1, function(x) {
    tapply(x, data$tt_obs, sum)
  })
  
  res  <- rstandard(lm(sumscores[-1] ~ head(sumscores, -1))) 
  acor <- acf(res, lag.max = 1, plot = FALSE)$acf[2, ,]
  
  acorrep <- apply(sumscoresrepy, 2, function(x) {
    res <- rstandard(lm(x[-1] ~ head(x, -1)))
    out <- acf(res, lag.max = 1, plot = FALSE)$acf[2, ,]
    return(out)
  })
  
  # Compute posterior predictive p-values
  out <- sum(acorrep <= acor) / length(acorrep)
  
  hist(acorrep, main = "Histogram Residuals Autocorrelation",
       xlab = "Autocorrelation")
  abline(v = acor, lwd = 3)
  mtext(paste0("PPP = ", round(out, 3)), line = -1.5, col = "red", 
        cex = 0.8, adj = 0)
  
  return(round(out, 3))
}

# Autocorrelation of the sumscores ----
ppc.acf <- function(object, data, lag.max = 5) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- apply(repy, 1, function(x) {
    tapply(x, data$tt_obs, sum)
  })
  
  acor <- acf(sumscores, lag.max = lag.max, plot = FALSE)$acf[2:(lag.max + 1), ,]
  acorrep <- apply(sumscoresrepy, 2, function(x) {
    acf(x, lag.max = lag.max, plot = FALSE)$acf[2:(lag.max + 1), ,]
  })
  
  # Compute posterior predictive p-values
  out <- apply(acorrep <= acor, 1, function(x) sum(x)/dim(acorrep)[2])
  
  for (i in 1:lag.max) {
    hist(acorrep[i, ], main = paste0("Histogram Autocorrelation Lag ", i),
         xlab = "Autocorrelation")
    abline(v = acor[i], lwd = 3)
    mtext(paste0("PPP = ", round(out[i], 3)), line = -1.5, col = "red", 
          cex = 0.8, adj = 0)
  }
  
  return(round(out, 3))
}

# Mean Square Successive Difference ----

ppc.mssd <- function(object, data) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y_obs
  
  # Compute sumscores
  sumscores <- tapply(y, data$tt_obs, sum)
  sumscoresrepy <- apply(repy, 1, function(x) {
    tapply(x, data$tt_obs, sum)
  })
  
  # Compute mssd
  mssd <- sum(diff(sumscores) ^ 2)/(nT - 1)
  
  mssdrep <- apply(sumscoresrepy, 2, function(x) {
    out <- sum(diff(x) ^ 2)/(nT - 1)
    return(out)
  })
  
  # Compute posterior predictive p-values
  out <- sum(mssdrep <= mssd) / length(mssdrep)
  
  hist(mssdrep, main = "Histogram MSSD",
       xlab = "MSSD")
  abline(v = mssd, lwd = 3)
  mtext(paste0("PPP = ", round(out, 3)), line = -1.5, col = "red", 
        cex = 0.8, adj = 0)
  
  return(round(out, 3))
}

# Item score time series ----
ppc.item.ts <- function(object, data, quiet = FALSE, items = NULL) {
  
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
            c(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.025), 
              rev(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.975))),
            border = NA,
            col = rgb(1, 0, 0, 0.15))
    polygon(c(data$tt_obs[data$ii_obs == items[i]], 
              rev(data$tt_obs[data$ii_obs == items[i]])),
            c(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.25), 
              rev(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.75))),
            border = NA,
            col = rgb(1, 0, 0, 0.25))
    lines(apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.5), 
          col = "red")
    
    # Count observations out of the 95 percent of the replications
    out95 <- apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.025) <= y[data$ii_obs == items[i]] &
      y[data$ii_obs == items[i]] <= apply(repy[, data$ii_obs == items[i]], 2, quantile, probs = 0.975)
    prop_out[items[i]] <- sum(!out95) / length(out95)
    
    mtext(paste0("Prop Out = ", prop_out[items[i]]), side = 3, line = -1.5)
  }
  
  return(prop_out)
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
ppc.itcor <- function(object, data, method = c("polyserial", "pearson"), items = NULL, quiet = FALSE) {
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
    
    repdat <- apply(repy, 1, function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        out[items[i]] <- polycor::polyserial(tapply(x[data$ii_obs != items[i]], 
                                                    data$tt_obs[data$ii_obs != items[i]], 
                                                    sum), 
                                             x[data$ii_obs == items[i]])
      }
      return(out)
    })
  }
  
  if (method == "pearson") {
    polcor <- rep(NA, length(items))
    for(i in 1:length(items)) {
      polcor[items[i]] <- cor(tapply(y[data$ii_obs != items[i]],
                                     data$tt_obs[data$ii_obs != items[i]], 
                                     sum),
                              y[data$ii_obs == items[i]])
    }
    
    repdat <- apply(repy, 1, function(x){
      out  <- rep(NA, length(items))
      for (i in 1:length(items)) {
        out[items[i]] <- cor(tapply(x[data$ii_obs != items[i]],
                                    data$tt_obs[data$ii_obs != items[i]], 
                                    sum),
                             x[data$ii_obs == items[i]])
      }
      return(out)
    })
  }
  
  # Compute posterior predictive p-values
  out <- apply(repdat <= polcor, 1, function(x) sum(x)/dim(repdat)[2])
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    hist(repdat[items[i], ], main = paste0("Histogram Item-Total Correlation Item ", items[i]),
         xlab = "Item-Total Correlation")
    abline(v = polcor[items[i]], lwd = 3)
    mtext(paste0("PPP = ", round(out[items[i]], 3)), line = -1.5, col = "red", 
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
ppc.itcor2 <- function(object, data, method = c("polyserial", "pearson"), items = NULL, quiet = FALSE) {
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
    
    repdat <- apply(repy, 1, function(x){
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
    })
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
    
    repdat <- apply(repy, 1, function(x){
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
    })
  }
  
  # Compute posterior predictive p-values
  out <- apply(repdat <= polcor, 1, function(x) sum(x)/dim(repdat)[2])
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    hist(repdat[items[i], ], main = paste0("Histogram Item-Total Correlation Item ", items[i]),
         xlab = "Item-Total Correlation")
    abline(v = polcor[items[i]], lwd = 3)
    mtext(paste0("PPP = ", round(out[items[i]], 3)), line = -1.5, col = "red", 
          cex = 0.8, adj = 0)
  }
  
  # return posterior predictive p-values
  return(round(out, 3))
}


# Item-total Correlation Version 3 ----
# In this version, we fit a autoregressive model to both the item scores and 
# the rescores. Then, the residuals from these model are correlated.
ppc.itcor3 <- function(object, data, method = "pearson", items = NULL, quiet = FALSE) {
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
    
    repdat <- apply(repy, 1, function(x){
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
    })
  }
  
  # Compute posterior predictive p-values
  out <- apply(repdat <= polcor, 1, function(x) sum(x)/dim(repdat)[2])
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    hist(repdat[items[i], ], main = paste0("Histogram Item-Total Correlation Item ", items[i]),
         xlab = "Item-Total Correlation")
    abline(v = polcor[items[i]], lwd = 3)
    mtext(paste0("PPP = ", round(out[items[i]], 3)), line = -1.5, col = "red", 
          cex = 0.8, adj = 0)
  }
  
  # return posterior predictive p-values
  return(round(out, 3))
}

# Yen's Q1 unmodified ----

# Function to compute the Yen's Q1 for polytomous models
# The arguments are the estimated theta and threshold parameters, and
# the standata list.
gpcm.Q1 <- function(data, theta, beta) {
  
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  tmpx <- order(theta)
  tmpx <- cbind(tmpx, cut(seq_along(tmpx), 10, labels = FALSE))
  tmpx <- tmpx[order(tmpx[, 1]), ]
  
  group       <- tmpx[, 2]
  group_index <- tmpx[data$tt_obs, 2]
  
  thresholds  <- matrix(beta, nrow = I, ncol = M, byrow = TRUE)
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (yy in 0:M) {
    probs.array[, , yy + 1] <- P.GPCM(y     = yy, 
                                     alpha = rep(1, I), 
                                     delta = delta, 
                                     taus  = taus, 
                                     theta = theta, 
                                     M     = M)
  }
  
  E <- apply(probs.array, c(3, 2), function(x) {
    tapply(x, group, mean)
    })
  
  O <-  tapply(y, list(group_index, data$ii_obs), function(x) {
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


ppc.Q1 <- function(object, data, items = NULL, quiet = FALSE) {
  
  tmp       <- list()
  tmp$beta  <- summary(object, pars = "beta")$summary
  tmp$theta <- summary(object, pars = "theta")$summary
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  times_obs <- intersect(data$tt, data$tt_obs)
  
  
  
  
  
  
  
  
  
  
  }





