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
# The arguments are the estimated theta, threshold, and discrimination 
# parameters, the number of items, the number of response categories,
# the grouping variable for each time point, and a grouping index for
# the data in long format.
gpcm.Q1 <- function(y, theta, thresholds, alpha, I, K, group, group_index) {
  
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
  
  # Create array to store the computed Yen's Q1.
  discrepancy <- array(NA, dim = c(nrow(repy), I, 2))
  
  for (r in 1:nrow(repy)) {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ][times_obs]
    
    # Create groups for Yen's Q1
    tmpx <- order(theta)
    tmpx <- cbind(tmpx, cut(seq_along(tmpx), 10, labels = FALSE))
    tmpx <- tmpx[order(tmpx[, 1]), ]
    
    group       <- tmpx[, 2]
    group_index <- tmpx[data$tt_obs, 2]
    rm(tmpx)
    
    # Yen's Q1 for the observed scores
    discrepancy[r, , 1] <- gpcm.Q1(y           = y, 
                                   theta       = theta, 
                                   thresholds  = thresholds, 
                                   alpha       = rep(1, I), 
                                   I           = I, 
                                   K           = K, 
                                   group       = group,
                                   group_index = group_index)
    # Yen's Q1 for the i-th replicated scores
    discrepancy[r, , 2] <- gpcm.Q1(y           = repy[r, ], 
                                   theta       = theta, 
                                   thresholds  = thresholds, 
                                   alpha       = rep(1, I), 
                                   I           = I, 
                                   K           = K, 
                                   group       = group,
                                   group_index = group_index)
  }
  rm(r)
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] > discrepancy[, , 1], 2, mean)
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    plot(discrepancy[, items[i], 1], discrepancy[, items[i], 2], las = 1,
         main = paste0("Scatterplot Yen's Q1 of item ", items[i]),
         ylab = expression(paste("Yen's ", Q[1], "(", y^rep, ";", Theta, ")")),
         xlab = expression(paste("Yen's ", Q[1], "(y;", Theta, ")")))
    abline(a= 0, b = 1, col = "red")
    mtext(paste0("PPP = ", round(out[items[i]], 3)), line = -1.5, col = "red", 
          cex = 0.8, adj = 0)
  }
  
  return(round(out, 3))
  }

# Yen's Q1 modified ----
# This is an alternative version of Yen's Q1 in which the groups are not made by
# ordering the estimated latent values. Instead, the groups are simply defined 
# by the time in which the observation were made. The purpose of this is to take 
# the time into account.

ppc.Q1.alt <- function(object, data, items = NULL, quiet = FALSE) {
  
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
  
  # Create array to store the computed Yen's Q1.
  discrepancy <- array(NA, dim = c(nrow(repy), I, 2))
  
  for (r in 1:nrow(repy)) {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ][times_obs]
    
    # Create groups for Yen's Q1 given the time order
    group       <- cut(seq_along(theta), 10, labels = FALSE)
    group_index <- group[data$tt_obs]
    
    # Yen's Q1 for the observed scores
    discrepancy[r, , 1] <- gpcm.Q1(y           = y, 
                                   theta       = theta, 
                                   thresholds  = thresholds, 
                                   alpha       = rep(1, I), 
                                   I           = I, 
                                   K           = K, 
                                   group       = group,
                                   group_index = group_index)
    # Yen's Q1 for the i-th replicated scores
    discrepancy[r, , 2] <- gpcm.Q1(y           = repy[r, ], 
                                   theta       = theta, 
                                   thresholds  = thresholds, 
                                   alpha       = rep(1, I), 
                                   I           = I, 
                                   K           = K, 
                                   group       = group,
                                   group_index = group_index)
  }
  rm(r)
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] > discrepancy[, , 1], 2, mean)
  names(out) <- paste0("Item_", items)
  
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    plot(discrepancy[, items[i], 1], discrepancy[, items[i], 2], las = 1,
         main = paste0("Scatterplot Alt. Yen's Q1 of item ", items[i]),
         ylab = expression(paste("Alt. Yen's ", Q[1], "(", y^rep, ";", Theta, ")")),
         xlab = expression(paste("Alt. Yen's ", Q[1], "(y;", Theta, ")")))
    abline(a= 0, b = 1, col = "red")
    mtext(paste0("PPP = ", round(out[items[i]], 3)), line = -1.5, col = "red", 
          cex = 0.8, adj = 0)
  }
  
  return(round(out, 3))
}

# Yen's Q3 ----

# First we create a function to compute Yen's Q3. The function requires the 
# vector of responses y, the estimated theta, thresholds, and discrimination
# parameters, the number of time points, the number of items, the number of 
# response categories, and two index variables t_index for time and i_index 
# for items given that the data is in long format.

gpcm.Q3 <- function(y, theta, thresholds, alpha, nT, I, K, t_index, i_index) {
  
  M <- K - 1
  
  # Restructure responses in a matrix
  Y <- matrix(NA, nT, I)
  
  for (t in 1:nT) {
    for (i in 1:I) {
      Y[t, i] <- y[t_index == t & i_index == i]
    }
  }
  
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
  
  D <- Y - E
  
  CorD <- cor(D, use = "pairwise.complete.obs")
  q3 <- CorD[lower.tri(CorD)]
  names(q3) <- apply(which(lower.tri(CorD), arr.ind = TRUE), 1, function(x)
    paste0("(", paste(x, collapse = ","), ")"))
  
  return(q3)
}

ppc.Q3 <- function(object, data, scatterplots = FALSE) {
  
  require(scatterpie)
  
  betasamples  <- extract(object)[["beta"]]
  thetasamples <- extract(object)[["theta"]]
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
  nT   <- data$nT
  y    <- data$y_obs
  
  # Create array to store the computed Yen's Q3.
  discrepancy <- array(NA, dim = c(nrow(repy), (I * (I - 1))/2, 2))
  
  for (r in 1:nrow(repy)) {
    # Get estimated parameters for the i-th iteration.
    thresholds <- betasamples[r, , ]
    theta      <- thetasamples[r, ]
    
    # Yen's Q3 for the observed scores
    discrepancy[r, , 1] <- gpcm.Q3(y          = y,
                                   theta      = theta,
                                   thresholds = thresholds,
                                   alpha      = rep(1, I),
                                   nT         = nT,
                                   I          = I,
                                   K          = K,
                                   t_index    = data$tt_obs,
                                   i_index    = data$ii_obs)
    
    # Yen's Q3 for the i-th replicated scores
    discrepancy[r, , 2] <- gpcm.Q3(y          = repy[r, ],
                                   theta      = theta,
                                   thresholds = thresholds,
                                   alpha      = rep(1, I),
                                   nT         = nT,
                                   I          = I,
                                   K          = K,
                                   t_index    = data$tt_obs,
                                   i_index    = data$ii_obs)
  }
  rm(r)
  
  # Compute posterior predictive p-values
  out  <- apply(discrepancy[, , 2] > discrepancy[, , 1], 2, mean)
  names(out) <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
    paste0("ppp-q3(", paste(x, collapse = ","), ")"))
  
  if (scatterplots) {
    for (i in 1:15) {
      plot(discrepancy[, i, 1], discrepancy[, i, 2], las = 1,
           main = paste0("Scatterplot Yen's Q3 of items ", substring(names(out)[i], 7)),
           ylab = expression(paste("Yen's ", Q[3], "(", y^rep, ";", Theta, ")")),
           xlab = expression(paste("Yen's ", Q[3], "(y;", Theta, ")")))
      abline(a= 0, b = 1, col = "red")
      mtext(paste0("PPP = ", round(out[i], 3)), line = -1.5, col = "red", 
            cex = 0.8, adj = 0)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out < 0.05 | tmp$out > 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                             cols=c("out","out2"), color = NA,
                             show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(round(out, 3))
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
odds.ratio <- function(y, nT, I, K, t_index, i_index) {
  
  # Restructure responses in a matrix
  Y <- matrix(NA, nT, I)
  
  for (t in 1:nT) {
    for (i in 1:I) {
      Y[t, i] <- y[t_index == t & i_index == i]
    }
  }
  rm(t, i)
  
  # Vector to store Odds ratio
  or <- rep(NA, (I * (I - 1))/2)
  
  count <- 1
  for (i in 1:(I - 1)) {
    j <- i + 1
    while (i<j) {
      n    <- rep(NA, 4)
      n[1] <- sum(Y[, i] == 1 & Y[, j] == 1, na.rm = TRUE)
      n[2] <- sum(Y[, i] == 0 & Y[, j] == 0, na.rm = TRUE)
      n[3] <- sum(Y[, i] == 1 & Y[, j] == 0, na.rm = TRUE)
      n[4] <- sum(Y[, i] == 0 & Y[, j] == 1, na.rm = TRUE)
      if (any(n[3:4] == 0)) {n <- n + 0.5}
      or[count] <- (n[1]*n[2])/(n[3]*n[4])
      names(or)[count] <- paste0('OR(', j, ',', i, ")")
      count <- count + 1
      if (j == I) {
        break
      } else {
        j <- j + 1
      }
    }
  }
  
  return(or)
}

ppc.OR <- function(object, data, cutoff = NULL, histograms = FALSE) {
  
  require(scatterpie)
  
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  M    <- data$K - 1
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
                      nT = nT, 
                      I  = I, 
                      K  = K, 
                      t_index = data$tt_obs, 
                      i_index = data$ii_obs)
  
  
  orrep <- apply(repy, 1, function(x) {
    repyd <- dicho(y = x, 
                   I = I, 
                   K = K, 
                   i_index = data$ii_obs, 
                   cutoff  = cutoff)
    out   <- odds.ratio(y  = repyd, 
                        nT = nT, 
                        I  = I, 
                        K  = K, 
                        t_index = data$tt_obs, 
                        i_index = data$ii_obs)
    return(out)
  })
  
  # Compute posterior predictive p-values
  out <- apply(orrep <= or, 1, function(x) sum(x)/dim(orrep)[2])
  
  
  if (histograms) {
    for (i in 1:length(or)) {
      hist(orrep[i, ], main = paste0("Histogram Odd Ratio of Items ", 
                                     substring(names(out)[i], 3)),
           xlab = "Odd Ratio")
      abline(v = or[i], lwd = 3)
      mtext(paste0("PPP = ", round(out[i], 3)), line = -1.5, col = "red", 
            cex = 0.8, adj = 0)
    }
  }
  
  tmp <- data.frame(which(lower.tri(diag(I)), arr.ind = TRUE), out)
  tmp$out2   <- 1 - out
  tmp$radius <- ifelse(tmp$out < 0.05 | tmp$out > 0.95, 0.4, 0.2)
  
  sp <- ggplot() + geom_scatterpie(aes(x=row, y=col, r = radius), data = tmp,
                                   cols=c("out","out2"), color = NA,
                                   show.legend = FALSE) + coord_equal()  +
    scale_fill_manual(values = c("black", gray(0.9))) +
    labs(y = "Item", x = "Item") + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  print(sp)
  return(round(out, 3))
}


# rm(list = setdiff(ls(), c(lsf.str(), "object", "data", "fit")))

