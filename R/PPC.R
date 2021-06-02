ppc.item.ts <- function(object, data, quiet = FALSE) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  nT   <- data$nT
  y    <- data$y
  for (i in 1:I) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    plot(y[((i - 1) * nT + 1):(i * nT)], type = "p", pch = 19, las = 1,
         ylim = c(min(y) - 0.5, max(y) + 0.5), xlab = "Time",
         ylab = paste0("Scores of Item ", i))
    #points(apply(repy[,((i - 1) * nT + 1):(i * nT)], 2, quantile, probs = 0.5), 
    #       col = "red", pch = 19)
    polygon(c(1:nT, rev(1:nT)),
            c(apply(repy[,((i - 1) * nT + 1):(i * nT)], 2, quantile, probs = 0.025), 
              rev(apply(repy[,((i - 1) * nT + 1):(i * nT)], 2, quantile, probs = 0.975))),
            border = NA,
            col = rgb(1, 0, 0, 0.15))
    polygon(c(1:nT, rev(1:nT)),
            c(apply(repy[,((i - 1) * nT + 1):(i * nT)], 2, quantile, probs = 0.25), 
              rev(apply(repy[,((i - 1) * nT + 1):(i * nT)], 2, quantile, probs = 0.75))),
            border = NA,
            col = rgb(1, 0, 0, 0.25))
    lines(apply(repy[,((i - 1) * nT + 1):(i * nT)], 2, quantile, probs = 0.5), 
          col = "red")
  }
}

ppc.sumscore.ts <- function(object, data) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y
  
  # Compute sumscores
  sumscores <- apply(matrix(y, ncol = I, nrow = nT), 1, sum)
  sumscoresrepy <- apply(repy, 1, function(x) {
    apply(matrix(x, ncol = I, nrow = nT), 1, sum)
  })
  
  plot(sumscores, type = "p", pch = 19, las = 1,
       ylim = c(I - 0.5, I * K + 0.5), xlab = "Time",
       ylab = "Sum Scores")
  polygon(c(1:nT, rev(1:nT)),
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
}

ppc.itcor <- function(object, data, method = c("polyserial", "pearson"), quiet = FALSE) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y
  
  if (method == "polyserial") {
    dat <- matrix(y, ncol = I, nrow = nT)
    polcor <- rep(NA, I)
    for(i in 1:I) {
      polcor[i] <- polycor::polyserial(apply(dat[, -i], 1, sum), dat[, i])
    }
    
    repdat <- apply(repy, 1, function(x){
      temp <- matrix(x, ncol = I, nrow = nT)
      out  <- rep(NA, I)
      for (i in 1:I) {
        out[i] <- polycor::polyserial(apply(temp[, -i], 1, sum), temp[, i])
      }
      return(out)
    })
  }
  
  if (method == "pearson") {
    dat <- matrix(y, ncol = I, nrow = nT)
    polcor <- rep(NA, I)
    for(i in 1:I) {
      polcor[i] <- cor(apply(dat[, -i], 1, sum), dat[, i])
    }
    
    repdat <- apply(repy, 1, function(x){
      temp <- matrix(x, ncol = I, nrow = nT)
      out  <- rep(NA, I)
      for (i in 1:I) {
        out[i] <- cor(apply(temp[, -i], 1, sum), temp[, i])
      }
      return(out)
    })
  }
  
  for (i in 1:I) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    hist(repdat[i, ], main = paste0("Histogram Item-Total Correlation Item ", i),
         xlab = "Item-Total Correlation")
    abline(v = polcor[i], lwd = 3)
  }
  
  # return posterior predictive p-values
  out <- apply(repdat <= polcor, 1, function(x) sum(x)/dim(repdat)[2])
  names(out) <- paste0("Item_", 1:I)
  return(round(out, 3))
}

ppc.itcor2 <- function(object, data, method = c("polyserial", "pearson"), quiet = FALSE) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y
  
  if (method == "polyserial") {
    dat <- matrix(y, ncol = I, nrow = nT)
    polcor <- rep(NA, I)
    for(i in 1:I) {
      sums <- apply(dat[, -i], 1, sum) 
      res <- rstandard(lm(sums[-1] ~ head(sums, -1)))
      polcor[i] <- polycor::polyserial(res, dat[-1, i])
    }
    
    repdat <- apply(repy, 1, function(x){
      temp <- matrix(x, ncol = I, nrow = nT)
      out  <- rep(NA, I)
      for (i in 1:I) {
        sums <- apply(temp[, -i], 1, sum)
        res <- rstandard(lm(sums[-1] ~ head(sums, -1)))
        out[i] <- polycor::polyserial(res, temp[-1, i])
      }
      return(out)
    })
  }

  if (method == "pearson") {
    dat <- matrix(y, ncol = I, nrow = nT)
    polcor <- rep(NA, I)
    for(i in 1:I) {
      sums <- apply(dat[, -i], 1, sum) 
      res <- rstandard(lm(sums[-1] ~ head(sums, -1)))
      polcor[i] <- cor(res, dat[-1, i])
    }
    
    repdat <- apply(repy, 1, function(x){
      temp <- matrix(x, ncol = I, nrow = nT)
      out  <- rep(NA, I)
      for (i in 1:I) {
        sums <- apply(temp[, -i], 1, sum)
        res <- rstandard(lm(sums[-1] ~ head(sums, -1)))
        out[i] <- cor(res, temp[-1, i])
      }
      return(out)
    })
  }
  
  for (i in 1:I) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    hist(repdat[i, ], main = paste0("Histogram Item-Total Correlation Item ", i),
         xlab = "Item-Total Correlation")
    abline(v = polcor[i], lwd = 3)
  }
  
  # return posterior predictive p-values
  out <- apply(repdat <= polcor, 1, function(x) sum(x)/dim(repdat)[2])
  names(out) <- paste0("Item_", 1:I)
  return(round(out, 3))
}

ppc.itcor3 <- function(object, data, method = "pearson", quiet = FALSE) {
  repy <- extract(object)[["rep_y"]]
  I    <- data$I
  K    <- data$K
  nT   <- data$nT
  y    <- data$y
  
  if (method == "pearson") {
    dat <- matrix(y, ncol = I, nrow = nT)
    polcor <- rep(NA, I)
    for(i in 1:I) {
      sums <- apply(dat[, -i], 1, sum) 
      res <- rstandard(lm(sums[-1] ~ head(sums, -1)))
      ires <- rstandard(lm(dat[-1, i] ~ head(dat[, i], -1)))
      polcor[i] <- cor(res, ires)
    }
    
    repdat <- apply(repy, 1, function(x){
      temp <- matrix(x, ncol = I, nrow = nT)
      out  <- rep(NA, I)
      for (i in 1:I) {
        sums <- apply(temp[, -i], 1, sum)
        res <- rstandard(lm(sums[-1] ~ head(sums, -1)))
        ires <- rstandard(lm(temp[-1, i] ~ head(temp[, i], -1)))
        out[i] <- cor(res, ires)
      }
      return(out)
    })
  }
  
  for (i in 1:I) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    hist(repdat[i, ], main = paste0("Histogram Item-Total Correlation Item ", i),
         xlab = "Item-Total Correlation")
    abline(v = polcor[i], lwd = 3)
  }
  
  # return posterior predictive p-values
  out <- apply(repdat <= polcor, 1, function(x) sum(x)/dim(repdat)[2])
  names(out) <- paste0("Item_", 1:I)
  return(round(out, 3))
}


ppc.itcor(fit, standata, method = "polyserial", quiet =TRUE)
ppc.itcor2(fit, standata, method = "polyserial", quiet =TRUE)
ppc.itcor3(fit, standata, quiet =TRUE)





