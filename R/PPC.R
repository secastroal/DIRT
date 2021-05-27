repy <- extract(fit)[["rep_y"]]
pdf(file = "ppc_pcm_data3.pdf")
for (i in 1:I) {
  #invisible(readline(prompt="Press [enter] to continue"))
  plot(c(responses)[((i - 1) * nT + 1):(i * nT)], type = "p", pch = 19, las = 1,
       ylim = c(min(c(responses)) - 0.5, max(c(responses)) + 0.5), xlab = "Time",
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


sumscores <- apply(matrix(c(responses), ncol = I, nrow = nT), 1, sum)
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


# item-total corr
dat <- matrix(c(responses), ncol = I, nrow = nT)
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
par(mfrow = c(3,2))
for (i in 1:I) {
  hist(repdat[i, ])
  abline(v = polcor[i], lwd = 3)
}

#pearson
dat <- matrix(c(responses), ncol = I, nrow = nT)
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
par(mfrow = c(3,2))
for (i in 1:I) {
hist(repdat[i, ])
abline(v = polcor[i], lwd = 3)
}

#length(repdat[1, repdat[1, ] <= polcor[1]])/dim(repdat)[2]

# Polyserial correlation variation with AR residuals of the sumscores

dat <- matrix(c(responses), ncol = I, nrow = nT)
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

par(mfrow = c(3,2))
for (i in 1:I) {
  hist(repdat[i, ])
  abline(v = polcor[i], lwd = 3)
}


#length(repdat[6, repdat[6, ] <= polcor[6]])/dim(repdat)[2]

# Correlate residuals of AR with sumscores and item

dat <- matrix(c(responses), ncol = I, nrow = nT)
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

par(mfrow = c(3,2))
for (i in 1:I) {
  hist(repdat[i, ])
  abline(v = polcor[i], lwd = 3)
}

#length(repdat[1, repdat[1, ] <= polcor[1]])/dim(repdat)[2]
dev.off()




