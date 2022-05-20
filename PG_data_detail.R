
for (i in 1:5) {
  for (j in 1:3) {
    matplot(PG_data[PG_data$phase == i, mu[-j]], 
            type = "l", lty = 1, las = 1,
            col = c("blue", "red", "darkgreen")[-j])
  }
}

for (i in 1:5) {
  print(cor(PG_data[PG_data$phase == i, mu]))
}


out <- c()

window <- 100

for (i in window:nrow(PG_data)) {
  tmp <- cor(PG_data[(i - (window - 1)):i, mu])[lower.tri(diag(3))] 
  out <- rbind(out, tmp)
}

plot(out[, 1], type = "l", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Irritated - Restless")
plot(out[, 2], type = "l", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Irritated - Agitated")
plot(out[, 3], type = "l", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Restless - Agitated")


out <- c()

for (i in 30:nrow(PG_data)) {
  tmp <- apply(PG_data[(i - 29):i, mu], 2, mean)
  out <- rbind(out, tmp)
}

plot(out[, 1], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Mean", ylim = c(1, 4),
     main = "Irritated")
plot(out[, 2], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Mean", ylim = c(1, 4),
     main = "Restless")
plot(out[, 3], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Mean", ylim = c(1, 4),
     main = "Agitated")


out <- c()

for (i in 30:nrow(PG_data)) {
  tmp <- apply(PG_data[(i - 29):i, mu], 2, var)
  out <- rbind(out, tmp)
}

plot(out[, 1], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Variance", ylim = c(0, 4),
     main = "Irritated")
plot(out[, 2], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Variance", ylim = c(0, 4),
     main = "Restless")
plot(out[, 3], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Variance", ylim = c(0, 4),
     main = "Agitated")

for (i in 1:5) {
  print(apply(PG_data[PG_data$phase == i, mu], 2, table))
}

for (i in 1:5) {
  print(sort(table(apply(PG_data[PG_data$phase == i, mu], 1, paste, collapse = ""))))
}


out <-  c()

for (i in 30:nrow(PG_data)) {
  tmp <- mean(apply(PG_data[(i - 29):i, mu], 1, mean))
  out <- c(out, tmp)
}

plot(out, type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Mean", ylim = c(1, 4),
     main = "Sumscore")


out <-  c()

for (i in 30:nrow(PG_data)) {
  tmp <- var(apply(PG_data[(i - 29):i, mu], 1, mean))
  out <- c(out, tmp)
}

plot(out, type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Variance", ylim = c(0, 2),
     main = "Sumscore")

out <- c()

window <- 100

for (i in window:nrow(responses)) {
  tmp <- cor(responses[(i - (window - 1)):i, ])[lower.tri(diag(ncol(responses)))] 
  out <- rbind(out, tmp)
}

plot(out[, 1], type = "l", ylim = c(0, 1))
plot(out[, 5], type = "l", ylim = c(-1, 1))
plot(out[, 3], type = "l", ylim = c(0, 1))
plot(out[, 4], type = "l", ylim = c(0, 1))
plot(out[, 5], type = "l", ylim = c(0, 1))

plot(out[, 1], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Irritated - Restless")
plot(out[, 2], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Irritated - Agitated")
plot(out[, 3], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Restless - Agitated")

# Restless seems to be a problematic item.
# Simulate less extreme cases of the change of meaning
# maybe the problem is collinearity... collinearity in IRT
# Our items are too correlated.
# See the Item fit for the other scales.
# Choose other items based on theory.
