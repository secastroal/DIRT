
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

for (i in 30:nrow(PG_data)) {
  tmp <- cor(PG_data[(i - 29):i, mu])[lower.tri(diag(3))] 
  out <- rbind(out, tmp)
}

plot(out[, 1], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Irritated - Restless")
plot(out[, 2], type = "p", col = PG_data$phase[1:1444], las = 1, 
     xlab = "Beep", ylab = "Correlation", ylim = c(0, 1),
     main = "Correlation Irritated - Agitated")
plot(out[, 3], type = "p", col = PG_data$phase[1:1444], las = 1, 
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
