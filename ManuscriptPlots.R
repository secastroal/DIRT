# Additional Plots for the Manuscript

# This script creates plots for the manuscript, which are unrelated with 
# the results of the simulation or the empirical example.

source("R/genTVDPCM.R")
source("R/IRT_models.R")

# B-splines ----
set.seed(1234)
x <- seq(0, 10, by = 0.1)
X <- t(bs(x, df = 10, degree = 3, intercept = TRUE))

num_data  <- length(x)
num_basis <- nrow(X)

# Define intercept and coefficients to generate Y based on X and the b-splines.
a  <- rnorm(num_basis, 0, 1) # coefficients of B-splines
a[c(1, num_basis)] <- a[c(1, num_basis)] * 0.1

# Compute predicted values of Y and add random error to generate Y. 
Y_true <- 5 + as.vector(a %*% X)    # generating the output
Y      <- Y_true + rnorm(length(x),0,.5) # adding noise

fit1 <- lm(Y ~ t(X) + 0)

# Underfit
X2 <- t(bs(x, df = 5, degree = 3, intercept = TRUE))

fit2 <- lm(Y ~ t(X2) + 0)

# Overfit
X3 <- t(bs(x, df = 30, degree = 3, intercept = TRUE))

fit3 <- lm(Y ~ t(X3) + 0)

pdf(file = "Figures/B-Splines.pdf", height = 4)
par(mfrow = c(1, 3), mar = c(0.6, 0.6, 0.6, 0.6), oma = c(4, 4, 2, 1), xpd = '')

plot(x, Y, pch = 20, ylim = c(0, 7), las = 1, xlab = "X", col = "darkgray")
lines(x, predict(fit2, data.frame(x)), lwd = 2)
matplot(x, t(X2), type = "l", add = TRUE, lty = 2, 
        col = gray((0:nrow(X2))/ (nrow(X2) + 3)))
#matplot(t(X2 * fit2$coefficients), type = "l", add = TRUE)

plot(x, Y, pch = 20, ylim = c(0, 7), xlab = "X", yaxt = 'n', ylab = "",
     col = "darkgray")
lines(x, predict(fit1, data.frame(x)), lwd = 2)
matplot(x, t(X), type = "l", add = TRUE, lty = 2, 
        col = gray((0:nrow(X))/ (nrow(X) + 3)))
#matplot(t(X * fit1$coefficients), type = "l", add = TRUE)

plot(x, Y, pch = 20, ylim = c(0, 7), xlab = "X", yaxt = 'n', ylab = "", 
     col = "darkgray")
lines(x, predict(fit3, data.frame(x)), lwd = 2)
matplot(x, t(X3), type = "l", add = TRUE, lty = 2, 
        col = gray((0:nrow(X3))/ (nrow(X3) + 3)))
#matplot(t(X3 * fit3$coefficients), type = "l", add = TRUE)

mtext("A", side = 3, line = 0.2, outer = TRUE, at = 0.025, font = 2, cex = 1.3)
mtext("B", side = 3, line = 0.2, outer = TRUE, at = 0.355, font = 2, cex = 1.3)
mtext("C", side = 3, line = 0.2, outer = TRUE, at = 0.685, font = 2, cex = 1.3)

dev.off()


# True Trend -----

gendata <- gen.TVDPCM(200, 10, 5, 
                      pop.param = list(lambda = 0), 
                      seed = 1234,
                      FUN = "sinusoidal", maxAbsValue = 1)

pdf(file = "Figures/Trend.pdf", height = 4)
plot(gendata$theta.gen, type = "l", las = 1, ylim = c(-2.7, 3.7), xaxt = 'n',
     xlab = "", ylab = expression(theta[t]), lwd = 0.5, col = "darkgray")
lines(gendata$attractor.gen, lwd = 2)
mtext("Time", side = 1, line = 1)
dev.off()
