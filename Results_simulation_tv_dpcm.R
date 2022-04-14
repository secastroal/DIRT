# Summary Simulations

# This file plots the results form the simulation study of the TV-DPCM model.

# 0.0 Prepare environment ----
library(plyr)

# Recreate matrix of conditions 
N.timepoints <- c(100, 200, 300, 500)                # Number of timepoints
N.items      <- c(3, 6)                  # Number of items 
S.lambda     <- c(0, 0.25, 0.5)          # Size of the autoregressive effect
M.prop       <- c(0, 0.3)            # Proportion of missing values

Cond        <- expand.grid(N.timepoints, N.items, S.lambda, M.prop)
names(Cond) <- c("nT", "I", "lambda", "NAprop")

#rm(N.timepoints, N.items, S.lambda, M.prop)

# Colors for plots
# If bnw = TRUE, plots are printed in grayscale.
bnw = TRUE #black and white
if (bnw) {
  color_plot <- gray((0:2)/2)
}else{
  color_plot <- rainbow(2)
}

# Axis and legend statements
# As we are plotting 6 plots in one figure, these matrices indicate for which of 
# these plots the y axis, the x axis, or the legend must be printed. 
yaxis <- matrix(rep(c(TRUE, FALSE), times = c(2, 4)), ncol = 3)
xaxis <- matrix(rep(c(FALSE, TRUE), each = 3), ncol = 3, byrow = TRUE)
laxis <- matrix(FALSE, ncol = 3, nrow = 2)
laxis[1, 3] <- TRUE

# 1.0 Read output files ----

# Read output files into R
files <- paste0(getwd(), "/Simulation/Sim_TV_DPCM_cond_", 1:nrow(Cond), ".txt")

data <- lapply(files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame
results <- ldply(data, data.frame)

rm(data, files)

Cond.ext <- Cond[results$cond, ]

# Plot percentage of analysys that converged
tmp <- 100 - tapply(results$corrupt, Cond.ext, mean) * 100

# Define plotting parameters.
par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)

for (i in 1:2) {
  for(l in 1:3) {
    
    matplot(1:4, tmp[, i, l,], type = "l", col = color_plot, 
            ylim = c(-0.5, 100.5), xlab = "", 
             ylab = "", xaxt = "n", yaxt = "n", lwd = 2, lty = 1)
    if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
    if (xaxis[i, l]) {axis(1, at = 1:4, labels = N.timepoints, cex.axis = 2)}
    if (laxis[i, l]) {legend("topright", legend = c("0%", "30%"),
                             col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                             seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
  }
}
mtext("Percentage of Divergent Analyses", 2, outer = TRUE, line = 4, cex = 1.8)
mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3/4, las = 1)
mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1/4, las = 1)
mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1/6)
mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3/6)
mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5/6)

dev.off()

# Plot percentage of analysys that converged
tmp <- tapply(results$beta.cover, Cond.ext, mean)

# Define plotting parameters.
par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)

for (i in 1:2) {
  for(l in 1:3) {
    
    matplot(1:4, tmp[, i, l,], type = "l", col = color_plot, 
            xlab = "", 
            ylab = "", xaxt = "n", yaxt = "n", lwd = 2, lty = 1)
    if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
    if (xaxis[i, l]) {axis(1, at = 1:4, labels = N.timepoints, cex.axis = 2)}
    if (laxis[i, l]) {legend("topright", legend = c("0%", "30%"),
                             col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                             seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
  }
}
mtext("Percentage of Divergent Analyses", 2, outer = TRUE, line = 4, cex = 1.8)
mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3/4, las = 1)
mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1/4, las = 1)
mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1/6)
mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3/6)
mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5/6)

dev.off()

tapply(results$corrupt, results$cond, sum)
tapply(results$efficiency, results$cond, sum)
tapply(results$run.time, results$cond, mean)

tapply(results$ndiv, results$cond, function(x) sum(x != 0))
tapply(results$nbfmi, results$cond, function(x) sum(x != 0))
tapply(results$nRhat, results$cond, function(x) sum(x != 0))
tapply(results$ntree, results$cond, function(x) sum(x != 0))
tapply(results$nbulk, results$cond, function(x) sum(x != 0))
tapply(results$ntail, results$cond, function(x) sum(x != 0))

# Exclude divergent analyses
results_conv <- results[results$corrupt == 0, ]

tapply(results_conv$ndiv, results_conv$cond, function(x) sum(x != 0))
tapply(results_conv$nbfmi, results_conv$cond, function(x) sum(x != 0))
tapply(results_conv$nRhat, results_conv$cond, function(x) sum(x != 0))
tapply(results_conv$ntree, results_conv$cond, function(x) sum(x != 0))
tapply(results_conv$nbulk, results_conv$cond, function(x) sum(x != 0))
tapply(results_conv$ntail, results_conv$cond, function(x) sum(x != 0))

#beta
tapply(results_conv$beta.cor, results_conv$cond, mean)
tapply(results_conv$beta.rmse, results_conv$cond, mean)
tapply(results_conv$beta.cover, results_conv$cond, mean)
tapply(results_conv$beta.CI, results_conv$cond, mean)

#theta
tapply(results_conv$theta.cor, results_conv$cond, mean)
tapply(results_conv$theta.bias, results_conv$cond, mean)
tapply(results_conv$theta.cover, results_conv$cond, mean)
tapply(results_conv$theta.CI, results_conv$cond, mean)

#attractor
tapply(results_conv$attra.cor, results_conv$cond, mean)
tapply(results_conv$attra.bias, results_conv$cond, mean)
tapply(results_conv$attra.cover, results_conv$cond, mean)
tapply(results_conv$attra.CI, results_conv$cond, mean)

#lambda
tapply(results_conv$lambda.bias, results_conv$cond, mean)
tapply(results_conv$lambda.rbias, results_conv$cond, mean)
tapply(results_conv$lambda.cover, results_conv$cond, mean)
tapply(results_conv$lambda.CI, results_conv$cond, mean)

#sigma2
tapply(results_conv$sigma2.bias, results_conv$cond, mean)
tapply(results_conv$sigma2.rbias, results_conv$cond, mean)
tapply(results_conv$sigma2.cover, results_conv$cond, mean)
tapply(results_conv$sigma2.CI, results_conv$cond, mean)

#pvar
tapply(results_conv$pvar.bias, results_conv$cond, mean)
tapply(results_conv$pvar.rbias, results_conv$cond, mean)
tapply(results_conv$pvar.cover, results_conv$cond, mean)
tapply(results_conv$pvar.CI, results_conv$cond, mean)

