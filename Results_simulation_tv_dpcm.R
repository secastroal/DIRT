# Summary Simulations

# This file plots the results from the simulation study of the TV-DPCM model.

# 0.0 Prepare environment ----
library(plyr)

# Recreate matrix of conditions 
N_timepoints <- c(100, 200, 300, 500)                # Number of timepoints
N_items      <- c(3, 6)                  # Number of items 
S_lambda     <- c(0, 0.25, 0.5)          # Size of the autoregressive effect
M_prop       <- c(0, 0.3)            # Proportion of missing values

Cond        <- expand.grid(N_timepoints, N_items, S_lambda, M_prop)
names(Cond) <- c("nT", "I", "lambda", "NAprop")

#rm(N.timepoints, N.items, S.lambda, M.prop)

# Colors for plots
# If bnw = TRUE, plots are printed in grayscale.
bnw = TRUE #black and white
if (bnw) {
  color_plot <- gray((0:1)/2)
} else {
  color_plot <- rainbow(2)
}

# Axis and legend statements
# As we are plotting 6 plots in one figure, these matrices indicate for which of 
# these plots the y axis, the x axis, or the legend must be printed. 
yaxis <- matrix(rep(c(TRUE, FALSE), times = c(2, 4)), ncol = 3)
xaxis <- matrix(rep(c(FALSE, TRUE), each = 3), ncol = 3, byrow = TRUE)
laxis <- matrix(FALSE, ncol = 3, nrow = 2)
# laxis[1, 3] <- TRUE

# 1.0 Read output files ----

# Read output files into R
files <- paste0(getwd(), "/Simulation/Sim_TV_DPCM_cond_", 1:nrow(Cond), ".txt")

data <- lapply(files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame
results <- ldply(data, data.frame)

rm(data, files)

# Extend condition matrix to match the total number of analyses performed.
Cond_ext <- Cond[results$cond, ]

# 2.0 Check Convergence per Condition ----
# For this, we take the variable corrupt, which is a variable indicator that
# identifies the analyses that had warnings about convergence. In other words,
# analyses that had Rhat statistics larger than 1.05, divergent transitions,
# or too low BMIF. 
# Plot percentage of analysys that converged
tmp <- 100 - tapply(results$corrupt, Cond_ext, mean) * 100

pdf(file = "Figures/Convergent_TVDPCM.pdf", width = 15)
# Define plotting parameters.
par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 6), xpd = NA)

for (i in 1:2) {
  for(l in 1:3) {
    matplot(1:4, tmp[, i, l,], type = "l", col = color_plot, 
            ylim = c(79, 100.5), xlab = "", 
             ylab = "", xaxt = "n", yaxt = "n", lwd = 2, lty = 1)
    if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
    if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
    if (laxis[i, l]) {legend("topright", legend = c("0%", "30%"),
                             col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                             seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
  }
}
mtext("Percentage of Convergent Analyses", 2, outer = TRUE, line = 4, cex = 1.8)
mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3/4, las = 1)
mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1/4, las = 1)
mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1/6)
mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3/6)
mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5/6)
rm(i, l, tmp)
dev.off()

# 3.0 Plotting Accuracy Measures ----
# In this section we plot the different accuracy measure of interest per 
# condition and for each parameter. However, before plotting these summaries, 
# we exclude the analyses that had convergence issues. 
# Filter the results as well as the Cond_ext to exclude divergent analyses.
results_conv <- results[results$corrupt == 0, ]
Cond_ext     <- Cond_ext[results$corrupt == 0, ]

# For loop to create all plots.
# Extract columns of interest
col_index <- grep("beta|theta|attra|lambda|pvar|sigma2", names(results_conv))

# Character vectors to define file name and y axis label.
par_name <- rep(c("Thresholds", "Theta", "Attractor", 
                  "Autoregression", "ProcessVar", "InnoVar"), each = 6)
lab_name <- rep(c("Thresholds", "State Disposition", "Attractor", 
                  "Autoregression", "Process Variance", 
                  "Innovation Variance"), each = 6)
sta_name <- sub(".*\\.","", names(results_conv)[col_index])
sta_name <- sub("cor", "Correlation", sta_name)
sta_name <- sub("rbias2", "Relative Bias 2.0", sta_name)
sta_name <- sub("rbias", "Relative Bias", sta_name)
sta_name <- sub("abbias", "Absolute Bias", sta_name)
sta_name <- sub("bias", "Bias", sta_name)
sta_name <- sub("rmse", "RMSE", sta_name)
sta_name <- sub("cover", "C.I. Coverage", sta_name)
sta_name <- sub("CI", "C.I. Width", sta_name)

plot_points <- cbind(1:4, 1:4 + 0.1)

for (p in 1:length(par_name)) {
  # Export Plot to PDF
  pdf(file = paste0("Figures/", par_name[p], "_",
                    gsub(" ", "_", sta_name[p]), ".pdf"), width = 15 )
  
  # Define plotting parameters.
  par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 8, 5, 6), xpd = NA)
  
  # Compute mean of the statistic across conditions
  tmp   <- tapply(results_conv[, col_index[p]], Cond_ext, mean, na.rm = TRUE)
  tmp_l <- tapply(results_conv[, col_index[p]], Cond_ext, quantile, 
                  probs = 0.25, na.rm = TRUE)
  tmp_u <- tapply(results_conv[, col_index[p]], Cond_ext, quantile, 
                  probs = 0.75, na.rm = TRUE)
  
  # Define lower and upper limits for the plots.
  lylim <- round(min(tmp_l, na.rm = TRUE) - 0.01, 2)
  uylim <- round(max(tmp_u, na.rm = TRUE) + 0.01, 2)
  if (sta_name[p] == "C.I. Coverage") {lylim <- 0}
  if (sta_name[p] == "Correlation" | sta_name[p] == "C.I. Coverage") {uylim <- 1}
  
  for (i in 1:2) {
    for (l in 1:3) {
      matplot(plot_points, tmp[, i, l,], type = "b", 
              col = color_plot, ylim = c(lylim, uylim), xlab = "", 
              ylab = "", xaxt = "n", yaxt = "n", lwd = 2, pch = 19, lty = 1)
      if (sta_name[p] == "Relative Bias 2.0") {
        abline(h = 1, lty = 2, col = gray(0.9), xpd = FALSE)
      } else {
        abline(h = 0, lty = 2, col = gray(0.9), xpd = FALSE)  
      }
      if (!(sta_name[p] == "C.I. Coverage" & 
          (par_name[p] %in% c("Autoregression", "ProcessVar", "InnoVar")))) {
        segments(x0 = c(plot_points),
                 y0 = c(tmp_l[, i, l, ]),
                 y1 = c(tmp_u[, i, l, ]),
                 lwd = 1, lty = 3, col = rep(color_plot, each = 4))
      }
      if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
      if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
      if (laxis[i, l]) {legend("topright", legend = c("0%", "30%"),
                               col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                               seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
    }
  }
  rm(i, l, tmp, tmp_l, tmp_u, lylim, uylim)
  
  mtext(paste("Average", lab_name[p], sta_name[p]), 2, outer = TRUE, 
        line = 5, cex = 1.8)
  mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
  mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
  mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
  mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
  mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
  mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)
  
  dev.off()
}
rm(p, plot_points, col_index, lab_name, par_name, sta_name)

# Additional plots of the median BULK and TAIL ESS.
# We run a very small version of the simulation including the conditions without 
# missing values and doing five replications per conditions as a suggestion of
# one anonymous reviewer. The purpose of this was to get information about the 
# typical BULK and TAIL ESS of the parameters of the TV-DPCM.

# Select the conditions without missing values.
Cond <- Cond[1:24, ]

# Read output files from small simulation into R
files <- paste0(getwd(), "/Simulation/Sim_TV_DPCM_testcond_", 1:nrow(Cond), ".txt")

data <- lapply(files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame
resultsESS <- ldply(data, data.frame)

rm(data, files)

# Extend condition matrix to match the total number of analyses performed.
Cond_ext <- Cond[resultsESS$cond, ]

# Get ESS column index
col_index <- grep("mbulk|mtail|neff", names(resultsESS))

# Character vectors to define file name and y axis label.
par_name <- c("BULKESS", "TAILESS", "neff")
lab_name <- c("Bulk ESS", "Tail ESS", "n_eff")
sta_name <- rep("Median", 3)

plot_points <- cbind(1:4, 1:4 + 0.1)

for (p in 1:length(par_name)) {
  pdf(file = paste0("Figures/", par_name[p], "_",
                    gsub(" ", "_", sta_name[p]), ".pdf"), width = 15 )
  
  # Define plotting parameters.
  par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 8, 5, 6), xpd = NA)
  
  # Compute mean of the statistic across conditions
  tmp   <- tapply(resultsESS[, col_index[p]], Cond_ext, mean, na.rm = TRUE)
  tmp_l <- tapply(resultsESS[, col_index[p]], Cond_ext, quantile, 
                  probs = 0.25, na.rm = TRUE)
  tmp_u <- tapply(resultsESS[, col_index[p]], Cond_ext, quantile, 
                  probs = 0.75, na.rm = TRUE)
  
  # Define lower and upper limits for the plots.
  lylim <- round(min(tmp_l, na.rm = TRUE) - 0.01, 2)
  uylim <- round(max(tmp_u, na.rm = TRUE) + 0.01, 2)
  if (sta_name[p] == "C.I. Coverage") {lylim <- 0}
  if (sta_name[p] == "Correlation" | sta_name[p] == "C.I. Coverage") {uylim <- 1}
  
  for (i in 1:2) {
    for (l in 1:3) {
      matplot(plot_points, tmp[, i, l,], type = "b", 
              col = color_plot, ylim = c(lylim, uylim), xlab = "", 
              ylab = "", xaxt = "n", yaxt = "n", lwd = 2, pch = 19, lty = 1)
      if (sta_name[p] == "Relative Bias 2.0") {
        abline(h = 1, lty = 2, col = gray(0.9), xpd = FALSE)
      } else {
        abline(h = 0, lty = 2, col = gray(0.9), xpd = FALSE)  
      }
      if (!(sta_name[p] == "C.I. Coverage" & 
            (par_name[p] %in% c("Autoregression", "ProcessVar", "InnoVar")))) {
        segments(x0 = c(plot_points),
                 y0 = c(tmp_l[, i, l, ]),
                 y1 = c(tmp_u[, i, l, ]),
                 lwd = 1, lty = 3, col = rep(color_plot, each = 4))
      }
      if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
      if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
      if (laxis[i, l]) {legend("topright", legend = c("0%", "30%"),
                               col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                               seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
    }
  }
  rm(i, l, tmp, tmp_l, tmp_u, lylim, uylim)
  
  mtext(paste("Average", lab_name[p], sta_name[p]), 2, outer = TRUE, 
        line = 5, cex = 1.8)
  mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
  mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
  mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
  mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
  mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
  mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)
  
  dev.off()
}
rm(p, plot_points, col_index, lab_name, par_name, sta_name)

