# Summary Simulations
# This file summarizes and plots the results from the final simulations in which
# the AR-PCM and the AR-GRM were tested.

# 0.0 Prepare environment ----
library(plyr) 

# Next, we create some vector and matrices that are used to create the plots 
# within the for loops.

# String vector to use in the names of the plot files
plot_name <- c("pcm", "grm")

# Recreate matrix of conditions to filter the results

N_timepoints <- c(100, 200, 350, 500) # Number of timepoints
N_items      <- c(3, 6)               # Number of items 
S_lambda     <- c(0, 0.25, 0.5)       # Size of the autoregressive effect
M_prop       <- c(0, 0.15, 0.25)      # Proportion of missing values

Cond        <- expand.grid(N_timepoints, N_items, S_lambda, M_prop)
names(Cond) <- c("nT", "I", "lambda", "NAprop")

# Colors for plots
# If bnw = TRUE, plots are printed in grayscale.
bnw = TRUE #black and white
if (bnw) {
        color_plot <- gray((0:2)/3)
}else{
        color_plot <- rainbow(3)
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
grm_files <- paste0(getwd(), "/Simulation/Sim_AR_GRM_cond_", 1:72, ".txt")
pcm_files <- paste0(getwd(), "/Simulation/Sim_AR_PCM_cond_", 1:72, ".txt")

data_grm <- lapply(grm_files, function(x) read.table(file = x, header = TRUE))
data_pcm <- lapply(pcm_files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame for each model
results_grm <- ldply(data_grm, data.frame)
results_pcm <- ldply(data_pcm, data.frame)

# Store results in a list to loop through it and use lapply()
results <- list(pcm = results_pcm, grm = results_grm)

rm(data_grm, data_pcm, grm_files, pcm_files, results_pcm, results_grm)

# 2.0 Analyses that failed ----
# To check for convergence issues in the analyses, we looked at the warning messages from stan.
# A summary of the warnings in Stan can be found here: https://mc-stan.org/misc/warnings.html
# We considered that an analysis did not converge or the estimation is not reliable enough if
# there were any rhat larger than 1.05, any divergent transition after warmup, or any BMFI was 
# too low. These analyses are identified with a 1 in the variable "corrupt" (e.g., results_pcm$corrupt).
# There are other warning messages that do not necessarily imply that the analyses failed or is 
# unreliable but they indicate that the estimation was inefficient. These warning messages 
# indicate that the maximum treedepth was exceeded or that the bulk or tail effective sample
# sizes were too low. These warning messages are not as serious as having divergent transitions
# and they are usually solved by increasing the number of iteration or the maximum treedepth 
# of the algorithm. Analyses that showed any of these warnings are identified with a 1 in the 
# variable "efficiency".

n_failed <- lapply(results, function(x) tapply(x$corrupt, x$cond, sum))

for (m in 1:2) {
  # Save plot to pdf
  pdf(file = paste0("Figures/Divergent_", plot_name[m], ".pdf"), width = 15)
  # Define plotting parameters.
  par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)

  for (i in 1:2) {
    for(l in 1:3) {
      plot(1:4, 
           n_failed[[m]][which(Cond$I == N_items[i] & Cond$lambda == S_lambda[l] &
                                 Cond$NAprop == M_prop[1])] / 5, 
           type = "l", col = color_plot[1], ylim = c(-0.5, 100.5), xlab = "", 
           ylab = "", xaxt = "n", yaxt = "n", lwd = 2)
      lines(1:4, 
            n_failed[[m]][which(Cond$I == N_items[i] & Cond$lambda == S_lambda[l] & 
                                  Cond$NAprop == M_prop[2])] / 5,
            col = color_plot[2], lwd = 2)
      lines(1:4, 
            n_failed[[m]][which(Cond$I == N_items[i] & Cond$lambda == S_lambda[l] & 
                                  Cond$NAprop == M_prop[3])] / 5,
            col = color_plot[3], lwd = 2)
      if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
      if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
      if (laxis[i, l]) {legend("topright", legend = c("0%", "15%", "25%"),
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
}
rm(i, l, m)

# 3.1 Analyses that converged ----

# To plot the average bias, average rmse, and other statistics across conditions,
# we filter the data to only include the analyses that converged.
results_conv <- lapply(results, function(x) x[x$corrupt == 0, ])
# From the analyses that converged, there was still one analysis that had warning
# messages. The 9th replication of the AR-GRM in condition 64 (nT = 500, I = 6,
# lambda = 0.25, and NAprop = 0.25) had 9 transitions that exceeded the maximum 
# treedepth after warm-up. Given that this problem is mainly an efficiency concern,
# this analysis is still included in the following summaries.
lapply(results_conv, function (x) tapply(x$efficiency, x$cond, sum))
results_conv$grm[results_conv$grm$cond == 64 & results_conv$grm$efficiency == 1, 1:10]

# 3.2 True vs Estimates Correlations ----

# In the following plots, the average correlations between the true and the estimated
# values for the alpha, beta, and theta parameters are displayed. The dots represent 
# the average correlation across replications per conditions and the lines represent
# the spread of the correlation given the 2.5 and 97.5 percentiles.

# Beta parameters
for (m in 1:2) {
  # Save plot to pdf
  pdf(file = paste0("Figures/Beta_Cor_", plot_name[m], ".pdf"), width = 15)
  # Define plotting parameters.
  par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)
  
  for (i in 1:2) {
    for (l in 1:3) {
      tmp <- results_conv[[m]][which(results_conv[[m]]$cond %in% 
                                       which(Cond$I == N_items[i] & Cond$lambda == S_lambda[l])),]
      
      tmp_up   <- tapply(tmp$beta.cor, tmp$cond, quantile, probs = 0.975)
      tmp_down <- tapply(tmp$beta.cor, tmp$cond, quantile, probs = 0.025)
      
      plot(rep(1:4, 3) + rep((0:2) / 10, each = 4),
           tapply(tmp$beta.cor, tmp$cond, mean), 
           ylim = c(0.75, 1), type = "p", pch = 19, 
           col = rep(color_plot, each = 4), xaxt = "n", yaxt = "n",
           xlab = "", ylab = "")
      segments(x0 = rep(1:4, 3) + rep((0:2) / 10, each = 4),
               y0 = tmp_down,
               y1 = tmp_up,
               lwd = 2, col = rep(color_plot, each = 4))
      if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
      if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
      if (laxis[i, l]) {legend("topright", legend = c("0%", "15%", "25%"),
                               col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                               seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
    }
  }
  
  mtext("Average Thresholds Correlation", 2, outer = TRUE, line = 4, cex = 1.8)
  mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
  mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
  mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
  mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
  mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
  mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)
  
  dev.off()
}
rm(tmp, tmp_down, tmp_up, i, m , l)

# Theta parameters
lylim <- c(0.6, 0.4)
  
for (m in 1:2) {
  # Save plot to pdf
  pdf(file = paste0("Figures/Theta_Cor_", plot_name[m], ".pdf"), width = 15)
  # Define plotting parameters.
  par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)
  
  for (i in 1:2) {
    for (l in 1:3) {
      tmp <- results_conv[[m]][which(results_conv[[m]]$cond %in% 
                                       which(Cond$I == N_items[i] & Cond$lambda == S_lambda[l])),]
      
      tmp_up   <- tapply(tmp$theta.cor, tmp$cond, quantile, probs = 0.975)
      tmp_down <- tapply(tmp$theta.cor, tmp$cond, quantile, probs = 0.025)
      
      plot(rep(1:4, 3) + rep((0:2) / 10, each = 4),
           tapply(tmp$theta.cor, tmp$cond, mean), 
           ylim = c(lylim[m], 1), type = "p", pch = 19, 
           col = rep(color_plot, each = 4), xaxt = "n", yaxt = "n",
           xlab = "", ylab = "")
      segments(x0 = rep(1:4, 3) + rep((0:2) / 10, each = 4),
               y0 = tmp_down,
               y1 = tmp_up,
               lwd = 2, col = rep(color_plot, each = 4))
      if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
      if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
      if (laxis[i, l]) {legend("topright", legend = c("0%", "15%", "25%"),
                               col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                               seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
    }
  }
  
  mtext("Average State Disposition Correlation", 2, outer = TRUE, line = 4, cex = 1.8)
  mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
  mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
  mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
  mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
  mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
  mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)
  
  dev.off()
}
rm(tmp, tmp_down, tmp_up, i, m , l, lylim)

# Alpha parameters
# Save plot to pdf
pdf(file = paste0("Figures/Alpha_Cor_grm.pdf"), width = 15)
# Define plotting parameters.
par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)

for (i in 1:2) {
  for (l in 1:3) {
    tmp <- results_conv$grm[which(results_conv$grm$cond %in% 
                                    which(Cond$I == N_items[i] & Cond$lambda == S_lambda[l])),]
    
    tmp_up   <- tapply(tmp$alpha.cor, tmp$cond, quantile, probs = 0.975)
    tmp_down <- tapply(tmp$alpha.cor, tmp$cond, quantile, probs = 0.025)
    
    plot(rep(1:4, 3) + rep((0:2) / 10, each = 4),
         tapply(tmp$alpha.cor, tmp$cond, mean), 
         ylim = c(-1, 1), type = "p", pch = 19, 
         col = rep(color_plot, each = 4), xaxt = "n", yaxt = "n",
         xlab = "", ylab = "")
    segments(x0 = rep(1:4, 3) + rep((0:2) / 10, each = 4),
             y0 = tmp_down,
             y1 = tmp_up,
             lwd = 2, col = rep(color_plot, each = 4))
    if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
    if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
    if (laxis[i, l]) {legend("topright", legend = c("0%", "15%", "25%"),
                             col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                             seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
  }
}

mtext("Average Discrimination Correlation", 2, outer = TRUE, line = 4, cex = 1.8)
mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)

dev.off()

rm(i, l, tmp, tmp_up, tmp_down)

# 3.3 True vs. Estimates Statistics ----

# To assess how close or distant were the estimates from the true parameters,
# we computed the average bias, absolute bias, and root mean squared error for 
# each set of parameters.  

# Plots: Beta, theta, and lambda.
# Column indexes of the statistics of interest
col_index <- sapply(results_conv, function(x) {
  c(grep("beta",   names(x))[-c(1, 5)],
    grep("theta",  names(x))[-c(1, 5)],
    grep("lambda", names(x))[-3])
  })

# Character vectors to use in the for loop to define the file name and y axis label.
par_name <- rep(c("Beta", "Theta", "Lambda"), times = c(3, 3, 2))
lab_name <- rep(c("Thresholds", "State Disposition", "Autoregression"), times = c(3, 3, 2))
sta_name <- rep(c("Bias", "Abbias", "RMSE"), length.out = length(par_name))

for (m in 1:2) {
  for (v in 1:length(par_name)) {
    # Save plot to pdf
    pdf(file = paste0("Figures/", par_name[v], "_", sta_name[v], "_", plot_name[m], ".pdf"), width = 15)
    # Define plotting parameters.
    par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)
    
    lylim <- min(tapply(results_conv[[m]][, col_index[v, m]], 
                        results_conv[[m]]$cond, quantile, probs = 0.025))
    lylim <- round(lylim - 0.05, 2)
    uylim <- max(tapply(results_conv[[m]][, col_index[v, m]], 
                        results_conv[[m]]$cond, quantile, probs = 0.975))
    uylim <- round(uylim + 0.05, 2)
    
    for (i in 1:2) {
      for (l in 1:3) {
        tmp <- results_conv[[m]][which(results_conv[[m]]$cond %in% 
                                         which(Cond$I == N_items[i] & 
                                                 Cond$lambda == S_lambda[l])),]
        
        tmp_up   <- tapply(tmp[, col_index[v, m]], tmp$cond, quantile, probs = 0.975)
        tmp_down <- tapply(tmp[, col_index[v, m]], tmp$cond, quantile, probs = 0.025)
        
        plot(rep(1:4, 3) + rep((0:2) / 10, each = 4),
             tapply(tmp[, col_index[v, m]], tmp$cond, mean), 
             ylim = c(lylim, uylim), type = "p", pch = 19, 
             col = rep(color_plot, each = 4), xaxt = "n", yaxt = "n",
             xlab = "", ylab = "")
        abline(h = 0, lty = 2, col = gray(0.9), xpd = FALSE)
        segments(x0 = rep(1:4, 3) + rep((0:2) / 10, each = 4),
                 y0 = tmp_down,
                 y1 = tmp_up,
                 lwd = 2, col = rep(color_plot, each = 4))
        if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
        if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
        if (laxis[i, l]) {legend("topright", legend = c("0%", "15%", "25%"),
                                 col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                                 seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
      }
    }
    
    mtext(paste("Average", lab_name[v], sta_name[v]), 2, outer = TRUE, line = 4, cex = 1.8)
    mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
    mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
    mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
    mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
    mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
    mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)
    
    dev.off()
  }
  }
rm(tmp, tmp_down, tmp_up, i, m, l, v, lylim, uylim,
   par_name, lab_name, sta_name, col_index)

# Plots: Alpha.

col_index <- grep("alpha", names(results_conv$grm))[-c(1, 5)]
sta_name  <- c("Bias", "Abbias", "RMSE")

for (v in 1:length(col_index)) {
  # Save plot to pdf
  pdf(file = paste0("Figures/Alpha_", sta_name[v], "_grm.pdf"), width = 15)
  # Define plotting parameters.
  par(mfrow = c(2, 3), mar = c(0.2, 0.2, 0.2, 1.2), oma = c(6, 7, 5, 10), xpd = NA)
  
  lylim <- min(tapply(results_conv[[2]][, col_index[v]], 
                      results_conv[[2]]$cond, quantile, probs = 0.025))
  lylim <- round(lylim - 0.05, 2)
  uylim <- max(tapply(results_conv[[2]][, col_index[v]], 
                      results_conv[[2]]$cond, quantile, probs = 0.975))
  uylim <- round(uylim + 0.05, 2)
  
  for (i in 1:2) {
    for (l in 1:3) {
      tmp <- results_conv$grm[which(results_conv$grm$cond %in% 
                                       which(Cond$I == N_items[i] & 
                                               Cond$lambda == S_lambda[l])),]
      
      tmp_up   <- tapply(tmp[, col_index[v]], tmp$cond, quantile, probs = 0.975)
      tmp_down <- tapply(tmp[, col_index[v]], tmp$cond, quantile, probs = 0.025)
      
      plot(rep(1:4, 3) + rep((0:2) / 10, each = 4),
           tapply(tmp[, col_index[v]], tmp$cond, mean), 
           ylim = c(lylim, uylim), type = "p", pch = 19, 
           col = rep(color_plot, each = 4), xaxt = "n", yaxt = "n",
           xlab = "", ylab = "")
      abline(h = 0, lty = 2, col = gray(0.9), xpd = FALSE)
      segments(x0 = rep(1:4, 3) + rep((0:2) / 10, each = 4),
               y0 = tmp_down,
               y1 = tmp_up,
               lwd = 2, col = rep(color_plot, each = 4))
      if (yaxis[i, l]) {axis(2, labels = TRUE, cex.axis = 2, las = 1)}
      if (xaxis[i, l]) {axis(1, at = 1:4, labels = N_timepoints, cex.axis = 2)}
      if (laxis[i, l]) {legend("topright", legend = c("0%", "15%", "25%"),
                               col = color_plot, lty = 1, lwd = 2, cex = 1.5,
                               seg.len = 3, title = "% of NA", inset = c(-1/3, 0))}
    }
  }
  
  mtext(paste("Average", "Discrimination", sta_name[v]), 2, outer = TRUE, line = 4, cex = 1.8)
  mtext("Number of Time Points", 1, outer = TRUE, line = 4, cex = 1.8)
  mtext("I = 3", 4, outer = TRUE, line = 1, cex = 1.5, at = 3 / 4, las = 1)
  mtext("I = 6", 4, outer = TRUE, line = 1, cex = 1.5, at = 1 / 4, las = 1)
  mtext("AR = 0",    3, outer = TRUE, line = 1, cex = 1.5, at = 1 / 6)
  mtext("AR = 0.25", 3, outer = TRUE, line = 1, cex = 1.5, at = 3 / 6)
  mtext("AR = 0.5",  3, outer = TRUE, line = 1, cex = 1.5, at = 5 / 6)
  
  dev.off()
}
rm(tmp, tmp_down, tmp_up, i, l, v, lylim, uylim,
   sta_name, col_index)

# continue here!
summary_grm <- tapply(results_grm$run.time, results_grm$cond, mean)
summary_grmf <- tapply(results_grmf$run.time, results_grmf$cond, mean)
summary_pcm <- tapply(results_pcm$run.time, results_pcm$cond, mean)
