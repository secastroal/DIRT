# Summary Simulations
# This file summarizes and plots the results from the final simulations in which
# the AR-PCM and the AR-GRM were tested.

# 0.0 Prepare environment ----
library(plyr)

# 1.0 Read output files ----

# Read output files into R
grm_files <- paste0(getwd(), "/Simulation/Sim_AR_GRM_cond_", 1:72, ".txt")
pcm_files <- paste0(getwd(), "/Simulation/Sim_AR_PCM_cond_", 1:72, ".txt")

data_grm <- lapply(grm_files, function(x) read.table(file = x, header = TRUE))
data_pcm <- lapply(pcm_files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame for each model
results_grm <- ldply(data_grm, data.frame)
results_pcm <- ldply(data_pcm, data.frame)

rm(data_grm, data_pcm, grm_files, pcm_files)

# 2.0 ----

bad_grm <- tapply(results_grm$corrupt, results_grm$cond, sum)
bad_pcm <- tapply(results_pcm$corrupt, results_pcm$cond, sum)

# continue here!
tapply(results_grm$alpha.cor > 0.9, results_grm$cond, sum)

pdf(file = "Figures/Sim_results.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:5, bad_grm[1:5], type = "l", col = "blue", ylim = c(-0.5, 25.5), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:5, bad_pcm[1:5], col = "red", lwd = 2)
lines(1:5, bad_grmf[1:5], col = "darkgreen", lwd = 2)
mtext("I = 3", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:5, bad_grm[6:10], type = "l", col = "blue", ylim = c(-0.5, 25.5), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
lines(1:5, bad_pcm[6:10], col = "red", lwd = 2)
lines(1:5, bad_grmf[6:10], col = "darkgreen", lwd = 2)
mtext("I = 6", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:5, bad_grm[11:15], type = "l", col = "blue", ylim = c(-0.5, 25.5), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:5, bad_pcm[11:15], col = "red", lwd = 2)
lines(1:5, bad_grmf[11:15], col = "darkgreen", lwd = 2)
axis(1, at=1:5, labels=c(60, 120, 200, 350, 500), cex.axis = 2)
mtext("I = 9", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:5, bad_grm[16:20], type = "l", col = "blue", ylim = c(-0.5, 25.5), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
lines(1:5, bad_pcm[16:20], col = "red", lwd = 2)
lines(1:5, bad_grmf[16:20], col = "darkgreen", lwd = 2)
axis(1, at=1:5, labels=c(60, 120, 200, 350, 500), cex.axis = 2)
mtext("I = 12", 3, outer = FALSE, line = -2.3, cex = 2.5)
legend("topright",c("AR-GRM","AR-PCM"), col = c("blue", "red"), lty = 1,
       seg.len = 2, pt.cex = 4, inset = 0)

mtext("Number of Non-Covergent Analyses", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

pdf(file = "theta.pdf", width = 17.2, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(1:nT, theta, type = "l", xlab = "Time Points", ylab = "Theta", ylim = c(-4, 4), 
     las  = 1, cex.axis = 2, cex.lab = 2.5, lwd = 2)
polygon(c(1:nT, rev(1:nT)),
        c(sum.fit$theta[, 4], rev(sum.fit$theta[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.15))
lines(1:nT, sum.fit$theta[, 1], col = "red", lwd = 2)
legend("bottomright", legend = c("True", "Estimated"), col = c("black", "red"), lty = 1, cex = 2)
dev.off()

results_grm_good <- results_grm[results_grm$corrupt == 0, ]
results_grmf_good <- results_grmf[results_grmf$corrupt == 0, ]
results_pcm_good <- results_pcm[results_pcm$corrupt == 0, ]

results_grm_good <- results_grm_good[results_grm_good$cond %in% (1:20)[-(seq(1, 20, by = 5))],]
results_grmf_good <- results_grmf_good[results_grmf_good$cond %in% (1:20)[-(seq(1, 20, by = 5))],]
results_pcm_good <- results_pcm_good[results_pcm_good$cond %in% (1:20)[-(seq(1, 20, by = 5))],]

# Discrimination plots ----
summary_grm <- tapply(results_grm_good$alpha.cor, results_grm_good$cond, mean)

pdf(file = "Figures/alpha_cor.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
mtext("I = 3", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
mtext("I = 6", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
mtext("I = 9", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
mtext("I = 12", 3, outer = FALSE, line = -2.3, cex = 2.5)

mtext("Mean Correlation Discrimination", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

summary_grm <- tapply(results_grm_good$alpha.rmse, results_grm_good$cond, mean)

pdf(file = "Figures/alpha_rmse.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 3", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 6", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 9", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 12", 3, outer = FALSE, line = -2.3, cex = 2.5)

mtext("Mean RMSE Discrimination", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

# Thresholds plots ----

summary_grm <- tapply(results_grm_good$beta.cor, results_grm_good$cond, mean)
summary_grmf <- tapply(results_grmf_good$beta.cor, results_grmf_good$cond, mean)
summary_pcm <- tapply(results_pcm_good$beta.cor, results_pcm_good$cond, mean)

pdf(file = "Figures/beta_cor.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[1:4], col = "red", lwd = 2)
lines(1:4, summary_grmf[1:4], col = "darkgreen", lwd = 2)
mtext("I = 3", 1, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[5:8], col = "red", lwd = 2)
lines(1:4, summary_grmf[5:8], col = "darkgreen", lwd = 2)
mtext("I = 6", 1, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[9:12], col = "red", lwd = 2)
lines(1:4, summary_grmf[9:12], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
mtext("I = 9", 1, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[13:16], col = "red", lwd = 2)
lines(1:4, summary_grmf[13:16], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
mtext("I = 12", 1, outer = FALSE, line = -2.3, cex = 2.5)
legend("bottomright",c("AR-GRM","AR-PCM"), col = c("blue", "red"), lty = 1,
       seg.len = 2, pt.cex = 4, inset = 0)

mtext("Mean Correlation Thresholds", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

summary_grm <- tapply(results_grm_good$beta.rmse, results_grm_good$cond, mean)
summary_grmf <- tapply(results_grmf_good$beta.rmse, results_grmf_good$cond, mean)
summary_pcm <- tapply(results_pcm_good$beta.rmse, results_pcm_good$cond, mean)

pdf(file = "Figures/beta_rmse.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[1:4], col = "red", lwd = 2)
lines(1:4, summary_grmf[1:4], col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 3", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[5:8], col = "red", lwd = 2)
lines(1:4, summary_grmf[5:8], col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 6", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[9:12], col = "red", lwd = 2)
lines(1:4, summary_grmf[9:12], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 9", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 0.5), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[13:16], col = "red", lwd = 2)
lines(1:4, summary_grmf[13:16], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 12", 3, outer = FALSE, line = -2.3, cex = 2.5)
legend("topright",c("AR-GRM","AR-PCM"), col = c("blue", "red"), lty = 1,
       seg.len = 2, pt.cex = 4, inset = 0)

mtext("Mean RMSE Thresholds", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

# Thetas plots ----

summary_grm <- tapply(results_grm_good$theta.cor, results_grm_good$cond, mean)
summary_grmf <- tapply(results_grmf_good$theta.cor, results_grmf_good$cond, mean)
summary_pcm <- tapply(results_pcm_good$theta.cor, results_pcm_good$cond, mean)

pdf(file = "Figures/theta_cor.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[1:4], col = "red", lwd = 2)
lines(1:4, summary_grmf[1:4], col = "darkgreen", lwd = 2)
mtext("I = 3", 1, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[5:8], col = "red", lwd = 2)
lines(1:4, summary_grmf[5:8], col = "darkgreen", lwd = 2)
mtext("I = 6", 1, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[9:12], col = "red", lwd = 2)
lines(1:4, summary_grmf[9:12], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
mtext("I = 9", 1, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[13:16], col = "red", lwd = 2)
lines(1:4, summary_grmf[13:16], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
mtext("I = 12", 1, outer = FALSE, line = -2.3, cex = 2.5)
legend("bottomright",c("AR-GRM","AR-PCM"), col = c("blue", "red"), lty = 1,
       seg.len = 2, pt.cex = 4, inset = 0)

mtext("Mean Correlation Theta", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

summary_grm <- tapply(results_grm_good$theta.rmse, results_grm_good$cond, mean)
summary_grmf <- tapply(results_grmf_good$theta.rmse, results_grmf_good$cond, mean)
summary_pcm <- tapply(results_pcm_good$theta.rmse, results_pcm_good$cond, mean)

pdf(file = "Figures/theta_rmse.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[1:4], col = "red", lwd = 2)
lines(1:4, summary_grmf[1:4], col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 3", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[5:8], col = "red", lwd = 2)
lines(1:4, summary_grmf[5:8], col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 6", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[9:12], col = "red", lwd = 2)
lines(1:4, summary_grmf[9:12], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 9", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 1), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[13:16], col = "red", lwd = 2)
lines(1:4, summary_grmf[13:16], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 12", 3, outer = FALSE, line = -2.3, cex = 2.5)
legend("topright",c("AR-GRM","AR-PCM"), col = c("blue", "red"), lty = 1,
       seg.len = 2, pt.cex = 4, inset = 0)

mtext("Mean RMSE Theta", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()

# Lambda plot ----
summary_grm <- tapply(results_grm_good$lambda.abbias, results_grm_good$cond, mean)
summary_grmf <- tapply(results_grmf_good$lambda.abbias, results_grmf_good$cond, mean)
summary_pcm <- tapply(results_pcm_good$lambda.abbias, results_pcm_good$cond, mean)

pdf(file = "Figures/lambda_abbias.pdf", width = 16.5, height = 8)
par(mfrow = c(2, 2), mar = c(0.2, 0.2, 0.2, 0.2), oma = c(6, 7, 2, 10), xpd = NA)
plot(1:4, summary_grm[1:4], type = "l", col = "blue", ylim = c(-0.05, 0.25), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[1:4], col = "red", lwd = 2)
lines(1:4, summary_grmf[1:4], col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 3", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[5:8], type = "l", col = "blue", ylim = c(-0.05, 0.25), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[5:8], col = "red", lwd = 2)
lines(1:4, summary_grmf[5:8], col = "darkgreen", lwd = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 6", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[9:12], type = "l", col = "blue", ylim = c(-0.05, 0.25), 
     xlab = "", ylab = "", xaxt="n", las = 1, cex.axis = 2, lwd = 2)
lines(1:4, summary_pcm[9:12], col = "red", lwd = 2)
lines(1:4, summary_grmf[9:12], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 9", 3, outer = FALSE, line = -2.3, cex = 2.5)

plot(1:4, summary_grm[13:16], type = "l", col = "blue", ylim = c(-0.05, 0.25), 
     xlab = "", ylab = "", xaxt="n", yaxt = "n", las = 1, lwd = 2)
lines(1:4, summary_pcm[13:16], col = "red", lwd = 2)
lines(1:4, summary_grmf[13:16], col = "darkgreen", lwd = 2)
axis(1, at=1:4, labels=c(120, 200, 350, 500), cex.axis = 2)
abline(h = 0, lty = 2, col = gray(0.5), xpd = FALSE)
mtext("I = 12", 3, outer = FALSE, line = -2.3, cex = 2.5)
legend("topright",c("AR-GRM","AR-PCM"), col = c("blue", "red"), lty = 1,
       seg.len = 2, pt.cex = 4, inset = 0)

mtext("Mean Abbias Autoregressive", 2, outer=TRUE, line=4, cex = 2.5)
mtext("Number of Time Points", 1, outer=TRUE, line=4, cex=2.5)
dev.off()


summary_grm <- tapply(results_grm$run.time, results_grm$cond, mean)
summary_grmf <- tapply(results_grmf$run.time, results_grmf$cond, mean)
summary_pcm <- tapply(results_pcm$run.time, results_pcm$cond, mean)
