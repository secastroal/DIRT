# Additional Plots for the Manuscript

# This script creates plots for the manuscript, which are unrelated with 
# the results of the simulation or the empirical example.

source("R/genTVDPCM.R")

# True Trend -----

y <- sinusoidal(200, maxAbsValue = 1)

pdf(file = "Figures/Trend.pdf", height = 4)
plot(y, type = "l", las = 1, ylim = c(-1.2, 1.2), xaxt = 'n',
     xlab = "", ylab = expression(theta[t]), lwd = 2)
mtext("Time", side = 1, line = 1)
dev.off()
