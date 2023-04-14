# Summary PPMC Simulations

# This file plots the results from the simulation study of the PPMC methods
# for the TV-DPCM model.

# 0.0 Prepare environment
# 1.0 Read output files
# 2.0 Check convergence per condition
# 3.0 Plot PPPs distributions
# 4.0 Tables of measures' power 

# 0.0 Prepare environment ----
library(plyr)
library(xtable)

## 0.1 Recreate simulation conditions ----

# Fix conditions

nT     <- 300   # Number of time points
I      <- 6     # Number of items
K      <- 5     # Number of categories per item
M      <- K - 1 # Number of thresholds per item
lambda <- 0.5   # Size of the autoregressive effect
in_var <- 1     # Variance of the innovations

# Manipulated conditions
# Generating model.
#   'TV-DPCM': No violation to assumptions of TV-DPCM
#   'BiDim': Unidimensionality is violated.
#   'GPCM': Discrimination is not equal for all items.
#   'Drift': Longitudinal measurement invariance does not hold. Parameter Drift.
#   'Meaning': The meaning of the items changes. 
#              This is also measurement non-invariance
gen.model <- c("TV-DPCM", "BiDim", "BiDim", "GPCM", "GPCM", 
               "DRIFT", "DRIFT", "Meaning")
# Parameters of the generating model. Correlation between latent variables for 
# 'BiDim' and proportion of items with discrimination different from 1 or with
# item parameter drift for 'GPCM' and 'DRIFT' respectively.
par.model <- c(NA, 0.3, 0.6, 1/3, 2/3, 1/3, 2/3, NA)

# Matrix of conditions to loop through
Cond <- data.frame(gen.model, par.model)
names(Cond) <- c("genModel", "parModel")
rm(gen.model, par.model)

# 1.0 Read output files ----

# Read output files into R
files <- paste0(getwd(), "/Simulation/Sim_PPMC_cond_", 1:nrow(Cond), ".txt")

data <- lapply(files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame
results <- ldply(data, data.frame)

rm(data, files)


results3_1 <- read.table("Simulation/Results_3_items.dat", header = TRUE)
results3_2 <- read.table("Simulation/Results_part2_3_items.dat", header = TRUE)
results6_1 <- read.table("Simulation/Results_6_items.dat", header = TRUE)
results6_2 <- read.table("Simulation/Results_part2_6_items.dat", header = TRUE)
results12_1 <- read.table("Simulation/Results_12_items.dat", header = TRUE)
results12_2 <- read.table("Simulation/Results_part2_12_items.dat", header = TRUE)

names(results3_1) <- names(results3_2)
results3 <- rbind(results3_1, results3_2)
results3 <- results3[order(results3$model), ]
names(results6_1) <- names(results6_2)
results6 <- rbind(results6_1, results6_2)
results6 <- results6[order(results6$model), ]
names(results12_1) <- names(results12_2)
results12 <- rbind(results12_1, results12_2)
results12 <- results12[order(results12$model), ]

rm(results3_1, results3_2, 
   results6_1, results6_2, 
   results12_1, results12_2)

sumresults3 <- apply(results3[, -(1:8)], 2, function(x) {
  tapply(x, results3$model, function(y) {
    sum(y <= 0.05 | y >= 0.95)
  })
})

sumresults3 <- t(sumresults3)

sumresults6 <- apply(results6[, -(1:8)], 2, function(x) {
  tapply(x, results6$model, function(y) {
    sum(y <= 0.05 | y >= 0.95)
  })
})

sumresults6 <- t(sumresults6)

sumresults12 <- apply(results12[, -(1:8)], 2, function(x) {
  tapply(x, results12$model, function(y) {
    sum(y <= 0.05 | y >= 0.95)
  })
})

sumresults12 <- t(sumresults12)

nmis3 <- apply(results3[, -(1:8)], 1, function(x) sum(x <= 0.05 | x >= 0.95))
nmis6 <- apply(results6[, -(1:8)], 1, function(x) sum(x <= 0.05 | x >= 0.95))
nmis12 <- apply(results12[, -(1:8)], 1, function(x) sum(x <= 0.05 | x >= 0.95))
pmis3 <- apply(results3[, -(1:8)], 1, function(x) round(sum(x <= 0.05 | x >= 0.95)/39, 2))
pmis6 <- apply(results6[, -(1:8)], 1, function(x) round(sum(x <= 0.05 | x >= 0.95)/117, 2))
pmis12 <- apply(results12[, -(1:8)], 1, function(x) round(sum(x <= 0.05 | x >= 0.95)/408, 2))

countresults3 <- data.frame(results3[, 1:2], nmis3, pmis3)
countresults6 <- data.frame(results6[, 1:2], nmis6, pmis6)
countresults12 <- data.frame(results12[, 1:2], nmis12, pmis12)

xlsx::write.xlsx(countresults3, file = "Simulation/Results3Items.xlsx")
xlsx::write.xlsx(countresults6, file = "Simulation/Results6Items.xlsx")


hist(unlist(results3[results3$model == "TV-DPCM", grep("itcor3.I", names(results3))]))
hist(unlist(results6[results6$model == "TV-DPCM", grep("itcor.I", names(results6))]))
hist(unlist(results12[results12$model == "TV-DPCM", grep("itcor.I", names(results12))]))

modeltype <- unique(results3$model)
discmeasures <- c("ACF.1", "ACF.2", "ACF.3", "PACF", "LPACF", "MSSD",
                  "itcor.I", "itcor2.", "itcor3.", "q1.I", "q1alt.",
                  "lpacf.", "q3.", "OR.", "ordiff.", "resid.", "rediff.")

out.base <- data.frame(matrix(NA, length(discmeasures), length(modeltype)))
names(out.base) <- modeltype
row.names(out.base) <- discmeasures

out3 <- out6 <- out12 <- out.base

for (i in 1:length(modeltype)) {
  for (j in 1:length(discmeasures)) {
    tmp.data <- results3
    tmp <- unlist(tmp.data[tmp.data$model == modeltype[i], 
                           grep(discmeasures[j], names(tmp.data))])
    out3[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
}
rm(tmp, tmp.data, i, j)

for (i in 1:length(modeltype)) {
  for (j in 1:length(discmeasures)) {
    tmp.data <- results6
    tmp <- unlist(tmp.data[tmp.data$model == modeltype[i], 
                           grep(discmeasures[j], names(tmp.data))])
    out6[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
}
rm(tmp, tmp.data, i, j)

for (i in 1:length(modeltype)) {
  for (j in 1:length(discmeasures)) {
    tmp.data <- results12
    tmp <- unlist(tmp.data[tmp.data$model == modeltype[i], 
                           grep(discmeasures[j], names(tmp.data))])
    out12[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
}
rm(tmp, tmp.data, i, j)

# Create tables for manuscript



Results <- list(results3, results6, results12)

# 2.0 Check convergence per condition ----
# 3.0 Plot PPPs distributions ----

# Plot the distribution of the PPP-values

discme_labels <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "RACF", "LRACF", "MSSD",
                   "Item-total Correlation", "Item-total Correlation (v2)", 
                   "Item-total Correlation (v3)", expression(paste("Yen's ", Q[1])), 
                   expression(paste("Yen's ", Q[1], " alt.")),
                   "Item LRACF", expression(paste("Yen's ", Q[3])), 
                   "OR", "OR difference", 
                   "RESID", "RESID  difference") 

pdf("Figures/ppp_dist_item_measures.pdf", height = 5)
par(mfrow = c(2, 3), mar = c(2, 4, 4, 1) + 0.1, oma = c(0, 1, 0, 1) + 0.1)
for (i in 7:12) {
  hist(unlist(results12[results12$model == "TV-DPCM", grep(discmeasures[i], names(results12))]),
       las = 1, xlim = c(0, 1), main = discme_labels[i], xlab = "", 
       freq = FALSE, breaks = seq(0, 1, by = 0.05), 
       col = c("black", rep("lightgray", 18), "black"), 
       border = c("black", rep("lightgray", 18), "black"))
}
rm(i)
dev.off()

pdf("Figures/ppp_dist_pair_measures.pdf", height = 5)
par(mfrow = c(2, 3), mar = c(2, 4, 4, 1) + 0.1, oma = c(0, 1, 0, 1) + 0.1)
for (i in 13:17) {
  hist(unlist(results12[results12$model == "TV-DPCM", grep(discmeasures[i], names(results12))]),
       las = 1, xlim = c(0, 1), main = discme_labels[i], xlab = "", 
       freq = FALSE, breaks = seq(0, 1, by = 0.05), 
       col = c("black", rep("lightgray", 18), "black"), 
       border = c("black", rep("lightgray", 18), "black"))
}
rm(i)
dev.off()

# 4.0 Tables of measures' power  ----

# TV-DPCM
outTVDPCM <- matrix(NA, length(discmeasures), 3)

for(i in 1:3) {
  for (j in 1:length(discmeasures)) {
    tmp.data <- Results[[i]][Results[[i]]$model == "TV-DPCM", ]
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outTVDPCM[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
}

outTVDPCM <- ifelse(outTVDPCM >= 0.7,
                    paste0("\\textbf{", round(outTVDPCM, 2), "}"),
                    round(outTVDPCM, 2))

outTVDPCM <- data.frame(outTVDPCM)
names(outTVDPCM) <- paste0("$I=", c(3, 6, 12), "$")
row.names(outTVDPCM) <- paste0("\\hspace{0.25in}",
                               c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", 
                                 "LPACF", "MSSD", "Item-total correlation", 
                                 "Item-total correlation (v2)", "Item-total correlation (v3)", 
                                 "Yen's $Q_1$", "Yen's $Q_1$ alt.", "Item LPACF", 
                                 "Yen's $Q_3$", "OR", "OR difference",
                                 "RESID", "RESID  difference"))

addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outTVDPCM))
addtorow$command <- c("& \\multicolumn{3}{c}{Number of Items} \\\\\n",
                      " & $I = 3$ & $I = 6$ & $I = 12$\\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} Test one.}\n\\end{tablenotes}\n")

print(xtable(outTVDPCM, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-DPCM as the Generating Model",
             label = "tab:tvdpcm", align = c("l", "r", "r", "r")),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outTVDPCM.tex")

# Bidimensional
nItems <- c(3, 6, 12)
pairsindex <- list(c(1, 2, 2),
                   rep(c(1, 2, 1, 2, 3), times = c(2, 3, 1, 6, 3)),
                   rep(c(rep(1:2, 5), 3), times = c(5, 6, 4, 6, 3, 6, 2, 6, 1, 12, 15)))

outBiDim <- matrix(NA, 6 + 6 * 2 + 5 * 3, 3)

for(i in 1:3) {
  tmp.data <- Results[[i]][Results[[i]]$model == "BiDim", ]
  for (j in 1:6) {
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outBiDim[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
  rm(tmp)
  for (j in 7:12) {
    tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[1:ceiling(nItems[i]/2)]])
    tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[(ceiling(nItems[i]/2) + 1):nItems[i]]])
    outBiDim[2 * j - 7, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outBiDim[2 * j - 6, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
  }
  rm(tmp1, tmp2)
  for (j in 13:17) {
    tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 1]])
    tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 2]])
    tmp3 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 3]])
    outBiDim[3 * j - 20, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outBiDim[3 * j - 19, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
    outBiDim[3 * j - 18, i] <- round(sum(tmp3 <= 0.05 | tmp3 >= 0.95)/length(tmp3), 2)
  }
  rm(tmp1, tmp2, tmp3)
}

outBiDim <- ifelse(outBiDim >= 0.7,
                    paste0("\\textbf{", round(outBiDim, 2), "}"),
                    round(outBiDim, 2))

outBiDim <- data.frame(outBiDim)
names(outBiDim) <- paste0("I=", c(3, 6, 12))
row.names(outBiDim) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", "LPACF", "MSSD",
                        paste("Item-total correlation", c("dim1", "dim2")),
                        paste("Item-total correlation (v2)", c("dim1", "dim2")),
                        paste("Item-total correlation (v3)", c("dim1", "dim2")),
                        paste("Yen's $Q_1$", c("dim1", "dim2")),
                        paste("Yen's $Q_1$ alt.", c("dim1", "dim2")),
                        paste("Item LPACF", c("dim1", "dim2")),
                        paste("Yen's $Q_3$", c("(dim1, dim1)", "(dim1, dim2)", "(dim2, dim2)")),
                        paste("OR", c("(dim1, dim1)", "(dim1, dim2)", "(dim2, dim2)")),
                        paste("OR difference", c("(dim1, dim1)", "(dim1, dim2)", "(dim2, dim2)")),
                        paste("RESID", c("(dim1, dim1)", "(dim1, dim2)", "(dim2, dim2)")),
                        paste("RESID  difference",c("(dim1, dim1)", "(dim1, dim2)", "(dim2, dim2)")))
row.names(outBiDim) <- paste0("\\hspace{0.25in}", row.names(outBiDim))

outBiDim <- outBiDim[-(1:6), ]

addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outBiDim))
addtorow$command <- c("& \\multicolumn{3}{c}{Number of Items} \\\\\n",
                      " & $I = 3$ & $I = 6$ & $I = 12$\\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} \\textit{dim1} denotes the items of dimension 1.
                      \\textit{dim2} denotes the items of dimension 2.}\n\\end{tablenotes}\n")

print(xtable(outBiDim, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-MDPCM as the Generating Model",
             label = "tab:tvmdpcm", align = c("l", "r", "r", "r")),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outTVMDPCM.tex")

# Default Responses

outDefault <- matrix(NA, length(discmeasures), 3)

for(i in 1:3) {
  for (j in 1:length(discmeasures)) {
    tmp.data <- Results[[i]][Results[[i]]$model == "Default", ]
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outDefault[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
}

outDefault <- ifelse(outDefault >= 0.7,
                   paste0("\\textbf{", round(outDefault, 2), "}"),
                   round(outDefault, 2))

outDefault <- data.frame(outDefault)
names(outDefault) <- paste0("I=", c(3, 6, 12))
row.names(outDefault) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", "LPACF", "MSSD",
                         "Item-total correlation", "Item-total correlation (v2)", 
                         "Item-total correlation (v3)", "Yen's $Q_1$", "Yen's $Q_1$ alt.",
                         "Item LPACF", "Yen's $Q_3$", "OR", "OR difference", 
                         "RESID", "RESID  difference")
row.names(outDefault) <- paste0("\\hspace{0.25in}", row.names(outDefault))

outDefault <- outDefault[-(1:6), ]

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{3}{c}{Number of Items} \\\\\n",
                      " & $I = 3$ & $I = 6$ & $I = 12$\\\\\n")

print(xtable(outDefault, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-DPCM-Meaning as the Generating Model",
             label = "tab:default", align = c("l", "r", "r", "r")),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outDefault.tex")

# Random Responses

outRandom <- matrix(NA, length(discmeasures), 3)

for(i in 1:3) {
  for (j in 1:length(discmeasures)) {
    tmp.data <- Results[[i]][Results[[i]]$model == "Random", ]
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outRandom[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
}

outRandom <- data.frame(outRandom)
names(outRandom) <- paste0("I=", c(3, 6, 12))
row.names(outRandom) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", "LPACF", "MSSD",
                           "Item-total correlation", "Item-total correlation (v2)", 
                           "Item-total correlation (v3)", "Yen's $Q_1$", "Yen's $Q_1$ alt.",
                           "Item LPACF", "Yen's $Q_3$", "OR", "OR difference", 
                           "RESID", "RESID  difference")

print(xtable(outRandom, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-DPCM-Random",
             label = "tab:random", align = c("l", "c", "c", "c")),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", sanitize.text.function = function(x){x},
      file = "Tables/outRandom.tex")

# GPCM different discrimination

pairsindex <- list(c(2, 2, 3),
                   rep(c(rep(1:2, 3), 3), times = c(3, 2, 2, 2, 1, 4, 1)),
                   rep(c(rep(1:2, 9), 3), times = c(9, 2, 8, 2, 7, 2, 6, 2, 5, 2,
                                                    4, 2, 3, 2, 2, 2, 1, 4, 1)))

outGPCM <- matrix(NA, 6 + 6 * 2 + 5 * 3, 3)

for(i in 1:3) {
  tmp.data <- Results[[i]][Results[[i]]$model == "GPCM", ]
  for (j in 1:6) {
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outGPCM[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
  rm(tmp)
  for (j in 7:12) {
    tmp1 <- unlist(tmp.data[, head(grep(discmeasures[j], names(tmp.data)), -2)])
    tmp2 <- unlist(tmp.data[, tail(grep(discmeasures[j], names(tmp.data)), 2)])
    outGPCM[2 * j - 7, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outGPCM[2 * j - 6, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
  }
  rm(tmp1, tmp2)
  for (j in 13:17) {
    tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 1]])
    tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 2]])
    tmp3 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 3]])
    outGPCM[3 * j - 20, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outGPCM[3 * j - 19, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
    outGPCM[3 * j - 18, i] <- round(sum(tmp3 <= 0.05 | tmp3 >= 0.95)/length(tmp3), 2)
  }
  rm(tmp1, tmp2, tmp3)
}

outGPCM <- ifelse(outGPCM >= 0.7,
                   paste0("\\textbf{", round(outGPCM, 2), "}"),
                   round(outGPCM, 2))

outGPCM <- data.frame(outGPCM)
names(outGPCM) <- paste0("I=", c(3, 6, 12))
row.names(outGPCM) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", "LPACF", "MSSD",
                        paste("Item-total correlation", c("disc1", "disc2")),
                        paste("Item-total correlation (v2)", c("disc1", "disc2")),
                        paste("Item-total correlation (v3)", c("disc1", "disc2")),
                        paste("Yen's $Q_1$", c("disc1", "disc2")),
                        paste("Yen's $Q_1$ alt.", c("disc1", "disc2")),
                        paste("Item LPACF", c("disc1", "disc2")),
                        paste("Yen's $Q_3$", c("(disc1, disc1)", "(disc1, disc2)", "(disc2, disc2)")),
                        paste("OR", c("(disc1, disc1)", "(disc1, disc2)", "(disc2, disc2)")),
                        paste("OR difference", c("(disc1, disc1)", "(disc1, disc2)", "(disc2, disc2)")),
                        paste("RESID", c("(disc1, disc1)", "(disc1, disc2)", "(disc2, disc2)")),
                        paste("RESID  difference",c("(disc1, disc1)", "(disc1, disc2)", "(disc2, disc2)")))

row.names(outGPCM) <- paste0("\\hspace{0.25in}", row.names(outGPCM))

outGPCM <- outGPCM[-(1:6), ]

addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outGPCM))
addtorow$command <- c("& \\multicolumn{3}{c}{Number of Items} \\\\\n",
                      " & $I = 3$ & $I = 6$ & $I = 12$\\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} \\textit{disc} denotes the items with discrimination parameter equal to 1.
                      \\textit{disc2} denotes the items with discrimination parameters different from 1.}\\end{tablenotes}\n")

print(xtable(outGPCM, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-DGPCM as the Generating Model",
             label = "tab:gpcm", align = c("l", "r", "r", "r")),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outGPCM.tex")

# Half Drift

pairsindex <- list(c(2, 2, 3),
                   rep(c(rep(1:2, 3), 3), times = c(3, 2, 2, 2, 1, 4, 1)),
                   rep(c(rep(1:2, 9), 3), times = c(9, 2, 8, 2, 7, 2, 6, 2, 5, 2,
                                                    4, 2, 3, 2, 2, 2, 1, 4, 1)))

outHDrift <- matrix(NA, 6 + 6 * 2 + 5 * 3, 3)

for(i in 1:3) {
  tmp.data <- Results[[i]][Results[[i]]$model == "HalfDrift", ]
  for (j in 1:6) {
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outHDrift[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
  rm(tmp)
  for (j in 7:12) {
    tmp1 <- unlist(tmp.data[, head(grep(discmeasures[j], names(tmp.data)), -2)])
    tmp2 <- unlist(tmp.data[, tail(grep(discmeasures[j], names(tmp.data)), 2)])
    outHDrift[2 * j - 7, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outHDrift[2 * j - 6, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
  }
  rm(tmp1, tmp2)
  for (j in 13:17) {
    tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 1]])
    tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 2]])
    tmp3 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 3]])
    outHDrift[3 * j - 20, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outHDrift[3 * j - 19, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
    outHDrift[3 * j - 18, i] <- round(sum(tmp3 <= 0.05 | tmp3 >= 0.95)/length(tmp3), 2)
  }
  rm(tmp1, tmp2, tmp3)
}

outHDrift <- ifelse(outHDrift >= 0.7,
                   paste0("\\textbf{", round(outHDrift, 2), "}"),
                   round(outHDrift, 2))

outHDrift <- data.frame(outHDrift)
names(outHDrift) <- paste0("I=", c(3, 6, 12))
row.names(outHDrift) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", "LPACF", "MSSD",
                        paste("Item-total correlation", c("drift1", "drift2")),
                        paste("Item-total correlation (v2)", c("drift1", "drift2")),
                        paste("Item-total correlation (v3)", c("drift1", "drift2")),
                        paste("Yen's $Q_1$", c("drift1", "drift2")),
                        paste("Yen's $Q_1$ alt.", c("drift1", "drift2")),
                        paste("Item LPACF", c("drift1", "drift2")),
                        paste("Yen's $Q_3$", c("(drift1, drift1)", "(drift1, drift2)", "(drift2, drift2)")),
                        paste("OR", c("(drift1, drift1)", "(drift1, drift2)", "(drift2, drift2)")),
                        paste("OR difference", c("(drift1, drift1)", "(drift1, drift2)", "(drift2, drift2)")),
                        paste("RESID", c("(drift1, drift1)", "(drift1, drift2)", "(drift2, drift2)")),
                        paste("RESID  difference",c("(drift1, drift1)", "(drift1, drift2)", "(drift2, drift2)")))
row.names(outHDrift) <- paste0("\\hspace{0.25in}", row.names(outHDrift))

outHDrift <- outHDrift[-(1:6), ]

addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outHDrift))
addtorow$command <- c("& \\multicolumn{3}{c}{Number of Items} \\\\\n",
                      " & $I = 3$ & $I = 6$ & $I = 12$\\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} \\textit{drift1} denotes the items without parameter drift.
                      \\textit{drift2} denotes the items with parameter drift.}\n\\end{tablenotes}\n")

print(xtable(outHDrift, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-DPCM-HIPD as the Generating Model",
             label = "tab:halfdrift", align = c("l", "r", "r", "r")),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outHalfDrift.tex")

# Change of Meaning

pairsindex <- list(c(1, 2, 2),
                   rep(rep(1:2, 4), times = c(4, 1, 3, 1, 2, 1, 1, 2)),
                   rep(rep(1:2, 10), times = c(10, 1, 9, 1, 8, 1, 7, 1, 6, 1, 
                                               5, 1, 4, 1, 3, 1, 2, 1, 1, 2)))

outMeaning <- matrix(NA, 6 + 6 * 2 + 5 * 2, 3)

for(i in 1:3) {
  tmp.data <- Results[[i]][Results[[i]]$model == "Meaning", ]
  for (j in 1:6) {
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outMeaning[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
  rm(tmp)
  for (j in 7:12) {
    tmp1 <- unlist(tmp.data[, head(grep(discmeasures[j], names(tmp.data)), -1)])
    tmp2 <- unlist(tmp.data[, tail(grep(discmeasures[j], names(tmp.data)), 1)])
    outMeaning[2 * j - 7, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outMeaning[2 * j - 6, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
  }
  rm(tmp1, tmp2)
  for (j in 13:17) {
    tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 1]])
    tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 2]])
    outMeaning[2 * j - 7, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
    outMeaning[2 * j - 6, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
  }
  rm(tmp1, tmp2)
}

outMeaning <- ifelse(outMeaning >= 0.7,
                   paste0("\\textbf{", round(outMeaning, 2), "}"),
                   round(outMeaning, 2))

outMeaning <- data.frame(outMeaning)
names(outMeaning) <- paste0("I=", c(3, 6, 12))
row.names(outMeaning) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "PACF", "LPACF", "MSSD",
                          paste("Item-total correlation", c("meaning1", "meaning2")),
                          paste("Item-total correlation (v2)", c("meaning1", "meaning2")),
                          paste("Item-total correlation (v3)", c("meaning1", "meaning2")),
                          paste("Yen's $Q_1$", c("meaning1", "meaning2")),
                          paste("Yen's $Q_1$ alt.", c("meaning1", "meaning2")),
                          paste("Item LPACF", c("meaning1", "meaning2")),
                          paste("Yen's $Q_3$", c("(meaning1, meaning1)", "(meaning1, meaning2)")),
                          paste("OR", c("(meaning1, meaning1)", "(meaning1, meaning2)")),
                          paste("OR difference", c("(meaning1, meaning1)", "(meaning1, meaning2)")),
                          paste("RESID", c("(meaning1, meaning1)", "(meaning1, meaning2)")),
                          paste("RESID  difference",c("(meaning1, meaning1)", "(meaning1, meaning2)")))
row.names(outMeaning) <- paste0("\\hspace{0.25in}", row.names(outMeaning))

outMeaning <- outMeaning[-(1:6), ]

addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outMeaning))
addtorow$command <- c("& \\multicolumn{3}{c}{Number of Items} \\\\\n",
                      " & $I = 3$ & $I = 6$ & $I = 12$\\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} \\textit{meaning1} denotes the items for which its meaning did not changed.
                      \\textit{meaning2} denotes the item for which its meaning changed.}\n\\end{tablenotes}\n")

print(xtable(outMeaning, type = "latex", caption = "Proportion of Extreme PPP-Values with the TV-DPCM-Meaning as the Generating Model",
             label = "tab:Meaning", align = c("l", "r", "r", "r")),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outMeaning.tex")



# END