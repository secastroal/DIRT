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

# Discrepancy measures variable names
discmeasures <- c("acf1", "acf2", "acf3", "racf", "lpacf", "mssd",
                  "itcor", "q1_", "q1alt_", "ilpacf_", "q3_", "or_", 
                  "ordiff_", "resid_", "rediff_")

# 1.0 Read output files ----

# Read output files into R
files <- paste0(getwd(), "/Simulation/Sim_PPMC_cond_", 1:nrow(Cond), ".txt")

data <- lapply(files, function(x) read.table(file = x, header = TRUE))

# Merge output into one data.frame
results <- ldply(data, data.frame)

rm(data, files)

# 2.0 Check convergence per condition ----

# Number of replications that did not converged
conv <- tapply(results$corrupt, results$cond, sum)

# In almost all of the cases it was due to divergent transitions.
# There is only one of these replications in which the Rhats were larger than 
# 1.05.

results[results$maxRhat > 1.05, c(1:13)]
results[results$ndiv != 0, c(1:13)]
table(results[results$ndiv != 0, c(1:13)]$ndiv)
results[results$nbfmi != 0, c(1:13)]

pdf("Figures/ppmc_conv.pdf", height = 4)
par(mar = c(5, 9, 0.25, 1) + 0.1)
bp <- barplot(rev(100 - conv), xlim = c(0, 100), las = 1, horiz = TRUE,
              names.arg = paste(rev(c("TV-DPCM", "TV-MDPCM", "TV-MDPCM", 
                                  "TV-DGPCM", "TV-DGPCM", "TV-DPCM-IPD",
                                  "TV-DPCM-IPD", "TV-DPCM-Meaning")), 
                                rev(c("", c("r = 0.3", "r = 0.6"), 
                                  rep(c("1/3", "2/3"), 2), ""))),
              xlab = "Convergence Rate (%)")
mtext(rev(paste0(100 - conv, "%")), side = 2, line = -2, las = 1,
      at = as.vector(bp))
rm(bp)
dev.off()

# Create and optional matrix to exclude replications that did not converged.

# results_conv <- results
results_conv <- results[results$corrupt == 0, ]

# 3.0 Plot PPPs distributions ----

# Plot the distribution of the PPP-values

discme_labels <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "RACF", "LRACF", "MSSD",
                   "Item-total Correlation", expression(paste("Yen's ", Q[1])), 
                   expression(paste("Yen's ", Q[1], " alt.")),
                   "Item LRACF", expression(paste("Yen's ", Q[3])), 
                   "OR", "OR difference", 
                   "RESID", "RESID  difference") 

pdf("Figures/ppp_dist_test.pdf", height = 5)
par(mfrow = c(2, 3), mar = c(2, 4, 4, 1) + 0.1, oma = c(0, 1, 0, 1) + 0.1)
for (i in 1:6) {
  hist(unlist(results_conv[results_conv$genModel == "TV-DPCM", 
                           grep(discmeasures[i], names(results_conv))]),
       las = 1, xlim = c(0, 1), main = discme_labels[i], xlab = "", 
       freq = FALSE, breaks = seq(0, 1, by = 0.05), 
       col = c("black", rep("lightgray", 18), "black"), 
       border = c("black", rep("lightgray", 18), "black"))
}
rm(i)
dev.off()

pdf("Figures/ppp_dist_item.pdf", height = 5)
par(mfrow = c(2, 2), mar = c(2, 4, 4, 1) + 0.1, oma = c(0, 1, 0, 1) + 0.1)
for (i in 7:10) {
  hist(unlist(results_conv[results_conv$genModel == "TV-DPCM", 
                           grep(discmeasures[i], names(results_conv))]),
       las = 1, xlim = c(0, 1), main = discme_labels[i], xlab = "", 
       freq = FALSE, breaks = seq(0, 1, by = 0.05), 
       col = c("black", rep("lightgray", 18), "black"), 
       border = c("black", rep("lightgray", 18), "black"))
}
rm(i)
dev.off()

pdf("Figures/ppp_dist_pair.pdf", height = 5)
par(mfrow = c(2, 3), mar = c(2, 4, 4, 1) + 0.1, oma = c(0, 1, 0, 1) + 0.1)
for (i in 11:15) {
  hist(unlist(results_conv[results_conv$genModel == "TV-DPCM", 
                           grep(discmeasures[i], names(results_conv))]),
       las = 1, xlim = c(0, 1), main = discme_labels[i], xlab = "", 
       freq = FALSE, breaks = seq(0, 1, by = 0.05), 
       col = c("black", rep("lightgray", 18), "black"), 
       border = c("black", rep("lightgray", 18), "black"))
}
rm(i)
dev.off()

# 4.0 Tables of measures' power  ----

# Create empty table to store results from all conditions

outresults <- matrix(NA, nrow = 6 + 4 * 2 + 5 * 3, ncol = nrow(Cond))

# Create indexes to compute power pair pairs of interest in the
# pairwise measures. For example, in the bidimiensional conditions,
# we would like to the the proportion of misfit when using items
# from the same dimension versus items from different dimensions.

# Pairs are built different depending on the condition. 
# In BiDim, there are two factors with three items each.
# In GPCM and Drift, either the two or the four
# last items violate the assumptions of the TV-DPCM.
# In Meaning, one item changes it meaning.

pairsindex <- list(rep(1, 15),
                   rep(c(1, 2, 1, 2, 3), times = c(2, 3, 1, 6, 3)),
                   rep(c(rep(1:2, 3), 3), times = c(3, 2, 2, 2, 1, 4, 1)),
                   rep(1:3, times = c(1, 8, 6)),
                   rep(rep(1:2, 4), times = c(4, 1, 3, 1, 2, 1, 1, 2)))

pairsindex <- pairsindex[c(1, 2, 2, 3, 4, 3, 4, 5)]

# It goes similar for the items.

itemsindex <- list(rep(1, 6),
                   rep(1:2, each = 3),
                   rep(1:2, times = c(4, 2)),
                   rep(1:2, times = c(2, 4)),
                   rep(1:2, times = c(5, 1)))

itemsindex <- itemsindex[c(1, 2, 2, 3, 4, 3, 4, 5)]

# For loop to compute measures' power and fill in the matrix outresults.

for (i in 1:nrow(Cond)) {
  # Select replications for condition i
  tmp.data <- results_conv[results_conv$cond == i, ]
  
  # Compute power of test-level discrepancy measures
  for (j in 1:6) {
    tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
    outresults[j, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
  }
  rm(j, tmp)
  
  # Compute power of item-level discrepancy measures
  for (j in 7:10) {
    if (i == 1) {
      tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
      outresults[2 * j - 7, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
      rm(tmp)
    } else {
      tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[itemsindex[[i]] == 1]])
      tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[itemsindex[[i]] == 2]])
      outresults[2 * j - 7, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
      outresults[2 * j - 6, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
      rm(tmp1, tmp2)
    }
  }
  rm(j)
  
  # Compute power of pair-wise discrepancy measures
  for (j in 11:15) {
    if (i == 1) {
      tmp <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))])
      outresults[3 * j - 18, i] <- round(sum(tmp <= 0.05 | tmp >= 0.95)/length(tmp), 2)
      rm(tmp)
    }
    
    if (i %in% 2:7) {
      tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 1]])
      tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 2]])
      tmp3 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 3]])
      outresults[3 * j - 18, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
      outresults[3 * j - 17, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
      outresults[3 * j - 16, i] <- round(sum(tmp3 <= 0.05 | tmp3 >= 0.95)/length(tmp3), 2)
      rm(tmp1, tmp2, tmp3)
    }
    
    if (i == 8) {
      tmp1 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 1]])
      tmp2 <- unlist(tmp.data[, grep(discmeasures[j], names(tmp.data))[pairsindex[[i]] == 2]])
      outresults[3 * j - 18, i] <- round(sum(tmp1 <= 0.05 | tmp1 >= 0.95)/length(tmp1), 2)
      outresults[3 * j - 17, i] <- round(sum(tmp2 <= 0.05 | tmp2 >= 0.95)/length(tmp2), 2)
      rm(tmp1, tmp2)
    }
  }
  rm(j)
}
rm(i, tmp.data)

# Highlight in boldface power larger than 0.7

outresults <- ifelse(outresults >= 0.7,
                     paste0("\\textbf{", round(outresults, 2), "}"),
                     round(outresults, 2))

# Turn into data frame and name rows and columns
outresults <- data.frame(outresults)
names(outresults) <- Cond$genModel
row.names(outresults) <- c("ACF lag 1", "ACF lag 2", "ACF lag 3", "RACF", "LRACF", "MSSD",
                           paste("Item-total correlation", c("set1", "set2")),
                           paste("Yen's $Q_1$", c("set1", "set2")),
                           paste("Yen's $Q_1$ alt.", c("set1", "set2")),
                           paste("Item LRACF", c("set1", "set2")),
                           paste("Yen's $Q_3$", c("(set1, set1)", "(set1, set2)", "(set2, set2)")),
                           paste("OR", c("(set1, set1)", "(set1, set2)", "(set2, set2)")),
                           paste("OR difference", c("(set1, set1)", "(set1, set2)", "(set2, set2)")),
                           paste("RESID", c("(set1, set1)", "(set1, set2)", "(set2, set2)")),
                           paste("RESID  difference",c("(set1, set1)", "(set1, set2)", "(set2, set2)")))
row.names(outresults) <- paste0("\\hspace{0.25in}", row.names(outresults))

outresultsfull <- outresults

# Export table for manuscript.
# Remove row with test-level discrepancy measures and alternative OR and RESID.
outresults <- outresults[-c(1:6, 21:23, 27:29), ]

# Adjust heading and table note.
addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outresults))
addtorow$command <- c("& TV-DPCM & \\multicolumn{2}{c}{TV-MDPCM} & \\multicolumn{2}{c}{TV-DGPCM} 
                      & \\multicolumn{2}{c}{TV-DPCM-IPD} & TV-DPCM \\\\\n",
                      " & & $r = 0.3$ & $r = 0.6$ & $1/3$ & $2/3$ & $1/3$ & $2/3$ & Meaning \\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} In conditions where the generating model was the TV-MDPCM, 
                      \\textit{set1} denotes the items of dimension 1 and \\textit{set2} denotes 
                      the items of dimension 2. In conditions where the generating model was the 
                      TV-DGPCM, \\textit{set1} denotes the items with discrimination parameters 
                      equal to 1 and \\textit{set2} denotes the items with discrimination
                      parameters different from 1. In conditions where the generating model was 
                      the TV-DPCM-IPD, \\textit{set1} denotes the items that do not present
                      item parameter drift and \\textit{set2} denotes the items that have item
                      parameter drift. In conditions where the generating model was the 
                      TV-DPCM-Meaning, \\textit{set1} denotes the items for which its meaning 
                      did not change and \\textit{set2} denotes the item for which its meaning 
                      changed. Proportions larger than 0.7 are highlighted in boldface.}\n
                      \\end{tablenotes}\n")

print(xtable(outresults, type = "latex", caption = "Proportion of Extreme PPP-Values across Conditions",
             label = "tab:outresults", align = c("l", rep("r", 8))),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outResults.tex")

# Export table for appendix.
# Select rows with ineffective measures.
outineffective <- outresultsfull[c(1:6, 21:23, 27:29), ]

# Adjust heading and table note.
addtorow <- list()
addtorow$pos <- list(0, 0, nrow(outineffective))
addtorow$command <- c("& TV-DPCM & \\multicolumn{2}{c}{TV-MDPCM} & \\multicolumn{2}{c}{TV-DGPCM} 
                      & \\multicolumn{2}{c}{TV-DPCM-IPD} & TV-DPCM \\\\\n",
                      " & & $r = 0.3$ & $r = 0.6$ & $1/3$ & $2/3$ & $1/3$ & $2/3$ & Meaning \\\\\n",
                      "\\begin{tablenotes}[para,flushleft]\n{\\small 
                      \\textit{Note.} In conditions where the generating model was the TV-MDPCM, 
                      \\textit{set1} denotes the items of dimension 1 and \\textit{set2} denotes 
                      the items of dimension 2. In conditions where the generating model was the 
                      TV-DGPCM, \\textit{set1} denotes the items with discrimination parameters 
                      equal to 1 and \\textit{set2} denotes the items with discrimination
                      parameters different from 1. In conditions where the generating model was 
                      the TV-DPCM-IPD, \\textit{set1} denotes the items that do not present
                      item parameter drift and \\textit{set2} denotes the items that have item
                      parameter drift. In conditions where the generating model was the 
                      TV-DPCM-Meaning, \\textit{set1} denotes the items for which its meaning 
                      did not change and \\textit{set2} denotes the item for which its meaning 
                      changed.}\n
                      \\end{tablenotes}\n")

print(xtable(outineffective, type = "latex", caption = "Proportion of Extreme PPP-Values across Conditions",
             label = "tab:outineffective", align = c("l", rep("r", 8))),
      include.colnames = FALSE, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      add.to.row = addtorow,
      file = "Tables/outIneffective.tex")

# END