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


hist(unlist(results3[results3$model == "BiDim", grep("itcor3.I", names(results3))]))
hist(unlist(results6[results6$model == "TV-DPCM", grep("itcor.I", names(results6))]))
hist(unlist(results12[results12$model == "TV-DPCM", grep("itcor.I", names(results12))]))

modeltype <- unique(results3$model)
discmeasures <- c("ACF.1", "ACF.2", "ACF.3", "PACF", "LPACF", "MSSD",
                  "itcor.", "itcor2.", "itcor3.", "q1.", "q1alt.",
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

