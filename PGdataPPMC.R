# Run PPMC methods on all the fits of MU or SE
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() * (1/2))
# options(mc.cores = 2)
library(bayesplot)
color_scheme_set("darkgray")

source("R/IRT_models.R")
source("R/IRT_plots.R")
source("R/PPMC.R")
# source("R/PPMC2.R")
source("R/genTVDPCM.R")
source("R/tvdpcm2stan.R")

I <- 3

pairsindex <- apply(which(lower.tri(diag(I)), arr.ind = TRUE), 1, function(x)
  paste0("(", paste(x, collapse = ","), ")"))

filenames <- paste0("Fits/SelfEwithNA_3items", 
                    c("", "_ph1-2", "_ph1-3", "_ph1-4", "_ph2-3",
                      "_ph2-4", "_ph2-5", "_ph3-4", "_ph3-5",
                      "_ph4-5", "_ph3", "_ph4", "_ph5") )

out <- as.data.frame(matrix(NA, nrow = length(filenames), 
              ncol = 6 + 6 + I * 6 + 5 * ((I * (I - 1))/2)))
names(out) <-  c("nRhat", "ndiv", "nbfmi", "ntree", "nbulk", "ntail",
                 paste0("ppp-ac", 1:3), "ppp-pac", "ppp-lpac", "ppp-mssd",
                 paste0("ppp-itcor", 1:I), paste0("ppp-itcor2", 1:I), 
                 paste0("ppp-itcor3", 1:I), paste0("ppp-q1", 1:I),
                 paste0("ppp-q1alt", 1:I), paste0("ppp-ilpac", 1:I),
                 paste0("ppp-q3", pairsindex), paste0("ppp-or", pairsindex),
                 paste0("ppp-ordiff", pairsindex), paste0("ppp-resid", pairsindex),
                 paste0("ppp-rediff", pairsindex))
row.names(out) <- filenames

tictoc::tic()
for (i in 1:length(filenames)) {
  cat(paste0(filenames[i], "\n"))
  standata <- readRDS(paste0(filenames[i], "_data.rds"))
  fit      <- readRDS(paste0(filenames[i], ".rds"))
  
  stan.diag <- monitor(extract(fit, permuted = FALSE, inc_warmup = FALSE,
                               pars = c("beta", "theta", "lambda", 
                                        "sigma2", "pvar", "attractor")), 
                       warmup = 0, print = FALSE)
  
  nRhat   <- sum(rhat(fit, pars = c("beta", "theta", "lambda",
                                    "sigma2", "pvar", "attractor")) > 1.05)
  ndiv  <- get_num_divergent(fit)     # number of divergent transitions
  nbfmi <- get_low_bfmi_chains(fit)   # number of chains with a low bayesian fraction missing information 
  ntree <- get_num_max_treedepth(fit) # number of transitions that exceeded the maximum treedepth
  nbulk <- sum(stan.diag$Bulk_ESS < 100 * dim(fit)[2]) # number of parameters with low bulk ESS
  ntail <- sum(stan.diag$Tail_ESS < 100 * dim(fit)[2]) # number of parameters with low tail ESS
  
  if (length(nbfmi) == 0) {nbfmi <- 0}
  
  ppmc01 <- ppmc.acf(object = fit, data = standata, lag.max = 3)
  ppmc02 <- ppmc.racf(object = fit, data = standata)
  ppmc03 <- ppmc.lpacf(object = fit, data = standata, quiet = TRUE, sumscores = TRUE)
  ppmc04 <- ppmc.mssd(object = fit, data = standata)
  ppmc05 <- ppmc.itcor(object = fit, data = standata, quiet = TRUE, method = "pearson")
  ppmc06 <- ppmc.itcor2(object = fit, data = standata, quiet = TRUE, method = "pearson")
  ppmc07 <- ppmc.itcor3(object = fit, data = standata, quiet = TRUE)
  ppmc08 <- ppmc.Q1(object = fit, data = standata, quiet = TRUE)
  ppmc09 <- ppmc.Q1.alt(object = fit, data = standata, quiet = TRUE)
  ppmc10 <- ppmc.lpacf(object = fit, data = standata, quiet = TRUE)
  ppmc11 <- ppmc.Q3(object = fit, data = standata)$ppp
  ppmc12 <- ppmc.OR(object = fit, data = standata)$ppp
  ppmc13 <- ppmc.ORDiff(object = fit, data = standata)$ppp
  ppmc14 <- ppmc.cov.resid(object = fit, data = standata)$ppp
  ppmc15 <- ppmc.cov.rediff(object = fit, data = standata)$ppp
  
  out[i, ] <- c(nRhat, ndiv, nbfmi, ntree, nbulk, ntail, ppmc01, ppmc02,
                ppmc03, ppmc04, ppmc05, ppmc06, ppmc07, ppmc08, ppmc09, 
                ppmc10, ppmc11, ppmc12, ppmc13, ppmc14, ppmc15)
  rm(nRhat, ndiv, nbfmi, ntree, nbulk, ntail, ppmc01, ppmc02,
     ppmc03, ppmc04, ppmc05, ppmc06, ppmc07, ppmc08, ppmc09, 
     ppmc10, ppmc11, ppmc12, ppmc13, ppmc14, ppmc15, 
     stan.diag, standata, fit)
}
tictoc::toc()

# Checking the number of misfit given the computed ppp

outMU <- read.table("Fits/MU_PPMC.dat")
outSE <- read.table("Fits/SelfE_PPMC.dat")

apply(outMU[, -(1:6)], 1, function(x) sum(x <= 0.05 | x >= 0.95))
apply(outMU[, 1:6], 1, function(x) sum(x != 0))
apply(outSE[, -(1:6)], 1, function(x) sum(x <= 0.05 | x >= 0.95))
apply(outSE[, 1:6], 1, function(x) sum(x != 0))

datanames <- c("Fits/MUwithNA_data.rds", "Fits/MUwithNA_ph2-3_data.rds",
               "Fits/SelfEwithNA_3items_data.rds", "Fits/SelfEwithNA_3items_ph1-2_data.rds")

pdf(file = "Figures/RawdataMUandSelfE.pdf")
for (i in 1:length(datanames)) {
  tmp.data <- readRDS(datanames[i])
  
  tmp.matrix <- matrix(tmp.data$y_obs, 
                       nrow = length(unique(tmp.data$tt_obs)), 
                       ncol = tmp.data$I)
  
  # matplot(tmp.matrix[round(seq(1, nrow(tmp.matrix), length.out = 100)),],
  #         lty = 1, type = "l",
  #         ylab = "Responses")
  
  par(mfrow = c(2, 2), mar = c(2, 4, 1, 1) + 0.1) 
  window <- max(round(nrow(tmp.matrix)/10), 50)
  out <- c()
  
  for (j in window:nrow(tmp.matrix)) {
    tmp <- cor(tmp.matrix[(j - (window - 1)):j, ])[lower.tri(diag(ncol(tmp.matrix)))] 
    out <- rbind(out, tmp)
  }
  rm(j, tmp)
  
  matplot(out, lty = 1, type = "l", 
          ylab = "Correlation", ylim = c(-1, 1))
  
  
  window <- max(round(nrow(tmp.matrix)/10), 50)
  out <- c()
  
  for (j in window:nrow(tmp.matrix)) {
    tmp <- apply(tmp.matrix[(j - (window - 1)):j, ], 2, mean) 
    out <- rbind(out, tmp)
  }
  rm(j, tmp)
  
  matplot(out, lty = 1, type = "l", 
          ylab = "Mean")
  
  
  window <- max(round(nrow(tmp.matrix)/10), 50)
  out <- c()
  
  for (j in window:nrow(tmp.matrix)) {
    tmp <- apply(tmp.matrix[(j - (window - 1)):j, ], 2, var) 
    out <- rbind(out, tmp)
  }
  rm(j, tmp)
  
  matplot(out, lty = 1, type = "l", 
          ylab = "Variance")
  
  rpatterns <- sort(table(apply(tmp.matrix, 1, paste, collapse="")), decreasing = TRUE)
  
  plot(1:length(rpatterns), 
       as.vector(rpatterns),
       ylab = "Number of Times the RP was Observed",
       xaxt = "n",
       xlab = "Response Patterns",
       pch = "")
  text(1:length(rpatterns), 
       as.vector(rpatterns), 
       labels= names(rpatterns), 
       cex=0.75, font=0, adj = 0.5, srt = 90)
  
  par(mfrow = c(1, 1))
  psych::pairs.panels(tmp.matrix,
                      smooth   = FALSE,
                      density  = FALSE,
                      ellipses = FALSE,
                      cex.cor  = 0.5,
                      lm       = TRUE,
                      hist.col = gray(1/2),
                      rug      = FALSE,
                      labels = paste0("Item", 1:tmp.data$I))
  rm(tmp.data, tmp.matrix, window, rpatterns, out)
}
rm(i)
dev.off()

fit      <- readRDS(gsub("_data", "", datanames[4]))
standata <- readRDS(datanames[4])

sumscores <- tapply(standata$y_obs, standata$tt_obs, sum)

threshold <- summary(fit, pars = "beta")$summary
theta     <- summary(fit, pars = "theta")$summary
attractor <- summary(fit, pars = "attractor")$summary
sigma2    <- summary(fit, pars = "sigma2")$summary
pvar      <- summary(fit, pars = "pvar")$summary
lambda    <- summary(fit, pars = "lambda")$summary

thresholdT <- threshold[, 6]/sqrt(pvar[, 6])
thetaT     <- theta[, 6]/sqrt(pvar[, 6])
attractorT <- attractor[, 6]/sqrt(pvar[, 6])

plot(theta[unique(standata$tt_obs), 6], sumscores)
plot(thetaT[unique(standata$tt_obs)], sumscores)


plot(theta[, 6], type = "l")
lines(attractor[, 6], col = "red", lwd = 2)

plot(thetaT, type = "l")
lines(attractorT, col = "red", lwd = 2)

par(mfrow = c(I, 1), mar = c(4, 4, 2, 1) + 0.1)
plot.ICC(fit, standata, quiet = TRUE)
plot.ICC(fit, standata, quiet = TRUE, scale = TRUE)

par(mfrow = c(1, 1))
plot.IIF(fit, standata, type = "IIF")
plot.IIF(fit, standata, type = "IIF", scale = TRUE)

plot.IIF(fit, standata, type = "TIF")
plot.IIF(fit, standata, type = "TIF", scale = TRUE)

