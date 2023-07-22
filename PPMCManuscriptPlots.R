# Plot for the Methods sections and the Empirical example.

# This script is used to make the plots of the manuscript about PPMC methods.
# This includes one figure that is used to introduce the PPMC methods and the 
# figures presented in the empirical example of the paper.

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()/4)
library(bayesplot)
color_scheme_set("darkgray")
library(xtable)
source("R/PPMC.R")
source("R/genTVDPCM.R")
source("R/tvdpcm2stan.R")
source("R/IRT_models.R")

# PPMC example ----

# Load Stan model.
model <- stan_model(file = "Stan/tv_dpcm_int_v5.1.stan", verbose = FALSE)

# Simulate TV-DPCM data and fit the TV-DPCM to the simulated data.
datasim <- gen.TVDPCM(nT = 200, I = 6, K = 5, seed = 2022, 
                      FUN = sinusoidal, maxAbsValue = 0.5)

standata <- tvdpcm2stan_data(resp = datasim$data,
                             I    = 6,
                             K    = 5,
                             nT   = 200, 
                             n_knots  = 8,
                             s_degree = 3)
tvdpcm_inits <- function() {
  list(lambda = runif(1, -1, 1),
       beta   = array(rnorm(6 * (5 - 1), 0, 3), dim = c(6, 5 - 1)),
       inno   = rnorm(200, 0, 3),
       sigma  = rlnorm(1, 1))
}

begin.time <- proc.time()
fit <- sampling(model,                  # Stan model. 
                data = standata,        # Data.
                iter = 2000,            # Number of iterations.
                chains  = 3,            # Number of chains.
                warmup  = 500,          # Burn-in samples.
                init    = tvdpcm_inits, # Initial values
                seed = 2022,            # Seed
                pars = c("beta", "theta", "lambda",
                         "sigma2", "pvar", "attractor", "rep_y"),
                control = list(adapt_delta   = 0.99,
                               max_treedepth = 15) # Other parameters to control sampling behavior.
) 
run.time <- proc.time() - begin.time
rm(begin.time)

# Verify model convergence.
stan.diag <- monitor(extract(fit,
                             pars = c("beta", "theta", "lambda", 
                                      "sigma2", "pvar", "attractor"),
                             permuted = FALSE, inc_warmup = FALSE), 
                     warmup = 0, print = FALSE)

ndiv  <- get_num_divergent(fit)     # number of divergent transitions
nbfmi <- get_low_bfmi_chains(fit)   # number of chains with a low bayesian fraction missing information 
ntree <- get_num_max_treedepth(fit) # number of transitions that exceeded the maximum treedepth
nbulk <- sum(stan.diag$Bulk_ESS < 100 * dim(fit)[2]) # number of parameters with low bulk ESS
ntail <- sum(stan.diag$Tail_ESS < 100 * dim(fit)[2]) # number of parameters with low tail ESS

#max Rhat
maxRhat <- round(max(rhat(fit, pars = c("beta", "theta", "lambda", 
                                        "sigma2", "pvar", "attractor"))), 4)
nRhat   <- sum(rhat(fit, pars = c("beta", "theta", "lambda",
                                  "sigma2", "pvar", "attractor")) > 1.05)

diag.table <- matrix(c(nRhat, ndiv, length(nbfmi), ntree, nbulk, ntail), nrow = 1)
row.names(diag.table)  <- "Fitted Model"
colnames(diag.table) <- c("Rhat>1.05", "N. Divergent", "N. Low BFMI", 
                          "Excedeed Treedepth", "Low Bulk ESS", "Low Tail ESS")

# Plot PPMC method by test statistic and by discrepancy measure.
pdf("Figures/PPMC_example.pdf", height = 4)
par(mfrow = c(1, 2), mar = c(4, 5, 0, 1) + 0.1, oma = c(0, 0, 2, 0) + 0.1)
ppmc.racf(object = fit, data = standata, 
          xlab = expression(T(y^rep)), ylim = c(0, 12), 
          col.yrep = rep(c("gray80", "gray30"), times = c(11, 3)))
ppmc.lpacf(object = fit, data = standata, quiet = TRUE, sumscores = TRUE,
           xlab = expression(D(y,omega)), ylab = expression(D(y^rep,omega)), 
           subtitle = FALSE, pch = 4, cex = 0.5, bty = "n",
           ylim = c(-0.3, 0.3), xlim = c(-0.3, 0.3), 
           col.abline = "black", col = c("gray30", "gray80"), split = TRUE)
dev.off()

# Example when the TV-DPCM does not fit the data ----

# Load saved fit and data.
standata <- readRDS("Fits/NegAMUwithNA_data.rds")
fit      <- readRDS("Fits/NegAMUwithNA.rds")

# MCMC Diagnostics
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

diag.table <- matrix(c(nRhat, ndiv, nbfmi, ntree, nbulk, ntail), nrow = 1)
row.names(diag.table)  <- "Fitted Model"
colnames(diag.table) <- c("Rhat>1.05", "N. Divergent", "N. Low BFMI", 
                          "Excedeed Treedepth", "Low Bulk ESS", "Low Tail ESS")

# Compute PPMC

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
ppmc11 <- ppmc.Q3(object = fit, data = standata)$sp
ppmc12 <- ppmc.OR(object = fit, data = standata)$sp
ppmc13 <- ppmc.ORDiff(object = fit, data = standata)$sp
ppmc14 <- ppmc.cov.resid(object = fit, data = standata)$sp
ppmc15 <- ppmc.cov.rediff(object = fit, data = standata)$sp

NAMUitemPPPs <- rbind(ppmc05, ppmc08, ppmc09, ppmc10)
colnames(NAMUitemPPPs) <- paste("Item ", 1:7)
row.names(NAMUitemPPPs) <- c("Item-Total Correlation", "Yen's $Q_1$", 
                             "Yen's $Q_1$ Alternative", "LPACF")
NAMUitemPPPs <- ifelse(NAMUitemPPPs <= 0.05 | NAMUitemPPPs >= 0.95,
                       paste0("\\textbf{", round(NAMUitemPPPs, 3), "}"),
                       round(NAMUitemPPPs, 3))

print(xtable(NAMUitemPPPs, type = "latex", caption = "PPPs Item-Level Measures Items of Negative Affect and Mental Unrest",
             label = "tab:namuitems", align = c("l", rep("c", 7))),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      file = "Tables/NAMUitemPPPs.tex")

pdf("Figures/PairsNAMU.pdf", heigh = 5)
ggpubr::ggarrange(ppmc11, ppmc12, ppmc13, ppmc14, ppmc15,
                  labels = c("Yen's Q3", "Odds Ratio", "OR Diff",
                             "RESID", "RESID Diff"), 
                  ncol = 3, nrow = 2, label.x = c(0.1, 0.05, 0.1, 0.125, 0.05))
dev.off()

# Example when the TV-DPCM fits the data ----

# Load saved fit and data.
standata <- readRDS("Fits/SelfEwithNA_3items_ph1-2_data.rds")
fit      <- readRDS("Fits/SelfEwithNA_3items_ph1-2.rds")

# MCMC Diagnostics
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

diag.table <- matrix(c(nRhat, ndiv, nbfmi, ntree, nbulk, ntail), nrow = 1)
row.names(diag.table)  <- "Fitted Model"
colnames(diag.table) <- c("Rhat>1.05", "N. Divergent", "N. Low BFMI", 
                          "Excedeed Treedepth", "Low Bulk ESS", "Low Tail ESS")

# Compute PPMC

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
ppmc11 <- ppmc.Q3(object = fit, data = standata)$sp
ppmc12 <- ppmc.OR(object = fit, data = standata)$sp
ppmc13 <- ppmc.ORDiff(object = fit, data = standata)$sp
ppmc14 <- ppmc.cov.resid(object = fit, data = standata)$sp
ppmc15 <- ppmc.cov.rediff(object = fit, data = standata)$sp

SelfEitemPPPs <- rbind(ppmc05, ppmc08, ppmc09, ppmc10)
colnames(SelfEitemPPPs) <- paste("Item ", 1:3)
row.names(SelfEitemPPPs) <- c("Item-Total Correlation", "Yen's $Q_1$", 
                             "Yen's $Q_1$ Alternative", "LPACF")
SelfEitemPPPs <- ifelse(SelfEitemPPPs <= 0.05 | SelfEitemPPPs >= 0.95,
                       paste0("\\textbf{", round(SelfEitemPPPs, 3), "}"),
                       round(SelfEitemPPPs, 3))

print(xtable(SelfEitemPPPs, type = "latex", caption = "PPPs Item-Level Measures Items of Self-Esteem",
             label = "tab:selfeitems", align = c("l", rep("c", 3))),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      file = "Tables/SelfEitemPPPs.tex")

pdf("Figures/PairsSelfE.pdf", heigh = 5)
ggpubr::ggarrange(ppmc11, ppmc12, ppmc13, ppmc14, ppmc15,
                  labels = c("Yen's Q3", "Odds Ratio", "OR Diff",
                             "RESID", "RESID Diff"), 
                  ncol = 3, nrow = 2, label.x = c(0.1, 0.05, 0.1, 0.125, 0.05))
dev.off()

# Example when the TV-DPCM fits the data ----

# Load saved fit and data.
standata <- readRDS("Fits/PosAwithNA_data.rds")
fit      <- readRDS("Fits/PosAwithNA.rds")

# MCMC Diagnostics
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

diag.table <- matrix(c(nRhat, ndiv, nbfmi, ntree, nbulk, ntail), nrow = 1)
row.names(diag.table)  <- "Fitted Model"
colnames(diag.table) <- c("Rhat>1.05", "N. Divergent", "N. Low BFMI", 
                          "Excedeed Treedepth", "Low Bulk ESS", "Low Tail ESS")

# Compute PPMC

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
ppmc11 <- ppmc.Q3(object = fit, data = standata)$sp
ppmc12 <- ppmc.OR(object = fit, data = standata)$sp
ppmc13 <- ppmc.ORDiff(object = fit, data = standata)$sp
ppmc14 <- ppmc.cov.resid(object = fit, data = standata)$sp
ppmc15 <- ppmc.cov.rediff(object = fit, data = standata)$sp

PosAitemPPPs <- rbind(ppmc05, ppmc08, ppmc09, ppmc10)
colnames(PosAitemPPPs) <- paste("Item ", 1:5)
row.names(PosAitemPPPs) <- c("Item-Total Correlation", "Yen's $Q_1$", 
                              "Yen's $Q_1$ Alternative", "LPACF")
PosAitemPPPs <- ifelse(PosAitemPPPs <= 0.05 | PosAitemPPPs >= 0.95,
                        paste0("\\textbf{", round(PosAitemPPPs, 3), "}"),
                        round(PosAitemPPPs, 3))

print(xtable(PosAitemPPPs, type = "latex", caption = "PPPs Item-Level Measures Items of Positive Affect",
             label = "tab:posaitems", align = c("l", rep("c", 5))),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", 
      sanitize.text.function = function(x){x}, booktabs = TRUE,
      file = "Tables/PosAitemPPPs.tex")

pdf("Figures/PairsPosA.pdf", heigh = 5)
ggpubr::ggarrange(ppmc11, ppmc12, ppmc13, ppmc14, ppmc15,
                  labels = c("Yen's Q3", "Odds Ratio", "OR Diff",
                             "RESID", "RESID Diff"), 
                  ncol = 3, nrow = 2, label.x = c(0.1, 0.05, 0.1, 0.125, 0.05))
dev.off()

