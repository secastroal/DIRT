# Here I test each function in 'IRT_models.R', to make sure they are all ok.

library(mirt)
source("R/IRT_models.R")

# Test gen.1PL() ----
set.seed(1)
N                      <- 1000
I                      <- 50
theta                  <- rnorm(N)
beta                   <- rnorm(I)
gen.data.1pl           <- gen.1PL(theta, beta)
colnames(gen.data.1pl) <- paste0("It", 1:I)
# 
mirt.1pl       <- mirt(gen.data.1pl, model = 1, itemtype = "Rasch")
mirt.1pl.items <- coef(mirt.1pl, IRTpars = TRUE, simplify = TRUE)
mirt.1pl.theta <- fscores(mirt.1pl)
# 
plot(beta, mirt.1pl.items$items[,"b"], pch = 4, 
     main = paste0("Difficulty; cor = ", round(cor(beta, mirt.1pl.items$items[,"b"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
plot(theta, mirt.1pl.theta, pch = 4, 
     main = paste0("Thetas; cor = ", round(cor(theta, mirt.1pl.theta), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
rm(list = setdiff(ls(), lsf.str()))

# Test gen.2PL() ----
set.seed(2)
N                      <- 1000
I                      <- 50
theta                  <- rnorm(N)
alpha                  <- runif(I, .5, 3)
beta                   <- rnorm(I)
gen.data.2pl           <- gen.2PL(theta, alpha, beta)
colnames(gen.data.2pl) <- paste0("It", 1:I)
# 
mirt.2pl       <- mirt(gen.data.2pl, model = 1, itemtype = "2PL")
mirt.2pl.items <- coef(mirt.2pl, IRTpars = TRUE, simplify = TRUE)
mirt.2pl.theta <- fscores(mirt.2pl)
# 
plot(alpha, mirt.2pl.items$items[,"a"], pch = 4, 
     main = paste0("Discrimination; cor = ", round(cor(alpha, mirt.2pl.items$items[,"a"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
plot(beta, mirt.2pl.items$items[,"b"], pch = 4, 
     main = paste0("Difficulty; cor = ", round(cor(beta, mirt.2pl.items$items[,"b"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
plot(theta, mirt.2pl.theta, pch = 4, 
     main = paste0("Thetas; cor = ", round(cor(theta, mirt.2pl.theta), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
rm(list = setdiff(ls(), lsf.str()))

# Test gen.3PL() ----
set.seed(3)
N                      <- 1000
I                      <- 50
theta                  <- rnorm(N)
alpha                  <- runif(I, .5, 3)
beta                   <- rnorm(I)
guess                  <- runif(I, 0, .25)
gen.data.3pl           <- gen.3PL(theta, alpha, beta, guess)
colnames(gen.data.3pl) <- paste0("It", 1:I)
# 
mirt.3pl       <- mirt(gen.data.3pl, model = 1, itemtype = "3PL")
mirt.3pl.items <- coef(mirt.3pl, IRTpars = TRUE, simplify = TRUE)
mirt.3pl.theta <- fscores(mirt.3pl)
# 
plot(alpha, mirt.3pl.items$items[,"a"], pch = 4, 
     main = paste0("Discrimination; cor = ", round(cor(alpha, mirt.3pl.items$items[,"a"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
plot(beta, mirt.3pl.items$items[,"b"], pch = 4, 
     main = paste0("Difficulty; cor = ", round(cor(beta, mirt.3pl.items$items[,"b"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
plot(guess, mirt.3pl.items$items[,"g"], pch = 4, 
     main = paste0("Guessing; cor = ", round(cor(guess, mirt.3pl.items$items[,"g"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
plot(theta, mirt.3pl.theta, pch = 4, 
     main = paste0("Thetas; cor = ", round(cor(theta, mirt.3pl.theta), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
rm(list = setdiff(ls(), lsf.str()))

# Test P.GRM() ----
set.seed(4)
N             <- 1000
I             <- 50
K             <- 5
theta         <- rnorm(N)
# I use Sebastian's code here, with some edits:
alpha         <- rlnorm(I, 0, 0.25)    # Discrimination parameters.
delta         <- matrix(NA, I, K - 1)  # Matrix to store difficulty parameters.
delta[, 1]    <- rnorm(I, -2)
for (i in 2:(K - 1)){
  delta[, i] <- delta[, i - 1] + runif(I, 0.4, 1.2) #!#!# edited!
}
delta <- delta - rowMeans(delta)                    #!#!# edited!
delta <- delta + rnorm(I)                           #!#!# edited!
IP           <- cbind(delta, alpha) # Matrix of item parameters to input in P.GRM
colnames(IP) <- c(paste0("delta", 1:(K-1)), "alpha")
# 
probs.array            <- P.GRM(K - 1, IP, theta)
gen.data.grm           <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
colnames(gen.data.grm) <- paste0("It", 1:I)
rm(probs.array, delta, alpha)
# 
mirt.grm       <- mirt(gen.data.grm, model = 1, itemtype = "graded")
mirt.grm.items <- coef(mirt.grm, IRTpars = TRUE, simplify = TRUE)
mirt.grm.theta <- fscores(mirt.grm)
# 
plot(IP[, "alpha"], mirt.grm.items$items[, "a"], pch = 4, 
     main = paste0("Discrimination; cor = ", round(cor(IP[, "alpha"], mirt.grm.items$items[, "a"]), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
for (i in 1:(K-1))
{
  plot(IP[, paste0("delta", i)], mirt.grm.items$items[, paste0("b", i)], pch = 4, 
       main = paste0("cor = ", round(cor(IP[, paste0("delta", i)], mirt.grm.items$items[, paste0("b", i)]), 4)))
  abline(0, 1, col = 2, lwd = 2)
}
# 
plot(theta, mirt.grm.theta, pch = 4, 
     main = paste0("Thetas; cor = ", round(cor(theta, mirt.grm.theta), 3)))
abline(0, 1, col = 2, lwd = 2)
# 
rm(list = setdiff(ls(), lsf.str()))

# Test P.GPCM() ----
# Not today.
