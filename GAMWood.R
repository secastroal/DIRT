# This script follows the examples and solve the questions of the book:
# Generalized Additive Models: An Introduction with R. by Simon N. Wood.

# Chapter 4 ----

# Using piecewise linear basis (pg. 165)

require(gamair)
data(engine)

plot(engine$size, engine$wear, xlab = "Engine capacity", ylab = "Wear Index")

# Write R function to compute b_j(x) tent function basis for a given j

tf <- function(x, xj, j){
  dj <- xj * 0
  dj[j] <- 1
  approx(xj, dj, x)$y
}

# Write R function to compute the matrix for the piecewise linear function

tf.X <- function(x, xj){
  nk <- length(xj)
  n  <- length(x)
  X  <- matrix(NA, n, nk)
  for (j in 1:nk) {
    X[,j] <- tf(x, xj, j)
  }
  return(X)
}

sj <- seq(min(engine$size), max(engine$size), length.out = 6)
X  <- tf.X(engine$size, sj)
b  <- lm(engine$wear ~ X - 1)
s  <- seq(min(engine$size), max(engine$size), length.out = 200)
Xp <- tf.X(s, sj)
plot(engine$size, engine$wear, ylim = c(0, 5))
lines(s, Xp %*% coef(b))
matplot(s, Xp, add = TRUE, type = "l", lty = 2, col = rainbow(length(sj)))
matplot(s, t(t(Xp) * coef(b)), add = TRUE, type = "l", lty = 1, col = rainbow(length(sj)))

# Using penalization (pg. 169)

# Make function that augment the data and the model matrix and fits the penalized model.
prs.fit <- function(y, x, xj, sp) {
  X <- tf.X(x, xj)
  D <- diff(diag(length(xj)), differences = 2)
  X <- rbind(X, sqrt(sp) * D)
  y <- c(y, rep(0, nrow(D)))
  fit <- lm(y ~ X - 1)
  return(fit)
}

# fit penalized model with 20 knots

sj <- seq(min(engine$size), max(engine$size), length.out = 20)
b  <- prs.fit(engine$wear, engine$size, sj, sp = 10)
s  <- seq(min(engine$size), max(engine$size), length.out = 200)
Xp <- tf.X(s, sj)
plot(engine$size, engine$wear, ylim = c(0, 5))
lines(s, Xp %*% coef(b))
matplot(s, Xp, add = TRUE, type = "l", lty = 2, col = rainbow(length(sj)))
matplot(s, t(t(Xp) * coef(b)), add = TRUE, type = "l", lty = 1, col = rainbow(length(sj)))

# Selecting the adequate smoothing (penalized) parameter based on generalized cross validation.

rho <- seq(-9, 11, length = 90)
n   <- length(engine$wear)
V   <- rep(NA, 90)

# fitting the model with all the smoothing parameters rho

for (i in 1:90) {
  b    <- prs.fit(engine$wear, engine$size, sj, sp = exp(rho[i])) # Fit model
  trF  <- sum(influence(b)$hat[1:n])                              # extract effective df
  rss  <- sum((residuals(b)[1:n]) ^ 2)                            # Residual ss
  V[i] <- (n * rss)/((n - trF) ^ 2)                               # GCV score
}

# plot GCV scores and optimal fit

plot(rho, V, type = "l",
     xlab = expression(log(lambda)), main = "GCV score")

sp <- exp(rho[V == min(V)])                          # get optimal lambda
b  <- prs.fit(engine$wear, engine$size, sj, sp = sp) # re fit the model
Xp <- tf.X(s, sj)

plot(engine$size, engine$wear, ylim = c(0, 5), main = "GCV optimal fit")
lines(s, Xp %*% coef(b))
matplot(s, Xp, add = TRUE, type = "l", lty = 2, col = rainbow(length(sj)))
matplot(s, t(t(Xp) * coef(b)), add = TRUE, type = "l", lty = 1, col = rainbow(length(sj)))

# Bayesian penalized reparameterization

X0      <- tf.X(engine$size, sj)                       # Original parameterization
D       <- rbind(0, 0, diff(diag(20), difference = 2)) # Augmented D
diag(D) <- 1
X       <- t(backsolve(t(D), t(X0)))                   # reparameterized X
Z       <- X[, -c(1, 2)]
X       <- X[, c(1, 2)]

# Create function to minimize in the gradient optim

llm <- function(theta, X, Z, y) {
  # Untransform parameters
  sigma.b <- exp(theta[1])
  sigma   <- exp(theta[2])
  # Extract dimensions
  n  <- length(y)
  pr <- ncol(Z)
  pf <- ncol(X)
  # obtain \hat \beta, \hat b
  X1   <- cbind(X, Z)
  ipsi <- c(rep(0, pf), rep(1/(sigma.b ^ 2), pr))
  b1   <-  solve(crossprod(X1)/(sigma ^ 2) + diag(ipsi),
                 t(X1)%*%(y/(sigma ^ 2)))
  # compute log|Z'Z
  ldet <- sum(log(diag(chol(crossprod(Z)/(sigma ^ 2) + diag(ipsi[-(1:pf)])))))
  # compute log profile likelihood
  l <- (-sum((y - X1 %*% b1) ^ 2)/(sigma ^ 2) - sum(b1 ^ 2 * ipsi) -
          n * log(sigma ^2) - pr * log(sigma.b ^ 2) - 2 * ldet -  
          n * log(2 * pi))/2
  attr(l, "b") <- as.numeric(b1)
  -l
}

# Estimating smoothing and variance parameters with ML and REML
m <- optim(c(0, 0), llm, method = "BFGS", X = X, Z = Z, y = engine$wear)
b <- attr(llm(m$par, X, Z, engine$wear), "b")

g <- factor(rep(1, nrow(X)))
m <- nlme::lme(wear ~ X - 1, random = list(g = nlme::pdIdent(~ Z - 1)), data = engine)

plot(engine$size, engine$wear, ylim = c(0, 5))
Xp1 <- t(backsolve(t(D), t(Xp)))
lines(s, Xp1 %*% as.numeric(b), col = "grey", lwd = 2)
lines(s, Xp1 %*% as.numeric(coef(m)), col = "black", lwd = 2)

# Additive models 
# y = a + f(x) + g(z) + e, with f() and g() smooth functions (Equation 4.8 pg. 174)

# When multiple smoothers, there might be identifiability issues, reason why
# constraints are needed.Most commonly sum-to-zero constraints (pg. 175).

# Function to create constrained version of X and D
tf.XD <- function(x, xk, cmx = NULL, m = 2) {
  # get X and D subject to constraint
  nk <- length(xk)
  X  <-  tf.X(x, xk)[, -nk]                    # constrained basis matrix
  D  <- diff(diag(nk), differences = m)[, -nk] # constrained D matrix
  if (is.null(cmx)) {cmx <- colMeans(X)}
  X <- sweep(X, 2, cmx)
  out <-  list(X = X, D = D, cmx = cmx)
  return(out)
}

# Function to fit model with 2 smoothers

am.fit <- function(y, x, v, sp, k = 10) {
  # setup bases and penalties
  xk  <- seq(min(x), max(x), length.out = k)
  xdx <- tf.XD(x, xk)
  vk  <- seq(min(v), max(v), length.out = k)
  xdv <- tf.XD(v, vk)
  # Create augmented model matrix and response
  nD <- nrow(xdx$D) * 2
  sp <- sqrt(sp)
  X  <- cbind(c(rep(1, nrow(xdx$X)), rep(0, nD)),
              rbind(xdx$X, sp[1] * xdx$D, xdv$D * 0),
              rbind(xdv$X, xdx$D * 0, sp[2] * xdv$D))
  y1 <- c(y, rep(0, nD))
  # fit model
  b <- lm(y1 ~ X - 1)
  # Compute useful quantities
  n   <- length(y)
  trA <- sum(influence(b)$hat[1:n])
  rds <- residuals(b)
  rss <- sum(rds ^ 2)
  sig.hat <- rss/(n - trA)
  gcv <- (sig.hat * n)/(n - trA)
  Vb  <- (vcov(b) * sig.hat)/(summary(b)$sigma ^ 2)
  # Return fitted model
  out <- list(b = coef(b), Vb = Vb, edf = trA, gcv = gcv, fitted = fitted(b)[1:n],
              rds = rds, xk = list(xk, vk), cmx = list(xdx$cmx, xdv$cmx))
  return(out)
}


# Fit additive model with two smoothers using optim

am.gcv <- function(lsp, y, x, v, k) {
  ## function suitable for GCV optimization with optim
  am.fit(y, x, v, exp(lsp), k)$gcv
}

fit <- optim(c(0, 0), am.gcv, y = trees$Volume, x = trees$Girth,
             v = trees$Height, k = 10)
sp <- exp(fit$par)
# get fit with optimal smoothing parameters
fit <- am.fit(trees$Volume, trees$Girth, trees$Height, sp, k = 10)

# function to plot smooth effects

am.plot <-  function(fit, xlab, ylab) {
  start <- 2 # where smooth coefficients start in beta
  for (i in 1:2) {
    # sequence of values to predict
    x <-  seq(min(fit$xk[[i]]), max(fit$xk[[i]]), length.out = 200)
    # prediction matrix
    Xp <- tf.XD(x, fit$xk[[i]], fit$cmx[[i]])$X
    # extract coefficients
    stop <- start + ncol(Xp) - 1
    ind  <- start:stop
    b    <- fit$b[ind]
    Vb   <- fit$Vb[ind, ind]
    # values for smooth
    fv <-  Xp %*% b
    # standard errors of smooth at x
    se <- sqrt(rowSums((Xp %*% Vb) * Xp))
    # 2 s.e. limits
    ul <- fv + 2 * se
    ll <- fv - 2 * se
    # plot smooth and limits
    plot(x, fv, type = "l", ylim = range(c(ll, ul)), 
         xlab = xlab[i], ylab = ylab[i])
    lines(x, ul, lty = 2)
    lines(x, ll, lty = 2)
    start <- stop + 1
  }
}

# plotting the fitted model for trees

par(mfrow = c(1, 3))

plot(fit$fitted, trees$Volume, xlab = "Fitted volume",
     ylab = "Observed volume")
am.plot(fit, xlab = c("Girth", "Height"), 
        ylab = c("s(Girth)", "s(Height)"))

# Fitting the generalized additive model (multiple predictors and link functions)

gam.fit <-  function(y, x, v, sp, k = 10) {
  # gamma error log link 2 term gam fit
  eta <-  log(y) # initial eta
  not.converged <-  TRUE
  old.gcv <- -100
  
  while (not.converged) {
    mu  <-  exp(eta)               # current mu estimate
    z   <-  (y - mu)/mu + eta      # pseudodata
    fit <-  am.fit(z, x, v, sp, k) # Penalized least squares
    if (abs(fit$gcv - old.gcv) < 1e-5 * fit$gcv) {
      not.converged <-  FALSE
    }
    old.gcv <-  fit$gcv
    eta     <-  fit$fitted # update eta
  }
  fit$fitted <- exp(fit$fitted) # mu
  return(fit)
}

gam.gcv <- function(lsp, y, x, v, k = 10) {
  gam.fit(y, x, v, exp(lsp), k = k)$gcv
}

fit <-  optim(c(0, 0), gam.gcv, y = trees$Volume, x = trees$Girth,
              v = trees$Height, k = 10)
sp  <- exp(fit$par)
fit <- gam.fit(trees$Volume, trees$Girth, trees$Height, sp)

par(mfrow = c(1, 3))
plot(fit$fitted, trees$Volume, xlab = "Fitted volume",
     ylab = "Observed volume")
am.plot(fit, xlab = c("Girth", "Height"), 
        ylab = c("s(Girth)", "s(Height)"))



theta <- rnorm(1)
K <- 10e5
beta  <- sort(rnorm(K))
beta  <- beta - mean(beta)



system.time({
  cumbeta <- cumsum(beta)
  
  probs <- rep(NA, K+1)
probs[1] <- 0
for (i in 2:(K + 1)) {
  probs[i] <- (i - 1) * theta - cumbeta[i - 1]
}}
)

system.time(probs2 <- c(0, cumsum(theta - beta)))

all.equal(probs, probs2)






