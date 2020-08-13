# IRT models to simulate data.
P.GRM <- function(C, IP, theta)
{
  N         <- length(theta)
  I         <- nrow(IP)
  alpha     <- IP[, ncol(IP)]
  betas     <- IP[, -ncol(IP)]
  res.cum   <- array(NA, c(N, I, C))
  for (i in 1:I)
  {
    for (c in 1:C)
    {
      arg              <- alpha[i] * (theta - betas[i, c])
      res.cum[ , i, c] <- exp(arg) / (1 + exp(arg))
    }
  }
  res.cum <- array(c(matrix(1, N, I), res.cum, matrix(0, N, I)), dim = c(N, I, C + 2))
  res     <- array(NA, c(N, I, C + 1))
  for (c in 1:(C + 1)) res[, , c] <- res.cum[, , c] - res.cum[, , c + 1]
  return(res)
}

P.GPCM <- function(y, alpha, delta, taus, theta, M)
{
  N         <- length(theta)
  I         <- length(delta)
  if (length(y) == 1) y <- rep(y, I)
  if (is.vector(taus)) taus <- matrix(rep(taus, I), nrow = I, byrow = TRUE)
  if (length(M) == 1) M <- rep(M, I)
  taus.zero <- cbind(0, taus)
  taus.cum  <- t(apply(taus.zero, 1, cumsum))
  part      <- function(i, w) 
  {
    if ((0 <= w) && (w <= M[i])) {
      exp(alpha[i] * (w * (theta - delta[i]) - taus.cum[i, ((max(M) - M[i])/2) + w + 1]))
    } else rep(0, N)
  }
  num       <- sapply(1:I,      function(i) part(i, y[i]), simplify = "array")
  tmp       <- sapply(0:max(M), function(w) sapply(1:I, function(i) part(i, w)))
  den       <- matrix(rowSums(tmp, na.rm = TRUE), ncol = I, byrow = FALSE)
  return(num / den)
}

gen.1PL <- function(theta, beta){
  probs <- t(sapply(theta, function(x) exp(x - beta)/(1 + exp(x - beta))))
  obs   <- t(apply(probs, 1, function(x) rbinom(length(beta), 1, x)))
  return(obs)
}

gen.2PL <- function(theta, alpha, beta){
  probs <- t(sapply(theta, function(x) exp(alpha * (x - beta))/(1 + exp(alpha * (x - beta)))))
  obs   <- t(apply(probs, 1, function(x) rbinom(length(beta), 1, x)))
  return(obs)
}

gen.3PL <- function(theta, alpha, beta, guess){
  probs <- t(sapply(theta, function(x) guess + (1 - guess) * exp(alpha * (x - beta))/(1 + exp(alpha * (x - beta)))))
  obs   <- t(apply(probs, 1, function(x) rbinom(length(beta), 1, x)))
  return(obs)
}