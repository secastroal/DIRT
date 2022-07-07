#JT# Please offer the option to only plot the ICCs for the desired items.
#JT# I actually wanted to plot one item at a time, to compare with ltm's ICCs.
plot.ICC <- function(object, data, range = c(-5, 5), items = NULL, 
                     item_labels = NULL, quiet = FALSE, scale = FALSE, ...) {
  
  tmp      <- list()
  tmp$beta <- summary(object, pars = "beta")$summary
  if (scale) {
    tmp$pvar <- summary(object, pars = "pvar")$summary
  }
  
  I <- data$I
  K <- data$K
  M <- data$K - 1
  
  if (is.null(item_labels)) {
    item_labels <- paste0("It", 1:I) 
  }
  
  if (is.null(items)) {
    items <- 1:I
  }
  
  theta <- seq(range[1], range[2], length = 300)
  
  if (scale) {
    thresholds <- matrix(tmp$beta[, 6]/sqrt(tmp$pvar[, 6]), nrow = I, ncol = M, byrow = TRUE)
    alpha      <- rep(sqrt(tmp$pvar[, 6]), I)
  } else {
    thresholds <- matrix(tmp$beta[, 6], nrow = I, ncol = M, byrow = TRUE)
    alpha      <- rep(1, I) 
  }
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (y in 0:M) {
    probs.array[, , y + 1] <- P.GPCM(y     = y, 
                                     alpha = alpha, 
                                     delta = delta, 
                                     taus  = taus, 
                                     theta = theta, 
                                     M     = M)
  }
    
  for (i in 1:length(items)) {
    if (!quiet) {invisible(readline(prompt="Press [enter] to continue"))}
    matplot(theta, probs.array[, items[i], ], type = "l", 
            ylim = c(0, 1), 
            ylab = expression(paste("P(", X[i], "|", theta[t], ")")), 
            xlab = expression(theta[t]),
            main = paste0("Item ", items[i], ": ", item_labels[items[i]]), ...)
  } 
}


plot.IIF <- function(object, data, range = c(-5, 5), item_labels = NULL, 
                     type = c("IIF", "TIF"), legend = TRUE, 
                     scale = FALSE,  ...) {
  
  tmp      <- list()
  tmp$beta <- summary(object, pars = "beta")$summary
  if (scale) {
    tmp$pvar <- summary(object, pars = "pvar")$summary
  }
  
  I <- data$I
  K <- data$K
  M <- data$K - 1
  
  if (is.null(item_labels)) {
   item_labels <- paste0("It", 1:I) 
  }
  
  theta <- seq(range[1], range[2], length = 300)
  
  if (scale) {
    thresholds <- matrix(tmp$beta[, 6]/sqrt(tmp$pvar[, 6]), nrow = I, ncol = M, byrow = TRUE)
    alpha      <- rep(sqrt(tmp$pvar[, 6]), I)
  } else {
    thresholds <- matrix(tmp$beta[, 1], nrow = I, ncol = M, byrow = TRUE)
    alpha      <- rep(1, I)
  }
  delta       <- rowMeans(thresholds)
  taus        <- thresholds - delta
  
  probs.array <- array(NA, dim = c(length(theta), I, K))
  
  for (y in 0:M) {
    probs.array[, , y + 1] <- P.GPCM(y     = y, 
                                     alpha = alpha, 
                                     delta = delta, 
                                     taus  = taus, 
                                     theta = theta, 
                                     M     = M)
  }
  
  Info <- matrix(NA, length(theta), I)
  
  for (i in 1:I) {
    T_bar     <- colSums(t(probs.array[, i,]) * 1:K)
    Info[, i] <- colSums((outer(1:K, T_bar, "-") ^ 2) * t(probs.array[, i,]))
  }
  
  if (type == "IIF") {
    matplot(theta, Info, type = "l", las = 1,
            ylab = "Information", 
            xlab = expression(theta[t]),
            main = "Item Information Functions", ...)
    if (legend) {
      legend("topright", item_labels, ...)
      }
  }
  
  if (type == "TIF") {
    plot(theta, rowSums(Info), type = "l",
         ylab = "Information",
         xlab = expression(theta[t]),
         main = "Test Information Function", ...)
    lines(theta, 1 / sqrt(rowSums(Info)), lty = 2)
  }
}



