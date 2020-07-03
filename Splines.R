# Splines tutorial: https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html

library("splines")
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# Splines is a continuous piece-wise polynomial regression. In short, we have a variable
# X that we use to predict Y, and b-splines are used to obtain the best fit to the data.
# The splines are a group of polynomial functions, which are defined for short intervals 
# of the variable X. The most common splines used are splines of degree 4. This means that
# the polynomial functions used to predict Y have and order 3 or lower. A function of each
# order is defined for each segment of X and the sum of these functions result in the 
# spline.

# In this script, I follow step by step the tutorial of Splines in Stan by Milad Kharratzadeh.
# The tutorial includes three spline stan models and shows how to generate data and fit
# the model. The first example just show how to generate the data in R and the stan code 
# is a simple linear regression because the b-splines are required as input for the stan
# model. The second example allows to compute the splines within the stan model. Finally, 
# the third example shows a penalized splines models that smooth the function in case too
# many knots are required.

# Example 1 ----

# Generating data and splines in R. The stan model requires X, Y, and the b-splines
# as input, hence, it assumes the b-splines and the number of knots are known. 

# Generate X (independent variable). Just a sequence between -5 and 5.
X <- seq(from = -5, to = 5, by = 0.1)
# Compute b-splines within the interval (-5, 5).
# Here, we generate b-splines of order fourth, which are cubic polynomial functions. These 
# splines are built based on 11 knots, which divides the variable X in 10 intervals. 
B <- t(bs(X, knots = seq(-5, 5, 1), degree = 3, intercept = TRUE)) # creating the B-splines

num_data  <- length(X)
num_basis <- nrow(B)

# Define intercept and coefficients to generate Y based on X and the b-splines.
a0 <- 0.2                    # intercept
a  <- rnorm(num_basis, 0, 1) # coefficients of B-splines

# Compute predicted values of Y and add random error to generate Y. 
Y_true <- as.vector(a0 * X + a %*% B)    # generating the output
Y      <- Y_true + rnorm(length(X),0,.2) # adding noise

# Now, we fit the model in stan. The file is "splines_1.stan". First, we create 
# the data for stan in a list, then we fit the model.
standata <- list(num_data  = num_data, 
                 num_basis = num_basis, 
                 X = X, 
                 Y = Y, 
                 B = B)

fit <- stan(file = "Stan/splines_1.stan",     # Stan model. 
            data = standata,                  # Data.
            iter = 1200,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 200,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

# Check traceplots of the predicted Y given the b-splines model.
traceplot(fit, pars = c("Y_hat[1]",  "Y_hat[13]", "Y_hat[25]",
                        "Y_hat[38]", "Y_hat[50]", "Y_hat[62]",
                        "Y_hat[75]", "Y_hat[87]", "Y_hat[100]"), 
          inc_warmup = FALSE)
# Extract predicted Y.
Y_stan <- summary(fit, pars = "Y_hat")$summary

# Scatter plot of X and Y. Including the true predicted Y values in blue and the estimated
# predicted values in red (With 95% credibility interval).
plot(X, Y, pch = 20)
polygon(c(X, rev(X)),
        c(Y_stan[, 4], rev(Y_stan[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
lines(X, Y_stan[, 1], col = "red", lwd = 2)
lines(X, Y_true, col = "blue", lwd = 2)

rm(list = ls())

# Example 2 ----

# In this example, the b-splines are computed in the stan model. Therefore, only the
# knots, the order of the splines - 1, and the data are needed to run the model in stan.
# For this analysis, we use the stan model "splines_2.stan".

# Generate X (independent variable). Just a sequence between -5 and 5.
X <- seq(from = -5, to = 5, by = 0.1)

num_knots     <- 11 # Define number of knots.
spline_degree <- 3  # Define degree of splines.
num_basis     <- num_knots + spline_degree - 1 # Compute number of basis.

num_data  <- length(X)
knots <- unname(quantile(X, probs = seq(0, 1, length.out = num_knots)))

# Compute b-splines within the interval (-5, 5).
# Here, we generate b-splines of order fourth, which are cubic polynomial functions. These 
# splines are built based on num_knots, num_basis, and spline_degree.
B_true <- t(bs(X, df = num_basis, degree = spline_degree, intercept = TRUE)) # creating the B-splines

# Define intercept and coefficients to generate Y based on X and the b-splines.
a0 <- 0.2                    # intercept
a  <- rnorm(num_basis, 0, 1) # coefficients of B-splines

# Compute predicted values of Y and add random error to generate Y. 
Y_true <- as.vector(a0 * X + a %*% B_true) # generating the output
Y      <- Y_true + rnorm(length(X),0,.2)   # adding noise

# Now, we fit the model in stan. The file is "splines_2.stan". First, we create 
# the data for stan in a list, then we fit the model.
standata <- list(num_data      = num_data, 
                 num_knots     = num_knots,
                 knots         = knots,
                 spline_degree = spline_degree,
                 X = X, 
                 Y = Y)

fit <- stan(file = "Stan/splines_2.stan",     # Stan model. 
            data = standata,                  # Data.
            iter = 1200,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 200,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

# Check traceplots of the predicted Y given the b-splines model.
traceplot(fit, pars = c("Y_hat[1]",  "Y_hat[13]", "Y_hat[25]",
                        "Y_hat[38]", "Y_hat[50]", "Y_hat[62]",
                        "Y_hat[75]", "Y_hat[87]", "Y_hat[100]"), 
          inc_warmup = FALSE)
# Extract predicted Y.
Y_stan <- summary(fit, pars = "Y_hat")$summary

# Scatter plot of X and Y. Including the true predicted Y values in blue and the estimated
# predicted values in red (With 95% credibility interval).
plot(X, Y, pch = 20)
polygon(c(X, rev(X)),
        c(Y_stan[, 4], rev(Y_stan[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
lines(X, Y_stan[, 1], col = "red", lwd = 2)
lines(X, Y_true, col = "blue", lwd = 2)

rm(list = ls())

# Example 3 ----

# The third example considers penalized splines. Usually, it is unclear how many knots are
# needed to best model the data. If too few/many knots are selected, one can underfit/overfit
# the data. To solve this, one can define some priors that enforce smoothness. The stan
# model for is in "splines_3.stan".

# Here, the idea is to generate data based on a few number of knots but to analyze it based 
# on a much larger number of knots and compare the results from the normal splines model
# with the penalized splines model.

# Basically, we generate the data as in Example 2.
# Generate X (independent variable). Just a sequence between -5 and 5.
X <- seq(from = -10, to = 10, by = 0.1)

num_knots     <- 11 # Define number of knots.
spline_degree <- 3  # Define degree of splines.
num_basis     <- num_knots + spline_degree - 1 # Compute number of basis.

num_data  <- length(X)
knots <- unname(quantile(X, probs = seq(0, 1, length.out = num_knots)))

# Compute b-splines within the interval (-5, 5).
# Here, we generate b-splines of order fourth, which are cubic polynomial functions. These 
# splines are built based on num_knots, num_basis, and spline_degree.
B_true <- t(bs(X, df = num_basis, degree = spline_degree, intercept = TRUE)) # creating the B-splines

# Define intercept and coefficients to generate Y based on X and the b-splines.
a0 <- 0.2                    # intercept
a  <- rnorm(num_basis, 0, 1) # coefficients of B-splines

# Compute predicted values of Y and add random error to generate Y. 
Y_true <- as.vector(a0 * X + a %*% B_true) # generating the output
Y      <- Y_true + rnorm(length(X),0,.2)   # adding noise

# Next, we overwrite the number of knots and basis to overfit the data in purpose. This
# will allow us to compare the simple splines model and the penalized splines model. 

num_knots <- 100
num_basis <- num_knots + spline_degree - 1
knots <- unname(quantile(X, probs = seq(0, 1, length.out = num_knots)))

# Now, we fit the splines and the penalized splines models in stan.
standata <- list(num_data      = num_data, 
                 num_knots     = num_knots,
                 knots         = knots,
                 spline_degree = spline_degree,
                 X = X, 
                 Y = Y)

fit <- stan(file = "Stan/splines_2.stan",     # Stan model. 
            data = standata,                  # Data.
            iter = 1200,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 200,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

fit_p <- stan(file = "Stan/splines_3.stan",     # Stan model. 
              data = standata,                  # Data.
              iter = 1200,                      # Number of iterations.
              chains  = 3,                      # Number of chains.
              warmup  = 200,                    # Burn-in samples.
              control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

# Check traceplots of the predicted Y given the b-splines model.
traceplot(fit, pars = c("Y_hat[1]",  "Y_hat[13]", "Y_hat[25]",
                        "Y_hat[38]", "Y_hat[50]", "Y_hat[62]",
                        "Y_hat[75]", "Y_hat[87]", "Y_hat[100]"), 
          inc_warmup = FALSE)

traceplot(fit_p, pars = c("Y_hat[1]",  "Y_hat[13]", "Y_hat[25]",
                          "Y_hat[38]", "Y_hat[50]", "Y_hat[62]",
                          "Y_hat[75]", "Y_hat[87]", "Y_hat[100]"), 
          inc_warmup = FALSE)

# Extract predicted Y.
Y_stan   <- summary(fit, pars = "Y_hat")$summary
Y_stan_p <- summary(fit_p, pars = "Y_hat")$summary


# Scatter plot of X and Y. Including the true predicted Y values in blue and the estimated
# predicted values in red (With 95% credibility interval).
plot(X, Y, pch = 20)
#polygon(c(X, rev(X)),
#        c(Y_stan[, 4], rev(Y_stan[, 8])),
#        border = NA,
#        col = rgb(1, 0, 0, 0.25))
#polygon(c(X, rev(X)),
#        c(Y_stan_p[, 4], rev(Y_stan_p[, 8])),
#        border = NA,
#        col = rgb(0, 1, 0, 0.25))
lines(X, Y_stan[, 1], col = "red", lwd = 2)
lines(X, Y_stan_p[, 1], col = "blue", lwd = 2)
lines(X, Y_true, col = "darkgray", lwd = 2, lty = 2)
legend("bottomright", c("True Y_hat", "Est Y_hat", "Est Y_hat + Smoothing"), 
       col = c("darkgray", "red", "blue"),
       lty = c(2, 1, 1),
       cex = 0.75,
       box.lty = 0,
       border = "white")

rm(list = ls())

# IRT-splines ----
source("R/IRT_models.R")

# In this section, I try to make an IRT-splines model. The idea is to assume that
# we have the time series of a test of a single individual. The relation between the 
# items' responses and the individual are modeled through the GRM. Moreover, the relation 
# between theta and time are modeled through a b-splines regression.

# First, let's assume the person was measured on 200 different occasions. Hence, our 
# dependent variable X is actually the time at which the person was measured.

time <- 1:200

# If the person was measured 10 times a day during 20 days, it makes sense to define 
# knots to separate measurements between days. This means 21 knots.
# We keep using cube splines.
num_knots     <- 21 # Define number of knots.
spline_degree <- 3  # Define degree of splines.
num_basis     <- num_knots + spline_degree - 1 # Compute number of basis.

num_data  <- length(time)
knots <- unname(quantile(time, probs = seq(0, 1, length.out = num_knots)))

# Compute b-splines within the time interval.
# Here, we generate b-splines of order fourth, which are cubic polynomial functions. These 
# splines are built based on num_knots, num_basis, and spline_degree.
B_true <- t(bs(time, df = num_basis, degree = spline_degree, intercept = TRUE)) # creating the B-splines

# Define intercept and coefficients to generate theta based on time and the b-splines.
a0 <- 0                      # intercept
a  <- rnorm(num_basis, 0, 1) # coefficients of B-splines
#!# How do we interpret this intercept? is it the trait theta?

# Compute theta based on time and b-splines. 
theta <- as.vector(a0 * time + a %*% B_true) # generating theta
#theta <- theta + rnorm(length(X),0,.2) 
#!# Should I add some noise on the relation between time and theta?

# Next, we generate data based on the GRM and the thetas we just created.

I <- 10 # Number of items.
K <- 5  # Number of categories per items.

# Create item parameters
alpha      <- rlnorm(I, 0, 0.25)    # Discrimination parameters.
delta      <- matrix(NA, I, K - 1)  # Matrix to store difficulty parameters.
delta[, 1] <- rnorm(I, -2)
for (i in 2:(K - 1)){
        delta[, i] <- delta[, i - 1] + runif(I, 0.4, 1.2)  #!#!# edited!
}
delta <- delta - rowMeans(delta)                     #!#!# edited!
delta <- delta + rnorm(I)                            #!#!# edited!
rm(i)

IP <- cbind(delta, alpha) # Matrix of item parameters to input in P.GRM

probs.array <- P.GRM(K - 1, IP, theta)
probs.array[probs.array < 0] <- 0
responses   <- apply(probs.array, 1:2, function(vec) {which( rmultinom(1, 1, vec) == 1) - 1 })
responses   <- responses + 1 # To fit the GRM in stan, items should be coded starting from 1. 
rm(probs.array, delta, alpha)

# Fitting an IRT-splines model in stan. This first try intends to combine the splines model 
# in example 2 with the GRM model. 

standata <- list(num_data      = num_data, 
                 num_knots     = num_knots,
                 knots         = knots,
                 spline_degree = spline_degree,
                 I = I,
                 K = K,
                 time = time, 
                 Y = responses)

fit <- stan(file = "Stan/splines_irt.stan",   # Stan model. 
            data = standata,                  # Data.
            iter = 1200,                      # Number of iterations.
            chains  = 3,                      # Number of chains.
            warmup  = 200,                    # Burn-in samples.
            control = list(adapt_delta=0.95)) # Other parameters to control sampling behavior.

traceplot(fit, pars = "alpha", 
          inc_warmup = FALSE)

for (i in 1:I) {
        print(traceplot(fit, pars = paste0("kappa[", i, ",", 1:4, "]"), 
                        inc_warmup = FALSE))
}

for (i in 1:I) {
        print(traceplot(fit, pars = paste0("beta[", i, ",", 1:4, "]"), 
                  inc_warmup = FALSE))
}

traceplot(fit, pars = paste0("theta[", c(1, 13, 25, 38, 50, 62, 75, 87, 100), "]"),
          inc_warmup = FALSE)

# Extract theta.
theta_est <- summary(fit, pars = "theta")$summary

# Scatter plot of time and theta. Including the estimated theta in red (With 95% credibility interval).
plot(time, theta, pch = 20)
polygon(c(time, rev(time)),
        c(theta_est[, 4], rev(theta_est[, 8])),
        border = NA,
        col = rgb(1, 0, 0, 0.25))
lines(time, theta_est[, 1], col = "red", lwd = 2)
