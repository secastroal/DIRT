data {
  int<lower=0> nr; // number of rows
  int<lower=0> n;  // number of subjects
  int<lower=0> p;  // number of items
  int<lower=0> K;  // number of categories per items
  
  int Y[nr];
  real TP[nr];           // time point
  int<lower=0> X1[nr];   // patient indicator
  int<lower=0> item[nr]; // item indicator

}
parameters {
  // Parameters
  vector[n] mu_i;          // Random intercept
  //vector[n] beta_i;        // Random slope for TP
  vector[n] r_i;           // Random amplitude
  vector[n] phi_i;         // Random phase shift
  real<lower = 0> sigma_2; // Residual variance
  vector[p - 1] lambda_f;  // Location parameters free
  vector[K - 2] delta_f;   // Tresholds parameters free
  //vector[p] lambda;        // Location parameters  
  //vector[K - 1] delta;     // Threshold parameters
  vector[nr] theta;        // Theta of person i at time t
  
  // Hyper parameters
  real mu_gamma0;         // Mean mu
  real<lower=0> v_gamma0; // Var mu
  //real mu_gamma1;         // Mean beta  
  //real<lower=0> v_gamma1; // Var beta
  real mu_gamma2;         // Mean r
  real<lower=0> v_gamma2; // Var r
  real mu_gamma3;         // Mean phi
  real<lower=0> v_gamma3; // Var phi
}
transformed parameters {
  real<lower=0> s_gamma0; // Sd mu
  //real<lower=0> s_gamma1; // Sd beta
  real<lower=0> s_gamma2; // Sd r
  real<lower=0> s_gamma3; // Sd phi
  real<lower=0> sigma;    // Sd residual variance
  vector[p] lambda;       // Location parameters
  vector[K - 1] delta;    // Threshold parameters
  vector[K] deltasum; // Cumulative sum threshold parameters
  vector[nr] theta_i;     // Predicted Theta of person i at time t
  vector[K] args [nr];        // Create vector argument for softmax function
  vector[K] probs [nr];   // Output vector of the softmax function 
  
  // Compute standard deviations
  s_gamma0 = sqrt(v_gamma0);
  //s_gamma1 = sqrt(v_gamma1);
  s_gamma2 = sqrt(v_gamma2);
  s_gamma3 = sqrt(v_gamma3);
  sigma    = sqrt(sigma_2);
  
  // Constrain sum of location and threshold parameters to 0
  for (i in 1:(p-1)) {
    lambda[i] = lambda_f[i]; 
  }
  lambda[p] = -sum(lambda_f);
  
  for (j in 1:(K-2)) {
    delta[j] = delta_f[j]; 
  }
  delta[K-1] = -sum(delta_f);
  
  // Compute cumulative sum of threshold parameters
  deltasum[1] = 0.0;
  deltasum[2:K] = cumulative_sum(delta);
  
  // Compute intraindividual means of theta
  
  for (rr in 1:nr) {
    theta_i[rr] = mu_i[X1[rr]] + r_i[X1[rr]] * (cos(((2*pi())/5))*TP[rr] + phi_i[X1[rr]]); //+ beta_i[X1[rr]] * TP[rr];
  }
  
  // Compute probabilities to input in likelihood
  for (s in 1:nr) {
    for (y in 0:(K - 1)) args[s, y + 1] = (y * (theta[s] - lambda[item[s]]) + deltasum[y + 1]);
    probs[s] = softmax(args[s]);
  }
}
model {
  // Define priors
  mu_i ~ normal(mu_gamma0, s_gamma0);          
  //beta_i ~ normal(mu_gamma1, s_gamma1);
  r_i ~ normal(mu_gamma2, s_gamma2);           
  phi_i ~ normal(mu_gamma3, s_gamma3);
  sigma_2 ~ inv_gamma(1, 1) ;// Residual variance
  lambda_f ~ normal(0, 1);  // Location parameters free
  delta_f ~ normal(0, 1);   // Tresholds parameters free
  //lambda ~ normal(0, 1);
  //delta ~ normal(0, 1);
  theta ~ normal(theta_i, sigma);        // Theta of person i at time t
  
  // Hyper parameters
  
  // Define Hyperpriors
  mu_gamma0 ~ normal(0, 1);
  v_gamma0 ~ inv_gamma(1, 1);
  //mu_gamma1 ~ normal(0, 1);   
  //v_gamma1 ~ inv_gamma(1, 1);
  mu_gamma2 ~ normal(0, 1);  
  v_gamma2 ~ inv_gamma(1, 1);
  mu_gamma3 ~ normal(0, 1);
  v_gamma3 ~ inv_gamma(1, 1); 
  
  // Likelihood
  for(rr in 1:nr){
    Y[rr] ~ categorical(probs[rr]);
  }
}

