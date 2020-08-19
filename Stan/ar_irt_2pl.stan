data {
  int<lower=0> nT;          // number of time points
  int<lower=0> I;           // number of items
  int<lower=0,upper=1> Y[nT, I]; // matrix with person responses
}

parameters {
  real lambda;                 // autoregressive effect
  real<lower=0> pr_inn;        // precision of the innovation    
  vector<lower=0>[I] alpha;    // Discrimination parameters
  vector[I] beta;              // Difficulty parameters
  real theta_init;             // Theta first value
  vector[nT-1] inno;           // innovations
  
}

transformed parameters {
  real<lower=0> var_inn;
  real<lower=0> sd_inn;
  vector[nT] theta;
  
  var_inn=1/pr_inn;
  sd_inn=1/sqrt(pr_inn);
  theta[1]=theta_init;
  
  for (i in 2:nT) {
     theta[i] = lambda * theta[i-1] + inno[i-1]; 
  }  

}

model {
  // Priors
  lambda ~ normal(0, 1);
  theta_init ~ normal(0, 1);
  pr_inn ~ gamma(1, 1);
  inno   ~ normal(0, sd_inn);
  alpha ~ lognormal(0, 1);
  beta ~ normal(0, 1);

  //Likelihood
  for (i in 1:nT) {
    for (j in 1:I) {
      Y[i, j] ~ bernoulli_logit(alpha[j] * (theta[i] - beta[j]));
    }
  }
}





