data {
  int<lower=0> nT;          // number of time points
  int<lower=0> I;           // number of items
  int<lower=1> N;               // number of observations = nT*I
  int<lower=1, upper=nT> tt[N]; // time point index
  int<lower=1, upper=I>  ii[N]; // item index
  int<lower=0,upper=1>    y[N]; // matrix with person responses
}

parameters {
  real lambda;                 // autoregressive effect
  //real<lower=0> pr_inn;        // precision of the innovation    
  vector<lower=0>[I] alpha;    // Discrimination parameters
  vector[I] beta;              // Difficulty parameters
  real theta_init;             // Theta first value
  vector[nT-1] inno;           // innovations
  
}

transformed parameters {
  //real<lower=0> var_inn;
  //real<lower=0> sd_inn;
  vector[nT] theta;
  
  //var_inn=1/pr_inn;
  //sd_inn=1/sqrt(pr_inn);
  theta[1] = theta_init;
  
  for (i in 2:nT) {
     theta[i] = lambda * theta[i-1] + inno[i-1]; 
  }  

}

model {
  vector[N] eta;
  
  // Priors
  lambda     ~ uniform(-1, 1);
  theta_init ~ normal(0, 1);
  //pr_inn ~ gamma(1, 1);
  inno       ~ normal(0, 1);
  alpha      ~ lognormal(0, 1);
  beta       ~ normal(0, 1);

  //Likelihood
  for (n in 1:N) {
    eta[n] = alpha[ii[n]] * (theta[tt[n]] - beta[ii[n]]);
  }
  y ~ bernoulli_logit(eta);
    
}





