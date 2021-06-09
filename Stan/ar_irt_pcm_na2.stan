data {
  int<lower=0> nT;          // number of time points planned
  int<lower=0> I;           // number of items
  int<lower=0> K;           // number of categories per item
  int<lower=1> N;           // number of planned observations = nT * I
  int<lower=1, upper=nT> tt[N]; // time point index planned
  int<lower=1, upper=I>  ii[N]; // item index planned
  int<lower=0, upper=1>  mm[N];// missing values index
  int y[N];                 // matrix with observed person responses
}

parameters {
  real lambda;                 // autoregressive effect
  vector[K-1] beta[I];         // Location parameters
  vector[nT] inno;             // innovations
}

transformed parameters {
  vector[nT] theta;
  vector[K-1] betasum[I];

  theta[1]=inno[1];
  // Here, we estimate a theta value even for 
  // unobserved time points.
  for (i in 2:nT) {
     theta[i] = lambda * theta[i-1] + inno[i]; 
  }
  
  for (k in 1:I) {
     betasum[k] = cumulative_sum(beta[k]);  
  }  

}

model {
  // Priors
  lambda ~ uniform(-1, 1);
  inno   ~ normal(0, 1);
  
  for(k in 1:I){
    beta[k] ~ normal(0, 1);
  }

  // Likelihood
  // Only for the observed data. 
  for (n in 1:N) {
      vector[K] prob;
      prob[1] = 0;
      for (k in 2:K) {
        prob[k] = (k - 1) * theta[tt[n]] - betasum[ii[n], k - 1];
      }
      target += mm[n] * categorical_logit_lpmf(y[n] | prob);
  }
}

generated quantities {
   int rep_y[N];      // posterior simulation
   int imp_y[N];          // posterior simulation and imputation of NAs
   vector[N] log_lik; // pointwise loglikelihood
   
   for (n in 1:N) {
     vector[K] prob;
     prob[1] = 0;
     for (k in 2:K) {
       prob[k] = (k - 1) * theta[tt[n]] - betasum[ii[n], k - 1];
     }
     rep_y[n] = mm[n] * categorical_logit_rng(prob);
     log_lik[n] = mm[n] * categorical_logit_lpmf(y[n] | prob);
   }
   
   for (n in 1:N) {
     vector[K] prob;
     prob[1] = 0;
     for (k in 2:K) {
       prob[k] = (k - 1) * theta[tt[n]] - betasum[ii[n], k - 1];
     }
     imp_y[n] = categorical_logit_rng(prob);
   }

}




