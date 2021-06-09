data {
  int<lower=0> nT;          // number of time points observed
  int<lower=0> I;           // number of items
  int<lower=0> K;           // number of categories per item
  int<lower=1> N;           // number of observations = nT * I
  int<lower=1, upper=nT> tt[N]; // time point index
  int<lower=1, upper=I>  ii[N]; // item index
  real time[nT];            // time as continuous 
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
  // Here, we include the time difference to 
  // adjust for the missing observations or 
  // uneven spaced over time.
  for (i in 2:nT) {
     theta[i] = lambda * (time[i] - time[i-1]) * theta[i-1] + inno[i]; 
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
      y[n] ~ categorical_logit(prob);
  }
}

generated quantities {
   int rep_y[N];      // posterior simulation
   vector[N] log_lik; // pointwise loglikelihood
   
   for (n in 1:N) {
     vector[K] prob;
     prob[1] = 0;
     for (k in 2:K) {
       prob[k] = (k - 1) * theta[tt[n]] - betasum[ii[n], k - 1];
     }
     rep_y[n] = categorical_logit_rng(prob);
     log_lik[n] = categorical_logit_lpmf(y[n] | prob);
   }
}




