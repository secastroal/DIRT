data {
  int<lower=0> nT;          // number of time points
  int<lower=0> I;           // number of items
  int<lower=0> K;           // number of categories per item
  int<lower=1> N;           // number of observations = nT * I
  int<lower=1, upper=nT> tt[N]; // time point index
  int<lower=1, upper=I>  ii[N]; // item index
  int y[N];                     // matrix with person responses
}

parameters {
  real lambda;                 // autoregressive effect
  vector[K-1] beta[I];         // Location parameters
  vector[nT] inno;           // innovations
  
}

transformed parameters {
  vector[nT] theta;
  vector[K-1] betasum[I];

  theta[1]=inno[1];
  
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

  //Likelihood
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




