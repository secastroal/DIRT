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
  real<lower=0> pr_inno;      // Precision of the innovations
  vector[K-1] beta[I];        // Location parameters
  vector[nT] inno;            // innovations
  
}

transformed parameters {
  vector[nT] theta;
  vector[K-1] betasum[I];
  real<lower=0> sd_inno;
  real<lower=0> var_inno;
  
  theta[1]=inno[1];
  
  for (i in 2:nT) {
     theta[i] = theta[i-1] + inno[i]; 
  }
  
  for (k in 1:I) {
     betasum[k] = cumulative_sum(beta[k]);  
  }
  
  var_inno=1/pr_inno;
  sd_inno=1/sqrt(pr_inno);

}

model {
  // Priors
  pr_inno ~ gamma(1, 1);
  inno   ~ normal(0, sd_inno);
  
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




