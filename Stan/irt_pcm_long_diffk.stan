data {
  int<lower=0> nT;          // number of time points planned
  int<lower=0> I;           // number of items
  int<lower=0> K[I];        // number of categories per item
  int<lower=1> N;           // number of planned observations = nT * I
  int<lower=1> N_obs;       // number of observations
  int<lower=1, upper=nT> tt[N]; // time point index planned
  int<lower=1, upper=I>  ii[N]; // item index planned
  int<lower=1, upper=nT>  tt_obs[N_obs]; // time point index observed
  int<lower=1, upper=I>  ii_obs[N_obs]; // item index observed
  int y_obs[N_obs];                 // matrix with observed person responses
}

transformed data {
  int M[I];
  int pos_I[I];
  int end_I[I];
  
  for (i in 1:I) M[i] = K[i] - 1; // number of treshold parameters per item
  pos_I[1] = 1;
  end_I[1] = M[1];
  // From DOI 10.1111/rssc.12385
  for (j in 2:I) {
    pos_I[j] = M[j - 1] + pos_I[j - 1];
    end_I[j] = M[j] + end_I[j - 1];
  }
}

parameters {
  vector[sum(M)] beta;         // Location parameters
  real mu_beta;
  real<lower=0> sigma2_beta;
  vector[nT] theta;             // innovations
}

transformed parameters {
  vector[sum(M)] betasum;
  real sigma_beta;

  for (k in 1:I) {
     betasum[pos_I[k]:end_I[k]] = cumulative_sum(beta[pos_I[k]:end_I[k]]);  
  }  
  sigma_beta = sqrt(sigma2_beta);
}

model {
  // Priors
  theta   ~ normal(0, 1);
  mu_beta ~ normal(0, 5);
  sigma2_beta ~ normal(0, 5);

  beta ~ normal(mu_beta, sigma_beta);
  

  // Likelihood
  // Only for the observed data. 
  for (n in 1:N_obs) {
      vector[K[ii_obs[n]]] prob;
      prob[1] = 0;
      for (k in 2:K[ii_obs[n]]) {
        prob[k] = (k - 1) * theta[tt_obs[n]] - betasum[pos_I[ii_obs[n]]:end_I[ii_obs[n]]][k - 1];
      }
      y_obs[n] ~ categorical_logit(prob);
  }
}

generated quantities {
   int rep_y[N_obs];      // posterior simulation
   int imp_y[N];          // posterior simulation and imputation of NAs
   vector[N_obs] log_lik; // pointwise loglikelihood
   
   for (n in 1:N_obs) {
     vector[K[ii_obs[n]]] prob;
     prob[1] = 0;
     for (k in 2:K[ii_obs[n]]) {
       prob[k] = (k - 1) * theta[tt_obs[n]] - betasum[pos_I[ii_obs[n]]:end_I[ii_obs[n]]][k - 1];
     }
     rep_y[n] = categorical_logit_rng(prob);
     log_lik[n] = categorical_logit_lpmf(y_obs[n] | prob);
   }
   
   for (n in 1:N) {
     vector[K[ii_obs[n]]] prob;
     prob[1] = 0;
     for (k in 2:K[ii_obs[n]]) {
       prob[k] = (k - 1) * theta[tt[n]] - betasum[pos_I[ii_obs[n]]:end_I[ii_obs[n]]][k - 1];
     }
     imp_y[n] = categorical_logit_rng(prob);
   }

}




