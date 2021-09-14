data {
  int<lower=0> nT;          // number of time points
  int<lower=0> I;           // number of items
  int<lower=0> K;           // number of categories per item
  int<lower=1> N;           // number of planned observations = nT * I
  int<lower=1> N_obs;       // number of observations
  int<lower=1, upper=nT> tt[N]; // time point index planned
  int<lower=1, upper=I>  ii[N]; // item index planned
  int<lower=1, upper=nT> tt_obs[N_obs]; // time point index observed
  int<lower=1, upper=I>  ii_obs[N_obs]; // item index observed
  int y_obs[N_obs];                     // matrix with observed person responses
}

parameters {
  real lambda;                 // autoregressive effect
  vector<lower=0>[I] alpha;    // Discrimination parameters
  ordered[K-1] kappa[I];       // Threshold parameters
  vector[nT] inno;             // innovations
  
}

transformed parameters {
  vector[nT] theta;
  
  theta[1]=inno[1];
  // Here, we estimate a theta value even for 
  // unobserved time points.
  for (i in 2:nT) {
     theta[i] = lambda * theta[i-1] + inno[i]; 
  }  

}

model {
  // Priors
  lambda ~ uniform(-1, 1);
  inno   ~ normal(0, 1);
  alpha ~ lognormal(0, 1);

  for(k in 1:I){
    kappa[k] ~ normal(0, 1);
  }

  //Likelihood
  for (n in 1:N_obs) {
    y_obs[n] ~ ordered_logistic(alpha[ii_obs[n]] * theta[tt_obs[n]], kappa[ii_obs[n]]); 
  }
}

generated quantities {
   ordered[K-1] beta[I]; //item category difficulty
   int rep_y[N_obs];     //posterior simulations 
   int imp_y[N];         //posterior simulations and imputation of NAs 
   vector[N_obs] log_lik;    //pointwise loglikelihood 

 for (i in 1:I){
                beta[i]=kappa[i]/alpha[i];
 }
   
 for (n in 1:N_obs) {
   rep_y[n] = ordered_logistic_rng(alpha[ii_obs[n]] * theta[tt_obs[n]], kappa[ii_obs[n]]);
   log_lik[n] = ordered_logistic_lpmf(y_obs[n] | alpha[ii_obs[n]] * theta[tt_obs[n]], kappa[ii_obs[n]]);
 }
 for (n in 1:N) {
   imp_y[n] = ordered_logistic_rng(alpha[ii[n]] * theta[tt[n]], kappa[ii[n]]);
 }
}




