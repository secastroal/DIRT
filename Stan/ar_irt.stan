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
  vector<lower=0>[I] alpha;    // Discrimination parameters
  ordered[K-1] kappa[I];       // Threshold parameters
  vector[nT] inno;             // innovations
  
}

transformed parameters {
  vector[nT] theta;
  
  theta[1]=inno[1];
  
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
  for (n in 1:N) {
    y[n] ~ ordered_logistic(alpha[ii[n]] * theta[tt[n]], kappa[ii[n]]); 
  }
}

generated quantities {
   ordered[K-1] beta[I]; //item category difficulty
   int rep_y[N];         //posterior simulations 
   vector[N] log_lik;    //pointwise loglikelihood 

 for (i in 1:I){
                beta[i]=kappa[i]/alpha[i];
 }
   
 for (n in 1:N) {
   rep_y[n] = ordered_logistic_rng(alpha[ii[n]] * theta[tt[n]], kappa[ii[n]]);
   log_lik[n] = ordered_logistic_lpmf(y[n] | alpha[ii[n]] * theta[tt[n]], kappa[ii[n]]);
 }   
}




