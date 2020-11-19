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
  real<lower=0> pr_inno;                // Precision of the innovations
  vector<lower=0>[I] alpha;    // Discrimination parameters
  ordered[K-1] kappa[I];       // Threshold parameters
  vector[nT] inno;             // innovations
  
}

transformed parameters {
  vector[nT] theta;
  real<lower=0> sd_inno;
  real<lower=0> var_inno;
  
  
  theta[1]=inno[1];
  
  for (i in 2:nT) {
     theta[i] = theta[i-1] + inno[i]; 
  }
  
  var_inno=1/pr_inno;
  sd_inno=1/sqrt(pr_inno);

}

model {
  // Priors
  pr_inno ~ gamma(1, 1);
  inno   ~ normal(0, sd_inno);
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




