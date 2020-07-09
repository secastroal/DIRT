functions {
  int[] seq(int x){
    int res[x];
    for (i in 1:x){
      res[i] = i;
    }
    return res;
  }
}

data {
  int<lower=0> n;  // number of subjects
  int<lower=0> nt; // number of time points
  int<lower=0> p;  // number of items
  int<lower=0> K;  // number of categories per item
  int<lower=0> nC; // number of covariates

  int Y[n, nt, p];
  matrix[n, nC] Z;

  real m_mu_gamma1;
  real<lower=0> sd_mu_gamma1;
  real m_alpha;
  real<lower=0> sd_alpha;
  real m_kappa;
  real<lower=0> sd_kappa;
  real<lower=0> a_pr_gamma1;
  real<lower=0> b_pr_gamma1;
}

transformed data {
  // int<lower=0> X1[n];   // patient indicator
  int<lower=0> TP[nt];  // time point
  // int<lower=0> item[p]; // item indicator
  
  // X1   = seq(n);
  TP   = seq(nt);
  // item = seq(p);
}

parameters {
  real mu_gamma1;
  real<lower=0> pr_gamma1;
  vector[n] gamma0;
  vector[n] gamma1s;
  vector<lower=0>[p] alpha;
  ordered[K-1] beta_i[p];
  vector[nC] beta;
}

transformed parameters {
  real<lower=0> sd_gamma1;
  real theta[n, nt];
  vector[n] gamma1;

  sd_gamma1 = 1/sqrt(pr_gamma1);
  gamma1 = mu_gamma1 + sd_gamma1*gamma1s + Z*beta;

  // for (r in 1:nr){  // r: row in dataset, k: value category
  //   theta[r] = gamma0[X1[r]]+gamma1[X1[r]]*TP[r];
  // }
  // Changed, see below:
  for (i in 1:n){
    for (t in 1:nt) {
      theta[i, t] = gamma0[i] + gamma1[i] * TP[t];
    }
  }
}

model {
  mu_gamma1 ~ normal(m_mu_gamma1, sd_mu_gamma1);
  pr_gamma1 ~ gamma(a_pr_gamma1, b_pr_gamma1);
  gamma0    ~ normal(0,1);
  gamma1s   ~ normal(0,1); // standardized gamma 1
  beta      ~ normal(0,1);
  alpha     ~ normal(m_alpha, sd_alpha);

  for(k in 1:p){
    beta_i[k] ~ normal(m_kappa, sd_kappa);
  }

  // for(r in 1:nr){
  //   Y[r] ~ ordered_logistic(alpha[item[r]]*theta[r],alpha[item[r]]*kappa[item[r]]);
  // }
  for (i in 1:n) {
    for (t in 1:nt) {
      for (k in 1:p) {
        Y[i, t, k] ~ ordered_logistic(alpha[k]*theta[i, t], alpha[k]*beta_i[k]);
        // Read stan manual on this function and cross we data generation.
      }
    }
  }
}

generated quantities {
   ordered[K-1] kappa[p]; //item category difficulty

 for (i in 1:p){
                kappa[i]=beta_i[i]*alpha[i];
 }
}



