data {
  int<lower=0> nr; // number of rows
  int<lower=0> n;  // number of subjects
  int<lower=0> p;  // number of items
  int<lower=0> K;  // number of categories per items
  int<lower=0> nC; // number of covariates

  int Y[nr];
  real TP[nr];           // time point
  int<lower=0> X1[nr];   // patient indicator
  int<lower=0> item[nr]; // item indicator

  matrix[n,nC] Z;
  int<lower=0> Tobs;
  int<lower=0> indoxp [Tobs];
  
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
  int<lower=0> indo[Tobs];
  int<lower=0> indr[nr];
  
  for (obs in 1:Tobs) indo[obs] = (obs - 1) * p + 1;
  for (r in 1:nr) indr[r] = ((r - 1)/ p) + 1;
}

parameters {
  real mu_gamma1;
  real<lower=0> pr_gamma1;
  real mu_gamma0;
  real<lower=0> pr_gamma0;
  vector[n] gamma0;
  vector[n] gamma1s;
  vector<lower=0>[p] alpha;
  ordered[K-1] kappa[p];
  vector[nC] beta;
}
transformed parameters {
  real<lower=0> sd_gamma1;
  real<lower=0> sd_gamma0;
  vector[Tobs] theta;
  vector[n] gamma1;
  //vector[K] P[nr];    // This is not defined neither used later on
  //vector[K] prob[nr]; // This is not defined neither used later on

  sd_gamma1 = 1/sqrt(pr_gamma1);
  sd_gamma0 = 1/sqrt(pr_gamma0);

  gamma1 = mu_gamma1+sd_gamma1*gamma1s+Z*beta;

  for (obs in 1:Tobs){  // obs: row in thetas
      //theta[obs] = gamma0[nn]+gamma1[nn]*TP[(obs - 1) * p + 1];
    theta[obs] = gamma0[indoxp[obs]]+gamma1[indoxp[obs]]*TP[indo[obs]];
  }
}
model {
  mu_gamma1 ~ normal(m_mu_gamma1, sd_mu_gamma1);
  pr_gamma1 ~ gamma(a_pr_gamma1, b_pr_gamma1);
  mu_gamma0 ~ normal(0, 1);
  pr_gamma0 ~ gamma(1, 1);
  gamma0 ~ normal(mu_gamma0, sd_gamma0);
  gamma1s ~ normal(0,1); // standardized gamma 1
  beta ~ normal(0,1);

  alpha ~ normal(m_alpha, sd_alpha);

  for(kk in 1:p){
    kappa[kk] ~ normal(m_kappa, sd_kappa);
  }

  for(r in 1:nr){
    //Y[r] ~ ordered_logistic(alpha[item[r]]*theta[((r - 1)/ p) + 1],kappa[item[r]]);
    Y[r] ~ ordered_logistic(alpha[item[r]]*theta[indr[r]],kappa[item[r]]);
  }
}

