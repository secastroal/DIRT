functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=0> nT;          // number of time points planned
  int<lower=0> n_knots;     // num of knots
  vector[n_knots] knots;    // the sequence of knots
  int<lower=0> s_degree;    // the degree of spline (is equal to order - 1)
  int<lower=0> I;           // number of items
  int<lower=0> K;           // number of categories per item
  int<lower=1> N;           // number of planned observations
  int<lower=1> N_obs;       // number of observations
  int<lower=1, upper=nT> tt[N]; // time point index planned
  int<lower=1, upper=I>  ii[N]; // item index planned
  int<lower=1, upper=nT> tt_obs[N_obs]; // time point index observed
  int<lower=1, upper=I>  ii_obs[N_obs]; // item index observed
  int y_obs[N_obs];         // vector with observed responses
  real time[nT];            // time of measurements planned
}

transformed data {
  int n_basis = n_knots + s_degree - 1; // total number of B-splines
  matrix[n_basis, nT] B;                // matrix of B-splines
  vector[s_degree + n_knots] ext_knots_temp;
  vector[2*s_degree + n_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], s_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[n_knots], s_degree));
  for (ind in 1:n_basis)
    B[ind,:] = to_row_vector(build_b_spline(time, to_array_1d(ext_knots), ind, s_degree + 1));
  B[n_knots + s_degree - 1, nT] = 1;
}

parameters {
  row_vector[n_basis] a_raw;   // initial intercepts
  real a0;                     // intercept
  real<lower=0> tau;           
  real<lower=-1, upper=1> lambda; // autoregressive effect
  vector[K-1] beta[I];            // Location parameters
  vector[nT] inno;                // innovations
  real<lower=0> sigma2;           // innovations variance
  
}

transformed parameters {
  row_vector[n_basis] a;  // spline coefficients
  vector[nT] tv_int;      // time-varying intercepts
  vector[nT] theta;       // latent state dispositions
  vector[K-1] betasum[I]; // thresholds cumulative sums
  real sigma;             // innovations sd 
  
  // compute time varying intercepts based on the B-splines
  // a = a_raw*tau;
  // for penalized splines use:
  a[1] = a_raw[1];
  for (s in 2:n_basis) {
   a[s] = a[s-1] + a_raw[s] * tau;
  }
  tv_int = a0 * to_vector(time) + to_vector(a * B);
  
  // compute thetas
  // first theta is just the first innovation #!# Discuss with JN and LB
  theta[1] = inno[1];
  // theta[1] = tv_int[1] + inno[1];
  // theta[1] = tv_int[1]/(1 - lambda) + inno[1];
  
  for (t in 2:nT) {
    theta[t] = tv_int[t] + lambda * theta[t - 1] + inno[t];
  }
  
  // thresholds cumulative sums
  for (i in 1:I) {
   betasum[i] = cumulative_sum(beta[i]);
   }
   
  // innovations sd
  sigma = sqrt(sigma2);
}

model {
  // Priors
  a_raw ~ normal(0, 1);
  a0    ~ normal(0, 1);
  tau   ~ normal(0, 1);
  sigma2 ~ normal(1, 1);
  lambda ~ uniform(-1, 1);
  inno   ~ normal(0, sigma);
  
  for(k in 1:I){
    beta[k] ~ normal(0, 1);
  }

  //Likelihood
  for (n in 1:N_obs) {
    vector[K] prob;
    prob[1] = 0;
    for (k in 2:K) {
      prob[k] = (k - 1) * theta[tt_obs[n]] - betasum[ii_obs[n], k - 1]; 
    }   
      y_obs[n] ~ categorical_logit(prob);
    }
}

generated quantities {
   // int rep_y[N_obs];      // posterior simulation
   // vector[N_obs] log_lik; // pointwise loglikelihood
   vector[nT] attractor;  // attractor
   real p_var;            // variance of the process 
   
   // for (n in 1:N_obs) {
   //   vector[K] prob;
   //   prob[1] = 0;
   //   for (k in 2:K) {
   //     prob[k] = (k - 1) * theta[tt_obs[n]] - betasum[ii_obs[n], k - 1];
   //   }
   //   rep_y[n] = categorical_logit_rng(prob);
   //   log_lik[n] = categorical_logit_lpmf(y_obs[n] | prob);
   // }
   
   attractor = tv_int / (1 - lambda);
   p_var     = sigma2 / (1 - square(lambda));

}




