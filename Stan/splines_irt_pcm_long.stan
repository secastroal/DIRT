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
  int num_data;             // number of data points
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  int<lower=0> I;           // number of items
  int<lower=0> K;           // number of categories per item
  int<lower=1> N;           // number of observations
  int<lower=1, upper=num_data> tt[N]; // time point index
  int<lower=1, upper=I>  ii[N];       // item index
  int y[N];                 // matrix with person responses
  real time[num_data];      // time of measurements
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_data] B;                 // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(time, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[num_knots + spline_degree - 1, num_data] = 1;
}

parameters {
  row_vector[num_basis] a_raw; 
  real a0;                     // intercept
  real<lower=0> tau;
  vector[K-1] beta[I];       // Location parameters
  
}

transformed parameters {
  row_vector[num_basis] a; // spline coefficients
  vector[num_data] theta;
  vector[K-1] betasum[I];
  a = a_raw*tau;
  theta = a0 * to_vector(time) + to_vector(a * B);
  for (i in 1:I) {
   betasum[i] = cumulative_sum(beta[i]);
   }
}

model {
  // Priors
  a_raw ~ normal(0, 1);
  a0    ~ normal(0, 1);
  tau   ~ normal(0, 1);
  
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
   int y_rep[N];      // posterior simulation
   vector[N] log_lik; // pointwise loglikelihood
   
   for (n in 1:N) {
     vector[K] prob;
     prob[1] = 0;
     for (k in 2:K) {
       prob[k] = (k - 1) * theta[tt[n]] - betasum[ii[n], k - 1];
     }
     y_rep[n] = categorical_logit_rng(prob);      
     log_lik[n] = categorical_logit_lpmf(y[n] | prob);
   }

}




