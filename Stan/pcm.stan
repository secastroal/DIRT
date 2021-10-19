data{
    int<lower=2, upper=7> K; //  number of categoies
    int <lower=0> n_student; //  number of individuals
    int <lower=0> n_item; //  number of items
    int<lower=1,upper=K> Y[n_student,n_item]; //array of responses
 }

 parameters {
    vector[K-1] beta[n_item]; //item difficulty parameter
    real mu_beta; // mean difficulty
    real<lower=0> sigma_beta; //difficulty sd
    vector[n_student] theta; // latent trait
 }

transformed parameters{
  vector[K-1] betasum[n_item];
  
  for (i in 1:n_item) {
    betasum[i] = cumulative_sum(beta[i]); 
  }
}

 model{
    theta ~ normal(0,1);
  for (i in 1:n_item){
    beta[i] ~ normal(mu_beta,sigma_beta);
                     }
    mu_beta ~ normal(0,5);
    sigma_beta ~ cauchy(0,5);

    for (i in 1:n_student){
       for (j in 1:n_item){
         vector[K] prob;
         prob[1] = 0;
         for (k in 2:K) {
           prob[k] = (k - 1) * theta[i] - betasum[j, k - 1]; 
         } 
         Y[i,j] ~ categorical_logit(prob);
 }}
 }


