data{
    int<lower=2, upper=5> K; // number of categories
    int <lower=0> n_student; //  number of individuals
    int <lower=0> n_item; //  number of items
    int<lower=1,upper=K> Y[n_student,n_item]; //array of responses
 }

 parameters {
    vector[n_student] theta; // latent trait
    real<lower=0> alpha [n_item]; // item discrimination
    ordered[K-1] kappa[n_item]; //item category intercept
    real theta_mean; //hyper mean for theta distribution
    real<lower=0> theta_var; //hyper variance for theta distribution 
 }

 model{
    alpha ~ lognormal(0,1);
    theta_mean ~ normal(0, 1);
    theta_var ~ gamma(1, 1);
    theta ~ normal(theta_mean, theta_var);

    for (i in 1: n_item){
       for (k in 1:(K-1)){ 
    kappa[i,k] ~ normal(0,1); //The prior of kappa is informative to be able to scale and identify the model
                        }}

    for (i in 1:n_student){
          for (j in 1:n_item){
            Y[i,j] ~ ordered_logistic(theta[i]*alpha[j],kappa[j]);
 }}
 }

 generated quantities {
   ordered[K-1] beta[n_item]; //item category difficulty

 for (i in 1:n_item){
                beta[i]=kappa[i]/alpha[i];
 }
 }


