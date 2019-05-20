data{
  int<lower=1> n; //number of observations
  int<lower=1> p; //number of parameters
  real<lower=0> s_max; //standard deviation largest principal component
  matrix[n,p] X; //covariates
  real y[n]; //target
}


parameters{
  real intercept;
  vector[p] beta; //regression weights
  real<lower=0> tau; //variance of regression ceofficient
  real<lower=0> sigma;
}

model{
  sigma ~ student_t(p,0,10);
  tau ~ student_t(p+1,0,s_max^(-2));
  beta ~ normal(0,tau^2);
  intercept ~ normal(0,tau^2);
  
  for(i in 1:n){
    target += normal_lpdf(y[i]|intercept + row(X,i)*beta,sigma);
  }
}

generated quantities {
  vector[n] log_lik;
  for (j in 1:n) {
    log_lik[j] = normal_lpdf(y[j] |intercept + row(X,j)*beta,sigma);
  }
}
