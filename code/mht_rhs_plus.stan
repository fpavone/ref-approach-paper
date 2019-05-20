data {
  int<lower=0> n; // number of observations 
  vector[n] y; // outputs
  real<lower=0> sigma;
}

// parameters {
//   vector[n] z;
//   real<lower=0,upper=1> tau; // global shrinkage parameter
//   vector<lower=0>[n] lambda2; // local shrinkage parameter
//   vector<lower=0>[n] nu; 
//   vector<lower=0>[n] phi2; // local shrinkage parameter
//   vector<lower=0>[n] eta;
// }
// 
// transformed parameters{
//   vector[n] theta;
// 
//   theta = z .* sqrt(lambda2);
// }
// 
// model {
//   lambda2 ~ inv_gamma(0.5,1 ./ nu);
//   nu ~ inv_gamma(0.5,1 ./ phi2);
//   phi2 ~ inv_gamma(0.5,1 ./ eta);
//   eta ~ inv_gamma(0.5,1);
//   z ~ normal(0,1);
//   
//   target += normal_lpdf(y|theta,sigma);
// }

parameters {
  vector[n] z;
  real<lower=0,upper=1> tau; // global shrinkage parameter
  vector<lower=0>[n] lambda; // local shrinkage parameter
  vector<lower=0>[n] eta;
}

transformed parameters {
  vector[n] theta;
  theta = z .* lambda;
}

model {
  eta ~ cauchy(0,1);
  tau ~ cauchy(0,1.0/n);
  z ~ normal(0,1);
  lambda ~ cauchy(0,eta*tau);

  target += normal_lpdf(y|theta,sigma);
}
