data {
  int<lower=0> n; // number of observations (indeed features) 
  vector[n] y; // outputs
  // real<lower=0> a; // Dirichlet hyperparameter Dir(a,..,a)
  // vector<lower=0>[n] a_vect;
  real<lower=0> sigma;  // noise std
}

parameters {
  simplex[n] phi; // Dirichlet distributed parameter
  vector<lower=0>[n] lambda; // local shrinkage parameter
  real<lower=0> a;
  
  real<lower=0> tau; // global shrinkage parameter
  vector[n] z; // helper for theta = weights
} 

transformed parameters {
  vector[n] theta; // weights
  vector<lower=0>[n] lambda_s;
  
  lambda_s = sqrt(lambda);
  theta = z .* lambda_s .* phi * tau;
}

model{
  a ~ uniform(1.0/n,0.5);
  tau ~ gamma(n*a,0.5); 
  lambda ~ exponential(0.5);
  // a ~ inv_gamma(2 + 1/(2*n^2),(1+1/(2*n^2))/n);
  phi ~ dirichlet(rep_vector(a,n));
  z ~ normal(0,1); // implies theta_j ~ normal(0, sqrt(lambda_j)*phi_j*tau)

  target += normal_lpdf(y|theta,sigma);
}

