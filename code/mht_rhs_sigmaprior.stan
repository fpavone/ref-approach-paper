data {
  int<lower=0> n; // number of observations 
  vector[n] y; // outputs

  real<lower=0> scale_global; // scale for the half-t prior for tau
  real<lower=1> nu_global; // degrees of freedom for the half-t prior for tau
  
  real<lower=1> nu_local; // degrees of freedom for the half-t priors for lambdas
  real<lower=0> slab_scale; // slab scale for the regularized horseshoe
  
  real<lower=0> slab_df; // slab degrees of freedom for the regularized horseshoe
  // real<lower=0> tau;
  
}

parameters {
  vector[n] z;
  real<lower=0,upper=1> tau; // global shrinkage parameter
  vector<lower=0>[n] lambda; // local shrinkage parameter
  real<lower=0> caux;
  real<lower=0> sigma;  // noise std
}

transformed parameters { 
  vector<lower=0>[n] lambda_tilde; // ’truncated’ local shrinkage parameter 
  real<lower=0> c;  // slab scale
  vector[n] theta; // regression coefficients
  // vector<lower=0,upper=1>[n] k; // shrinkage factor
  // real<lower=0> m_eff; // number of effective parameters different from zero for RHS
  
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)) ); 
  theta = z .* lambda_tilde*tau;
  // for(i in 1:n) {
  //   k[i] = 1/(1 + n*sigma^(-2)*tau^2*lambda[i]^2); 
  // }
  // m_eff = (n - sum(k))*(1-1/(1+n*sigma^(-2)*c)); 
}


model {
// half-t priors for lambdas and tau, and inverse-gamma for c^2 
  sigma ~ lognormal(0,1);
  z ~ normal(0, 1); 
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global*sigma); 
  caux ~ inv_gamma (0.5* slab_df, 0.5*slab_df);
  target += normal_lpdf(y|theta,sigma);
  // for(i in 1:n){
  //   target += normal_lpdf(y[i]|theta[i],sigma);
  // }
}
