library(tidyverse)
library(dimreduce)
library(locfdr)
library(EbayesThresh)
library(rstan)
library(rstanarm)
options(mc.cores = parallel::detectCores())
set.seed(453876)

args = commandArgs(trailingOnly=TRUE) # n rho nc

## Simulation parameters
n <- as.numeric(args[1])  # number of observations
rho <- as.numeric(args[2])  # correlation level
p <- 1000   # total number of features
k <- 100    # number of relevant features

times <- 100

## Data simulation mechanism
simulate_data <- function(){ # n, rho, p and k from global environment
  f <- rnorm(n) # the true underlying function
  
  data <- tibble(y = f + rnorm(n)) %>%   # observed targe variable y
    bind_cols(as_tibble(sqrt(1-rho)*matrix(rnorm(n*k), ncol=k) + sqrt(rho)*f)) %>%  # set of relevant covariates
    bind_cols(as_tibble(matrix(rnorm(n*(p-k)), ncol=p-k))) %>%  # set of spurious covariates
    set_names(c("y",paste("r.",1:k,sep=""),paste("s.",(k+1):p,sep="")))
  
  return(data)
}

results <- tibble(error = numeric(),
                  error_type = character(),
                  prior = character(),
                  approach = character(),
                  n = numeric(),
                  rho = numeric())

## Stan models
fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))
rhs_model <- stan_model("mht_rhs.stan", model_name = "rhs_model",
                        auto_write = rstan_options("auto_write"=TRUE))
dl_model <- stan_model("dl_shrinkage.stan", model_name = "dl_model",
                        auto_write = rstan_options("auto_write"=TRUE))

k_RHS_ref <- numeric(p)
k_RHS_data <- numeric(p)

for(i in 1:times){
  data <- simulate_data() 
  y <- dplyr::pull(data, var="y")
  x <- dplyr::select(data, -contains("y"))
  
  
  ## Fitting the reference model
  nc <- 5
  dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F)
  z <- predict(dr, x)
  
  model_data <- list(y = y,
                     X = z,
                     p = nc,
                     n = n,
                     s_max = dr$sdev[1])
  fit <- sampling(fit_model, data = model_data, 
                  chains = 1, iter = 2000, seed = 45342)
  
  draws <- as.matrix(fit) # posterior draws
  sigma <- draws[,'sigma'] # noise std
  beta <- draws[,2:(nc+1)] # regression coefficients
  alpha <- draws[,'intercept'] # intercept
  
  predfun <- function(zt){
    colMns <- colMeans(x)
    colSds <- apply(x,2,sd)
    zt <- t(apply(zt,1,function(x){(x-colMns)/colSds}))%*%dr$w
    return(t( beta %*% t(zt) + alpha ))   # n x s
  }
  
  mufit <- predfun(x)
  yfit <- apply(mufit,1,function(x){ # predict with the posterior predictive mean
    mean(rnorm(1000,x,sigma))
  })
  
  ## Computing sample correlation and z-values using reference model's predictions and observed data
  r_ref <- as.vector(cor(yfit,x))
  z_ref <- sqrt(n-3) * 0.5*log((1+r_ref)/(1-r_ref)) # Fisher transformation * (n-3)
  
  r_data <- as.vector(cor(y,x))
  z_data <- sqrt(n-3) * 0.5*log((1+r_data)/(1-r_data)) # Fisher transformation * (n-3)
  
  true <- sqrt(n-3)*0.5*log((1+sqrt(rho))/(1-sqrt(rho)))
  
  # ## Horseshoe plus prior
  # rhs_plus_ref <- stan("mht_rhs_plus.stan", data = list(n=p, y=z_ref, sigma=1), 
  #                      control = list(adapt_delta = 0.95),
  #                      chain = 1, seed = 143246)
  
  
  ## Regularized horseshoe prior
  data_stan <- list()
  p0 <- 1 # prior guess for the number of relevant variables
  tau0 <- p0/(p-p0)
  rhs <- hs(global_scale=tau0)
  data_stan$n <- p
  data_stan$scale_global <- p0/(p-p0) # p0/(p-p0) * 1/sqrt(p)
  data_stan$nu_global <- 1
  data_stan$nu_local <- 1
  data_stan$slab_scale <- rhs$slab_scale
  data_stan$slab_df <- 4
  data_stan$sigma <- 1 # Features are already transformed to have scale = 1 
  
  rhs_ref <- sampling(rhs_model, data = c(data_stan,list(y = z_ref)), 
                      chain = 1, iter = 2000, seed = 38925)  # No problems
  rhs_data <- sampling(rhs_model, data =  c(data_stan,list(y = z_data)),
                       chain = 1, iter = 2000, seed = 38925)  # No problems
  
  theta_RHS_ref <- extract(rhs_ref,pars=paste("theta[",1:p,"]",sep=""))
  theta_RHS_data <- extract(rhs_data,pars=paste("theta[",1:p,"]",sep=""))
  
  SESE_RHS_ref <- theta_RHS_ref %>% 
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    map(~mean(.)) %>%
    reduce(sum)
  SESE_RHS_data <- theta_RHS_data %>% 
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    map(~mean(.)) %>%
    reduce(sum)
  
  SSE_RHS_ref <- theta_RHS_ref %>%
    map(~median(.)) %>%
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    reduce(sum)
  SSE_RHS_data <- theta_RHS_data %>%
    map(~median(.)) %>%
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    reduce(sum)
  
  lambda_RHS_ref <- extract(rhs_ref,pars=paste("lambda[",1:p,"]",sep=""))
  lambda_RHS_data <- extract(rhs_data,pars=paste("lambda[",1:p,"]",sep=""))
  
  tau_RHS_ref <- extract(rhs_ref,par="tau")
  tau_RHS_data <- extract(rhs_data,par="tau")
  
  
  k_RHS_ref <- k_RHS_ref + (lambda_RHS_ref %>% 
                              map(~.*tau_RHS_ref$tau) %>%
                              map(~.^2) %>%
                              map(~. + 1) %>%
                              map(~.^(-1)) %>%
                              map_dbl(~mean(.)))
  
  k_RHS_data <- k_RHS_data + (lambda_RHS_data %>% 
                                map(~.*tau_RHS_data$tau) %>%
                                map(~.^2) %>%
                                map(~. + 1) %>%
                                map(~.^(-1)) %>%
                                map_dbl(~mean(.)))
  
  
  ## Dirichlet-Laplace prior
  DL_ref <- sampling(dl_model, data = list(n = p, y = z_ref, sigma = 1),
                     control = list(adapt_delta = 0.9),
                     chain = 1, iter = 2000, seed = 143246)
  DL_data <- sampling(dl_model, data = list(n = p, y = z_data, sigma = 1),
                      control = list(adapt_delta = 0.9),
                      chain = 1, iter = 2000, seed = 143246)
  
  theta_DL_ref <- extract(DL_ref,pars=paste("theta[",1:p,"]",sep=""))
  theta_DL_data <- extract(DL_data,pars=paste("theta[",1:p,"]",sep=""))
  
  SESE_DL_ref <- theta_DL_ref %>% 
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    map(~mean(.)) %>%
    reduce(sum)
  SESE_DL_data <- theta_DL_data %>% 
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    map(~mean(.)) %>%
    reduce(sum)
  
  SSE_DL_ref <- theta_DL_ref %>%
    map(~median(.)) %>%
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    reduce(sum)
  SSE_DL_data <- theta_DL_data %>%
    map(~median(.)) %>%
    map2(as.list(c(rep(true,k),rep(0,p-k))), ~ ..1 - ..2) %>%
    map(~.^2) %>%
    reduce(sum)
  
  ## Save results
  results <- results %>%
    add_row(error = c(SESE_RHS_ref,
                   SESE_RHS_data,
                   SSE_RHS_ref,
                   SSE_RHS_data),
           error_type = rep(c("SESE","SSE"),each=2),
           prior = rep("RHS",4),
           approach = rep(c("ref","data"),2),
           n = rep(n,4),
           rho = rep(rho,4)) %>%
    add_row(error = c(SESE_DL_ref,
                     SESE_DL_data,
                     SSE_DL_ref,
                     SSE_DL_data),
           error_type = rep(c("SESE","SSE"),each=2),
           prior = rep("DL",4),
           approach = rep(c("ref","data"),2),
           n = rep(n,4),
           rho = rep(rho,4))
  
}

k_result <- tibble(k.value = c(k_RHS_ref/times, k_RHS_data/times),
                   n = rep(n,p*2),
                   rho = rep(rho,p*2),
                   label = rep(c(rep('r',k),rep('s',p-k)),2),
                   approach = c(rep(c('ref','data'),each=p)))

save.image(paste("fullBayes_n",n,"rho",rho,".Rdata",sep=""))

