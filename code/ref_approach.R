library(dplyr)
library(purrr)
library(loo)
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

  data <- as_tibble(f + rnorm(n)) %>%   # observed targe variable y
    bind_cols(as_tibble(sqrt(1-rho)*matrix(rnorm(n*k), ncol=k) + sqrt(rho)*f)) %>%  # set of relevant covariates
    bind_cols(as_tibble(matrix(rnorm(n*(p-k)), ncol=p-k))) %>%  # set of spurious covariates
    set_names(c("y",paste("r.",1:k,sep=""),paste("s.",(k+1):p,sep="")))

  return(data)
}


## Control of the local false discovery rate
lfdr_X_ref <- matrix(0,nrow=times,ncol=p)
lfdr_X_data <- matrix(0,nrow=times,ncol=p)

lfdr_mean_ref <- numeric(p)
lfdr_mean_data <- numeric(p)

## Empirical Bayes median thresholding
ebmt_X_ref <- matrix(0,nrow=times,ncol=p)
ebmt_X_data <- matrix(0,nrow=times,ncol=p)

## Confidence intervals inclusion probabilities
ci80_X_ref <- matrix(0,nrow=times,ncol=p)
ci90_X_ref <- matrix(0,nrow=times,ncol=p)
ci95_X_ref <- matrix(0,nrow=times,ncol=p)
ci80_X_data <- matrix(0,nrow=times,ncol=p)
ci90_X_data <- matrix(0,nrow=times,ncol=p)
ci95_X_data <- matrix(0,nrow=times,ncol=p)

## Control of the Q-value
rho_hat_seq <- c(0.05,0.1,0.2,0.3)
q_avg_ref <- matrix(0,nrow=p,ncol=length(rho_hat_seq))
q_avg_data  <- matrix(0,nrow=p,ncol=length(rho_hat_seq))

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
                     p = ncol(z),
                     n = n,
                     s_max = dr$sdev[1])
  fit <- stan("model.stan", data = model_data,
              chains=1, seed=45342)

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

  ## Control of of the local false discovery rate
  lfdr_ref <- locfdr(z_ref, nulltype = 0, plot = 0)
  lfdr_data <- locfdr(z_data, nulltype = 0, plot = 0)

  lfdr_mean_ref <- lfdr_mean_ref + sort(lfdr_ref$fdr,decreasing=F)/times
  lfdr_mean_data <- lfdr_mean_data + sort(lfdr_data$fdr,decreasing=F)/times

  lfdr_sel_ref <- which(lfdr_ref$fdr<0.2)
  lfdr_sel_data <- which(lfdr_data$fdr<0.2)

  lfdr_X_ref[i,lfdr_sel_ref] <- 1
  lfdr_X_data[i,lfdr_sel_data] <- 1

  ## Empirical Bayes median thresholding
  z.est_ref <- ebayesthresh(z_ref, prior = "laplace", a = NA,
                            sdev = 1, verbose = TRUE, threshrule = "median")
  z.est_data <- ebayesthresh(z_data, prior = "laplace", a = NA,
                             sdev = 1, verbose = TRUE, threshrule = "median")

  ebmt_sel_ref <- which(z.est_ref$muhat!=0)
  ebmt_sel_data <- which(z.est_data$muhat!=0)

  ebmt_X_ref[i,ebmt_sel_ref] <- 1
  ebmt_X_data[i,ebmt_sel_data] <- 1

  ## Credibility intervals inclusion probabilities with regularized horseshoe prior
  data_stan <- list()
  p0 <- 1 # prior guess for the number of relevant variables
  tau0 <- p0/(p-p0) * 1/sqrt(p)
  rhs <- hs(global_scale=tau0)
  data_stan$n <- p
  data_stan$scale_icept <- rhs$scale
  data_stan$scale_global <- rhs$global_scale
  data_stan$nu_global <- rhs$global_df
  data_stan$nu_local <- rhs$df
  data_stan$slab_scale <- rhs$slab_scale
  data_stan$slab_df <- rhs$slab_df
  data_stan$sigma <- 1 # Features are already transformed to have scale = 1

  rhs_ref <- stan(file="mht_rhs.stan", data = c(data_stan,list(y=z_ref)), chain = 1)  # No problems
  rhs_data <- stan(file="mht_rhs.stan", data =  c(data_stan,list(y=z_data)), chain = 1) # Divergent transitions (4)

  theta_ref <- as.matrix(extract(rhs_ref,pars="theta")$theta)
  colnames(theta_ref) <- paste("theta",1:p,sep=".")
  theta_data <- as.matrix(extract(rhs_data,pars="theta")$theta)
  colnames(theta_data) <- paste("theta",1:p,sep=".")

  rho_ref <- apply(theta_ref/sqrt(n-3),2,tanh)  # transformed back in correlation scale
  rho_data <- apply(theta_data/sqrt(n-3),2,tanh) # transofrmed back in correlation scale

  ci80_l_ref <- apply(rho_ref,2,function(x){quantile(x,probs=0.1)})
  ci80_u_ref <- apply(rho_ref,2,function(x){quantile(x,probs=0.9)})
  ci90_l_ref <- apply(rho_ref,2,function(x){quantile(x,probs=0.05)})
  ci90_u_ref <- apply(rho_ref,2,function(x){quantile(x,probs=0.95)})
  ci95_l_ref <- apply(rho_ref,2,function(x){quantile(x,probs=0.025)})
  ci95_u_ref <- apply(rho_ref,2,function(x){quantile(x,probs=0.975)})

  ci80_l_data <- apply(rho_data,2,function(x){quantile(x,probs=0.1)})
  ci80_u_data <- apply(rho_data,2,function(x){quantile(x,probs=0.9)})
  ci90_l_data <- apply(rho_data,2,function(x){quantile(x,probs=0.05)})
  ci90_u_data <- apply(rho_data,2,function(x){quantile(x,probs=0.95)})
  ci95_l_data <- apply(rho_data,2,function(x){quantile(x,probs=0.025)})
  ci95_u_data <- apply(rho_data,2,function(x){quantile(x,probs=0.975)})

  rhs_sel80_ref <- which(ci80_l_ref*ci80_u_ref>0)
  rhs_sel90_ref <- which(ci90_l_ref*ci90_u_ref>0)
  rhs_sel95_ref <- which(ci95_l_ref*ci95_u_ref>0)
  rhs_sel80_data <- which(ci80_l_data*ci80_u_data>0)
  rhs_sel90_data <- which(ci90_l_data*ci90_u_data>0)
  rhs_sel95_data <- which(ci95_l_data*ci95_u_data>0)

  ci80_X_ref[i,rhs_sel80_ref] <- 1
  ci90_X_ref[i,rhs_sel90_ref] <- 1
  ci95_X_ref[i,rhs_sel95_ref] <- 1
  ci80_X_data[i,rhs_sel80_data] <- 1
  ci90_X_data[i,rhs_sel90_data] <- 1
  ci95_X_data[i,rhs_sel95_data] <- 1

  ## Control of the Q-value with regularized horseshoe prior
  q_ref <- NULL
  q_data <- NULL

  fdr_ref <- NULL
  fdr_data <- NULL

  for(rho_hat in rho_hat_seq) {
    PEP_ref <- apply(rho_ref,2,function(x){sum(abs(x)<=rho_hat)/length(x)})
    PEP_data <- apply(rho_data,2,function(x){sum(abs(x)<=rho_hat)/length(x)})

    sorted_PEP_ref <- sort(PEP_ref, decreasing = F, index.return=T)
    sorted_PEP_data <- sort(PEP_data, decreasing = F, index.return=T)

    q_ref <- cbind(q_ref,cumsum(sorted_PEP_ref$x)/(1:p))
    q_data <- cbind(q_data,cumsum(sorted_PEP_data$x)/(1:p))

    fdr_ref <- c(fdr_ref,cumsum(sorted_PEP_ref$ix>k)/(1:p))
    fdr_data <- c(fdr_data,cumsum(sorted_PEP_data$ix>k)/(1:p))
  }

  q_avg_ref <- q_avg_ref + q_ref/times
  q_avg_data <- q_avg_data + q_data/times
}

save.image(file = paste("ref_approach_n",n,"rho",rho,".RData"))
