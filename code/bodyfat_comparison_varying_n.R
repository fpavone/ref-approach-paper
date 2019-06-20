## BODYFAT COMPARISON ##
rm(list=ls())
args = commandArgs(trailingOnly=TRUE) 
library(rstan)
library(rstanarm)
library(locfdr)
library(dimreduce)
library(EbayesThresh)
set.seed(3428954)

df <- read.table("bodyfat.txt", header = T, sep = ";")
df[,4:19] <- scale(df[,4:19])
df <- as.data.frame(df)
n_full <- nrow(df)
colnames(df[c("weight_kg", "height")]) <- c("weight", "height")
pred <- c("age", "weight", "height", "neck", "chest", "abdomen", "hip", 
          "thigh", "knee", "ankle", "biceps", "forearm", "wrist")
target <- "siri"
formula <- paste("siri~", paste(pred, collapse = "+"))
p <- length(pred)
df <- df[,c(target,pred)]

p <- 1000
n <- as.numeric(args[1])
times <- 100

## Control of the local false discovery rate
lfdr_X_ref <- matrix(0,nrow=times,ncol=p)
lfdr_X_data <- matrix(0,nrow=times,ncol=p)

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

## Stan models
fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))
rhs_model <- stan_model("mht_rhs_sigmaprior.stan", model_name = "rhs_model",
                        auto_write = rstan_options("auto_write"=TRUE))

noise <- array(rnorm((p-13)*n_full), c(n_full,(p-13)))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
pred2 <- c(pred,colnames(dfr[,15:(p+1)]))

for(i in 1:times){
  sel  <- sample(1:n_full,n,replace = T)
  
  x <- dfr[sel,-1]
  y <- dfr[sel,1]

  nc <- 5
  dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F)
  z <- predict(dr, x)
  
  # fit the model
  model_data <- list(y=y,X=z,p=ncol(z),n=n,s_max=dr$sdev[1])
  fit <- sampling(fit_model, data = model_data, 
                  chains = 1, iter = 2000, seed = 45342)  
  predfun <- function(zt){
    colMns <- colMeans(x)
    colSds <- apply(x,2,sd)
    zt <- t(apply(zt,1,function(x){(x-colMns)/colSds}))%*%dr$w
    return(t( beta %*% t(zt) + alpha ))   # n x s
  }
  
  draws <- as.matrix(fit) # posterior draws
  sigma <- draws[,'sigma'] # noise std
  beta <- draws[,2:(ncol(z)+1)] # regression coefficients
  alpha <- draws[,'intercept'] # intercept
  
  # Mht
  mufit <- predfun(x)
  yfit <- apply(mufit,1,function(x){ # predict with the posterior predictive mean
    mean(rnorm(1000,x,sigma))
  })
  
  r_ref <- as.vector(cor(yfit,x))
  z_ref <- sqrt(n-3) * 0.5*log((1+r_ref)/(1-r_ref)) # Fisher transformation * (n-3)
  
  r_data <- as.vector(cor(y,x))
  z_data <- sqrt(n-3) * 0.5*log((1+r_data)/(1-r_data)) # Fisher transformation * (n-3)
  
  ## Control of of the local false discovery rate
  lfdr_ref <- locfdr(z_ref, nulltype = 1, plot = 0, bre = 120, type = 1)
  lfdr_data <- locfdr(z_data, nulltype = 1, plot = 0, bre = 120, type = 1)
  
  lfdr_sel_ref <- which(lfdr_ref$fdr<=0.2)
  lfdr_sel_data <- which(lfdr_data$fdr<=0.2)
  
  lfdr_X_ref[i,lfdr_sel_ref] <- 1
  lfdr_X_data[i,lfdr_sel_data] <- 1
  
  ## Empirical Bayes median thresholding
  z.est_ref <- ebayesthresh(z_ref, prior = "laplace", a = NA,
                            sdev = NA, verbose = TRUE, threshrule = "median") 
  z.est_data <- ebayesthresh(z_data, prior = "laplace", a = NA,
                             sdev = NA, verbose = TRUE, threshrule = "median") 
  
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
 # data_stan$sigma <- 1 # Features are already transformed to have scale = 1 
  
  rhs_ref <- sampling(rhs_model, data = c(data_stan,list(y = z_ref)), 
                      chain = 1, iter = 2000, seed = 38925,
                      control = list(adapt_delta = 0.95))  # No problems
  rhs_data <- sampling(rhs_model, data =  c(data_stan,list(y = z_data)),
                       chain = 1, iter = 2000, seed = 38925,
                       control = list(adapt_delta = 0.95))  # No problems
  
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
}


save.image(file=paste("bodyfat_type1_varyN_n",n,".Rdata",sep=""))