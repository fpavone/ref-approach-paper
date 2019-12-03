## BODYFAT COMPARISON ##
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
library(rstan)
library(rstanarm)
library(projpred)
library(dimreduce)
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

projpred_X <- matrix(0,nrow=times,ncol=p)
alpha_lev <- 0.16
## Stan model
fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))
noise <- array(rnorm((p-13)*n_full), c(n_full,(p-13)))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
pred2 <- c(pred,colnames(dfr[,15:(p+1)]))

for(i in 1:times){
    print(paste('ITERATION NUMBER: ',i))
    sel  <- sample(1:n_full,n,replace = T)
    x <- dfr[sel,-1]
    y <- dfr[sel,1]
    nc <- 5
    dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p,
               sup.only = T, verbose = F)
    z <- predict(dr, x)
  ## fit the model
    cvfun <- function(folds) {
      lapply(1:max(folds), function (k) {
        w <- which(folds!=k)
        dr <- spca(x[w,],y[w], nctot=nc, screenthresh = 0.6, window=p,
                   sup.only = T, verbose = F) # screenthresh = 0.6, sup.only=T
        zz <- predict(dr, x[w,])
        # fit the model
        model_data <- list(y=y[w],X=zz,p=ncol(zz),n=length(w),s_max=dr$sdev[1])
        fit_k <- sampling(fit_model, data = model_data,
                  chains = 1, iter = 2000, seed = 4345342+k,refresh=0)
        draws <- as.matrix(fit_k) # posterior draws
        sigma <- draws[,'sigma'] # noise std
        beta <- draws[,2:(ncol(zz)+1)] # regression coefficients
        alpha <- draws[,'intercept'] # intercept
        predfun <- function(zt){
          colMns <- colMeans(x[w,])
          colSds <- apply(x[w,],2,sd)
          zt <- t(apply(zt,1,function(x){(x-colMns)/colSds}))%*%dr$w
          return(t( beta %*% t(zt) + alpha ))   # n x s
        }
        list(predfun=predfun, dis=sigma)
      })
  }
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
  ## Iterative projection
  ref <- init_refmodel(z=as.matrix(x), y,
                       gaussian(), x=as.matrix(x),
                       predfun=predfun, dis=sigma,
                       cvfun=cvfun)
  ref.cv <- kfold.refmodel(ref, seed = 5438913, K=10)
  vs_proj <- list()
  suggested_proj <- numeric()
  relevants_proj <- NULL
  costs_proj <- rep(1,p)
  j <- 0
  while(length(relevants_proj)<p){
      j <- j + 1
      vs_proj[[j]] <- cv_varsel(ref.cv, method = 'l1',
                                penalty = costs_proj,
                                cv_method = 'kfold',
                                seed = 5438913, K=10)
      suggested_proj[j] <- suggest_size(vs_proj[[j]], stat="mlpd",
                                        alpha = alpha_lev, baseline='best')
      if(suggested_proj[j]==0) break
      relevants_proj <- c(relevants_proj,vs_proj[[j]]$vind[1:suggested_proj[j]])
      costs_proj[relevants_proj] <- Inf
  }
  projpred_X[i,relevants_proj] <- 1
}


save.image(file=paste("bodyfat_type1_varyN_iteratedproj_n",n,".Rdata",sep=""))
