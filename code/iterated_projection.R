library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(projpred)
library(dimreduce)
SEED=1513306866
source("getStability.R")
args = commandArgs(trailingOnly=TRUE)

n <- as.numeric(args[1])
p <- 1000
k <- 100
rho <- as.numeric(args[2])
nc <- 5

times <- 100

fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))

X_proj_test_list <- list()
X_lasso_test_list <- list()
X_proj_cv_list <- list()
X_lasso_cv_list <- list()

alpha_vec <- c(0.05,0.16,0.25) ##c(0.05,0.1,0.16,0.2,0.25)

for(aaa in alpha_vec){
  alpha_lev <- 2*aaa
  X_proj_test <- matrix(0,nrow=times,ncol=p)
  X_lasso_test <- matrix(0,nrow=times,ncol=p)
  X_proj_cv <- matrix(0,nrow=times,ncol=p)
  X_lasso_cv <- matrix(0,nrow=times,ncol=p)
  for(j in 1:times){
    print(j)

    f <- rnorm(n) # the true underlying function
    x <- sqrt(1-rho)*matrix(rnorm(n*k), ncol=k) + sqrt(rho)*f # features are noisy observations from f
    x <- cbind(x,matrix(rnorm(n*(p-k)), ncol=p-k))
    x <- as.data.frame(x)
    y <- f + rnorm(n) # target variable
    names(x) <- c(paste("r",1:k,sep="."),paste("s",(k+1):p,sep="."))

    ntest <- 1000
    ftest <- rnorm(ntest) # the true underlying function
    xtest <- sqrt(1-rho)*matrix(rnorm(ntest*k), ncol=k) + sqrt(rho)*ftest # features are noisy observations from f
    xtest <- cbind(xtest,matrix(rnorm(ntest*(p-k)), ncol=p-k))
    ytest <- ftest + rnorm(ntest) # target variable

    cvfun <- function(folds) {
      lapply(1:max(folds), function (k) {
        w <- which(folds!=k)

        dr <- spca(x[w,],y[w], nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F) # screenthresh = 0.6, sup.only=T
        zz <- predict(dr, x[w,])

        # fit the model
        model_data <- list(y=y[w],X=zz,p=ncol(zz),n=length(w),s_max=dr$sdev[1])
        fit_k <- sampling(fit_model, data = model_data,
                  chains = 1, iter = 2000, seed = 4345342+k,refresh=0)
        draws <- as.matrix(fit_k) # posterior draws
        sigma <- draws[,'sigma'] # noise std
        beta <- draws[,2:(nc+1)] # regression coefficients
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

    # nc <- 5
    dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F) # screenthresh = 0.6, sup.only=T
    z <- predict(dr, x)

    # fit the model
    model_data <- list(y=y,X=z,p=ncol(z),n=n,s_max=dr$sdev[1])
    fit <- sampling(fit_model, data = model_data,
                  chains = 2, iter = 2000, seed = 45342)
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

###### LASSO TEST #######
    ref_lasso <- init_refmodel(z=as.matrix(x), y, gaussian() ,x=as.matrix(x))
    vs_lasso <- list()
    suggested_lasso <- numeric()
    relevants_lasso <- NULL
    costs_lasso <- rep(1,p)

    i <- 0
    while(length(relevants_lasso)<p){
        i <- i + 1
        vs_lasso[[i]] <- varsel(ref_lasso, method = 'l1',
                                penalty = costs_lasso,
                                d_test = list(z=xtest,x=xtest,y=ytest))

        suggested_lasso[i] <- suggest_size(vs_lasso[[i]],stat="rmse",alpha=alpha_lev,baseline='best')

        if(suggested_lasso[i]==0) break

        relevants_lasso <- c(relevants_lasso,vs_lasso[[i]]$vind[1:suggested_lasso[i]])
        costs_lasso[relevants_lasso] <- Inf
    }
    X_lasso_test[j,relevants_lasso] <- 1

###### PROJ TEST #######
    ref <- init_refmodel(z=as.matrix(x), y,
                         gaussian(), x=as.matrix(x),
                         predfun=predfun, dis=sigma,
                         cvfun=cvfun)
    vs_proj <- list()
    suggested_proj <- numeric()
    relevants_proj <- NULL
    costs_proj <- rep(1,p)

    i <- 0
    while(length(relevants_proj)<p){
        i <- i + 1
        vs_proj[[i]] <- varsel(ref, method = 'l1',
                               penalty = costs_proj,
                               d_test = list(z=xtest,x=xtest,y=ytest))

        suggested_proj[i] <- suggest_size(vs_proj[[i]], stat="mlpd", alpha = alpha_lev, baseline='best')

        if(suggested_proj[i]==0) break

        relevants_proj <- c(relevants_proj,vs_proj[[i]]$vind[1:suggested_proj[i]])
        costs_proj[relevants_proj] <- Inf
    }
    X_proj_test[j,relevants_proj] <- 1

###### LASSO CV ######
    ref_lasso <- init_refmodel(z=as.matrix(x), y, gaussian() ,x=as.matrix(x))
    vs_lasso <- list()
    suggested_lasso <- numeric()
    relevants_lasso <- NULL
    costs_lasso <- rep(1,p)

    i <- 0
    while(length(relevants_lasso)<p){
        i <- i + 1
        vs_lasso[[i]] <- cv_varsel(ref_lasso, method = 'l1',
                                penalty = costs_lasso,
                                cv_method = 'kfold',
                                K = 10,
                                seed = 589435)

        suggested_lasso[i] <- suggest_size(vs_lasso[[i]],stat="rmse",alpha=alpha_lev,baseline='best')

        if(suggested_lasso[i]==0) break

        relevants_lasso <- c(relevants_lasso,vs_lasso[[i]]$vind[1:suggested_lasso[i]])
        costs_lasso[relevants_lasso] <- Inf
    }
    X_lasso_cv[j,relevants_lasso] <- 1

###### PROJ CV #######
    ref <- init_refmodel(z=as.matrix(x), y,
                         gaussian(), x=as.matrix(x),
                         predfun=predfun, dis=sigma,
                         cvfun=cvfun)
    ref.cv <- kfold.refmodel(ref, seed = 5438913, K=10)
    vs_proj <- list()
    suggested_proj <- numeric()
    relevants_proj <- NULL
    costs_proj <- rep(1,p)

    i <- 0
    while(length(relevants_proj)<p){
        i <- i + 1
        vs_proj[[i]] <- cv_varsel(ref.cv, method = 'l1',
                                  penalty = costs_proj,
                                  cv_method = 'kfold',
                                  seed = 5438913, K=10)

        suggested_proj[i] <- suggest_size(vs_proj[[i]], stat="mlpd", alpha = alpha_lev, baseline='best')

        if(suggested_proj[i]==0) break

        relevants_proj <- c(relevants_proj,vs_proj[[i]]$vind[1:suggested_proj[i]])
        costs_proj[relevants_proj] <- Inf
    }
    X_proj_cv[j,relevants_proj] <- 1
  }
  X_proj_test_list <- c(X_proj_test_list,list(X_proj_test))
  X_lasso_test_list <- c(X_lasso_test_list,list(X_lasso_test))
  X_proj_cv_list <- c(X_proj_cv_list, list(X_proj_cv))
  X_lasso_cv_list <- c(X_lasso_cv_list, list(X_lasso_cv))
}


save.image(file=paste("iterated_n",n,"_rho",rho,"_3alphas_suponlyT.Rdata",sep=""))
