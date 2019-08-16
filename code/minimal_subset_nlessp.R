library(tidyverse)
library(projpred)
library(dimreduce)
library(rstan)
library(rstanarm)
library(loo)
options(mc.cores = parallel::detectCores())
set.seed(453876)

args = commandArgs(trailingOnly=TRUE) # n rho nc

## Simulation parameters
n <- as.numeric(args[1])  # number of observations
rho <- as.numeric(args[2])  # correlation level
p <- 100   # total number of features
k <- 20    # number of relevant features

times <- 100

## Data simulation mechanism
simulate_data <- function(n){ # n, rho, p and k from global environment
  f <- rnorm(n) # the true underlying function
  data <- as_tibble(f + rnorm(n)) %>%   # observed targe variable y
    bind_cols(as_tibble(sqrt(1-rho)*matrix(rnorm(n*k), ncol=k) + sqrt(rho)*f)) %>%  # set of relevant covariates
    bind_cols(as_tibble(matrix(rnorm(n*(p-k)), ncol=p-k))) %>%  # set of spurious covariates
    set_names(c("y",paste("r.",1:k,sep=""),paste("s.",(k+1):p,sep="")))

  return(data)
}
## Outputs
## X.step.ref <- matrix(0,nrow=times,ncol=p)
## X.step.data <- matrix(0,nrow=times,ncol=p)
X.bayes.step <- matrix(0,nrow=times,ncol=p)
X.projpred <- matrix(0,nrow=times,ncol=p)
rmse.projpred <- numeric(times)
## rmse.step.data <- numeric(times)
## rmse.step.ref <- numeric(times)
rmse.bayes.step <- numeric(times)
NA_count <- 0
## Experiment
for(i in 1:times){
    print(paste('iteration:',i))
    data <- simulate_data(n)
    test <- simulate_data(1000)
    y <- dplyr::pull(data, var="y")
    x <- dplyr::select(data, -contains("y"))
    y.test <- dplyr::pull(test, var="y")
    x.test <- dplyr::select(test, -contains("y"))
    ## Fitting the reference model
    nc <- 5
    dr <- spca(x,y, nctot=nc, screenthresh = 0.6,
               window=p, sup.only = T, verbose = F)
    z <- predict(dr, x)
    model_data <- list(y = y,
                       X = z,
                       p = nc,
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
    cvfun <- function(folds) {
        lapply(1:max(folds), function (k) {
            w <- which(folds!=k)
            dr <- spca(x[w,],y[w], nctot=nc, screenthresh = 0.6,
                       window=p, sup.only = T, verbose = F)
            zz <- predict(dr, x[w,])
            ## fit the model
            model_data <- list(y=y[w],X=zz,p=nc,n=length(w),s_max=dr$sdev[1])
            fit_k <- stan("model.stan", data = model_data, chains=1, # form the model fit using only data from indices folds!=k
                          seed=(9834574+k), refresh=0)
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
    ref <- init_refmodel(z=as.matrix(x), y,
                         gaussian(), x=as.matrix(x),
                         predfun=predfun, dis=sigma,
                         cvfun=cvfun)
    ## features <- names(x)
    ## formula.model <- paste(features, collapse='+')
    ## dfr <- as.data.frame(cbind(y,x))
    ## p0 <- 5 # prior guess for the number of relevant variables
    ## tau0 <- p0/(p-p0) * 1/sqrt(n)
    ## rhs_prior <- hs(global_scale=tau0)
    ## fit <- stan_glm(formula=paste('y',formula.model,sep='~'), family = gaussian(),
    ##                 data = dfr, prior = rhs_prior,
    ##                QR = TRUE, chain = 2)
    ## ref <- fit
    ## projpred
    vs <- cv_varsel(ref, method='forward', cv_method = 'kfold', K=10)
    ## vs <- cv_varsel(ref, method='forward', cv_method = 'loo')
    suggested <- suggest_size(vs, stat='mlpd', alpha=0.16, baseline='ref')
    if(is.na(suggested)){
        suggested <- suggest_size(vs, stat='mlpd', alpha=0.16, baseline='best')
        NA_count <- NA_count + 1
    }
    sub.proj <- proj_predict(vs,
                             xnew = as.matrix(x.test[,vs$vind[1:suggested]]),
                             vind = vs$vind[1:suggested])
    X.projpred[i,vs$vind[1:suggested]] <- 1
    rmse.projpred[i] <- sqrt(mean((apply(sub.proj,2,mean)-y.test)^2))
    ## steplm data/ref
    ## mufit <- predfun(x)
    ## yfit <- apply(mufit,1,function(x){ # predict with the posterior predictive mean
    ##    mean(rnorm(1000,x,sigma))
    ## })
    ## ## yfit <- apply(posterior_predict(fit),2,mean)
    dfr <- as.data.frame(cbind(y,x))
    ## dfr.ref <- dfr
    ## dfr.ref$y <- yfit
    ## step.ref <- step(lm('y~.',data=dfr.ref),direction='backward')
    ## step.data <- step(lm('y~.',data=dfr),direction='backward')
    ## y.test.ref <- predict(step.ref,newdata=x.test)
    ## y.test.data <- predict(step.data,newdata=x.test)
    ## rmse.step.ref[i] <- sqrt(mean((y.test.ref - y.test)^2))
    ## rmse.step.data[i] <- sqrt(mean((y.test.data - y.test)^2))
    ## sel.step.ref <- as.numeric(gsub('[^0-9]', '',
    ##                                 names(step.ref$coefficients))[-1])
    ## sel.step.data <- as.numeric(gsub('[^0-9]', '',
    ##                                  names(step.data$coefficients))[-1])
    ## X.step.ref[i, sel.step.ref] <- 1
    ## X.step.data[i, sel.step.data] <- 1
    ## bayes.steplm
    p0 <- 5 # prior guess for the number of relevant variables
    tau0 <- p0/(p-p0) * 1/sqrt(n)
    rhs_prior <- hs(global_scale=tau0)
    features <- names(dfr)[-1]
    formula.model <- paste(features, collapse='+')
    current <- stan_glm(formula=paste('y',formula.model,sep='~'),
                        family = gaussian(),
                        data = dfr, prior = rhs_prior,
                        QR = TRUE, chain = 2)
    loo.current <- loo(current)
    escape <- 0
    while(escape<k){
        escape <- escape + 1
        ## computing bayesian pvalues
        bpval <- extract(current$stanfit, pars=features) %>%
            map_dbl(function(x) min(sum(x<0),sum(x>0))/length(x))
        exclude <- which.max(bpval)
        location.exclude <- which(features==names(exclude))
        features <- features[-location.exclude]
        formula.model <- paste(features, collapse='+')
        proposed <- stan_glm(formula=paste('y',formula.model,sep='~'),
                             family = gaussian(),
                             data = dfr, prior = rhs_prior,
                             QR = TRUE, chain = 2)
        loo.proposed <- loo(proposed)
        comparison <- loo_compare(loo.current,loo.proposed)
        if(comparison['model1','elpd_diff']==0) break
        current <- proposed
    }
    predict_dist <- posterior_predict(current, newdata = x.test)
    rmse.bayes.step[i] <- sqrt(mean((apply(predict_dist,2,mean)-y.test)^2))
    sel.bayes.step <- as.numeric(gsub('[^0-9]', '', c(names(exclude),features)))
    X.bayes.step[i, sel.bayes.step] <- 1
}

save.image(file = paste("minimal_subset_SPC_nlessp_n",n,"rho",rho,".RData",sep=""))
