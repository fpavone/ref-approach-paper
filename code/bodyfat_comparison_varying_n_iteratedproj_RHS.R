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

p <- 100
n <- as.numeric(args[1])
times <- 100
projpred_X <- matrix(0,nrow=times,ncol=p)
alpha_lev <- 0.16
## Stan model
p0 <- 5 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
hs_prior <- hs(global_scale=tau0)
noise <- array(rnorm((p-13)*n_full), c(n_full,(p-13)))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
pred2 <- c(pred,colnames(dfr[,15:(p+1)]))

for(i in 1:times){
    print(paste('ITERATION NUMBER: ',i))
    sel  <- sample(1:n_full,n,replace = T)
    fit <- stan_glm(formula2, data = dfr[sel,], prior = hs_prior,
                    seed=34523, refresh=0)
  ## Iterative projection
  #ref <- get_refmodel(fit)
  #ref.cv <- kfold.refmodel(ref, seed = 5438913, K=10)
  vs_proj <- list()
  suggested_proj <- numeric()
  relevants_proj <- NULL
  costs_proj <- rep(1,p)
  j <- 0
  while(length(relevants_proj)<p){
      j <- j + 1
      vs_proj[[j]] <- cv_varsel(fit, method = 'l1',
                                penalty = costs_proj,
                                cv_method = 'loo',
                                seed = 5438913)
      suggested_proj[j] <- suggest_size(vs_proj[[j]], stat="mlpd",
                                        alpha = alpha_lev, baseline='best')
      if(suggested_proj[j]==0) break
      relevants_proj <- c(relevants_proj,vs_proj[[j]]$vind[1:suggested_proj[j]])
      costs_proj[relevants_proj] <- Inf
  }
  projpred_X[i,relevants_proj] <- 1
}


save.image(file=paste("bodyfat_type1_varyN_p100_iteratedproj_RHS_n",n,".Rdata",sep=""))
