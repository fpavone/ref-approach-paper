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

lasso_X <- matrix(0,nrow=times,ncol=p)
alpha_lev <- 0.16

noise <- array(rnorm((p-13)*n_full), c(n_full,(p-13)))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
pred2 <- c(pred,colnames(dfr[,15:(p+1)]))

for(i in 1:times){
  print(paste('ITERATION NUMBER: ',i))
  sel  <- sample(1:n_full,n,replace = T)
  x <- dfr[sel,-1]
  y <- dfr[sel,1]
  
  ## Iterative lasso
  ref_lasso <- init_refmodel(z=as.matrix(x), y, gaussian() ,x=as.matrix(x))
  vs_lasso <- list()
  suggested_lasso <- numeric()
  relevants_lasso <- NULL
  costs_lasso <- rep(1,p)
  
  j <- 0
  while(length(relevants_lasso)<p){
    j <- j + 1
    vs_lasso[[j]] <- cv_varsel(ref_lasso, method = 'l1',
                               penalty = costs_lasso,
                               cv_method = 'kfold',
                               K = 15,
                               seed = 589435)
    
    suggested_lasso[j] <- suggest_size(vs_lasso[[j]],stat="rmse",alpha=alpha_lev,baseline='best')
    
    if(suggested_lasso[j]==0) break
    
    relevants_lasso <- c(relevants_lasso,vs_lasso[[j]]$vind[1:suggested_lasso[j]])
    costs_lasso[relevants_lasso] <- Inf
  }
  lasso_X[i,relevants_lasso] <- 1
}


save.image(file=paste("bodyfat_type1_varyN_p100_iteratedlasso_n",n,".Rdata",sep=""))
