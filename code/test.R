rm(list=ls())
library(rstan)
library(rstanarm)
library(loo)
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
n <- 251
noise <- array(rnorm((p-13)*n_full), c(n_full,(p-13)))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
pred2 <- c(pred,colnames(dfr[,15:(p+1)]))
fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))

##########################
## SPCA reference model ##
x <- dfr[,-1]
y <- dfr[,1]
nc <- 5
dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p,
           sup.only = T, verbose = F)
z <- predict(dr, x)
model_data <- list(y=y,X=z,p=ncol(z),n=n,s_max=dr$sdev[1])
fit_SPCA <- sampling(fit_model, data = model_data,
                chains = 1, iter = 2000, seed = 45342)
loglik_SPCA <- extract_log_lik(fit_SPCA)
loo_SPCA <- loo(loglik_SPCA)
  
#########################
## RHS reference model ##
p0 <- 5 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
hs_prior <- hs(global_scale=tau0)
fit_RHS <- stan_glm(formula2, data = dfr[sel,], prior = hs_prior,
                    seed=34523, refresh=0)
loo_RHS <- loo(fit_RHS)

compare(loo_SPCA,loo_RHS) # elpd_diff 53.0 (14.5 se), SPCA is better
