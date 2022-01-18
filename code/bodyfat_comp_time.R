library(tidyverse)
library(rstanarm)
library(projpred)
set.seed(453874)

df <- read.table("bodyfat.txt", header = T, sep = ";")
df[,4:19] <- scale(df[,4:19])
df <- as.data.frame(df)
n <- nrow(df)
colnames(df[c("weight_kg", "height")]) <- c("weight", "height")
pred <- c("age", "weight", "height", "neck", "chest", "abdomen", "hip",
          "thigh", "knee", "ankle", "biceps", "forearm", "wrist")
target <- "siri"
formula <- paste("siri~", paste(pred, collapse = "+"))
p <- length(pred)
df <- df[,c(target,pred)]

p0 <- 2 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
rhs_prior <- hs(global_scale=tau0)

## Analysis with original covariates
system.time({fit_ref <- stan_glm(formula, data = df,
                                 prior=rhs_prior, QR=TRUE, seed=8786, refresh=0)})

system.time({proj_path <- cv_varsel(fit_ref, method='forward', cv_method='LOO', nloo=n,
                     verbose = FALSE) 
proj_size <- suggest_size(proj_path, alpha=0.1)})
  
system.time({steplm <- step(lm(formula, data = df), direction = 'backward')})

## Analysis with additional noisy covariates
p <- 100
noise <- array(rnorm((p-13)*n), c(n,p-13))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
noisy_names <- colnames(dfr[,15:(p+1)])

p0 <- 5 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
hs_prior <- hs(global_scale=tau0)

system.time({fit_ref2 <- stan_glm(formula2, data = dfr, prior = hs_prior, QR = TRUE,
                                  seed=34523, refresh=0)})

system.time({proj_path2 <- cv_varsel(fit_ref2, method = 'forward', cv_method = 'LOO',
                                     nloo=n, verbose = FALSE)
nv2 <- suggest_size(proj_path2, alpha=0.1)})

system.time({steplm2 <- step(lm(formula2, data = dfr), direction = 'backward')})
