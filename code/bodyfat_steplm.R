library(tidyverse)
library(rstan)
library(dimreduce)
set.seed(5345334)

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

p <- 100

fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))

noisy_ref <- numeric(100)
noisy_data <- numeric(100)

rmse_ref <- numeric(100)
rmse_data <- numeric(100)

for(boot in 1:100){
    noise <- array(rnorm((p-13)*n), c(n,(p-13)))
    dfr<-cbind(df,noise=noise)
    formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))
    pred2 <- c(pred,colnames(dfr[,15:(p+1)]))

    sel <- sample(1:n,n,replace=TRUE)
    dfr.train <- dfr[sel,]
    dfr.test  <- dfr[-unique(sel),]

    step_data <- step(lm(formula2,data=dfr.train),direction='backward')

    nc <- 5
    x <- dfr.train[,-1]
    y <- dfr.train$siri
    dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F)
    z <- predict(dr, x)

    ## fit the model

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

    ## Mht
    mufit <- predfun(x)
    yfit <- apply(mufit,1,function(x){ # predict with the posterior predictive mean
        mean(rnorm(1000,x,sigma))
    })

    dfr_ref <- dfr.train
    dfr_ref$siri <- yfit
    step_ref <- step(lm(formula2,data=dfr_ref),direction='backward')

    names_ref <- names(step_ref$coefficients)
    names_data <- names(step_data$coefficients)

    noisy_ref[boot] <- length(grep(pattern='noise',names_ref))
    noisy_data[boot] <- length(grep(pattern='noise',names_data))

    ## testing predictive performance -> RMSE
    y.test.ref <- predict(step_ref,newdata=dfr.test[,-1])
    y.test.data <- predict(step_data,newdata=dfr.test[,-1])

    rmse_ref[boot] <- sqrt(mean((y.test.ref - dfr.test$siri)^2))
    rmse_data[boot] <- sqrt(mean((y.test.data - dfr.test$siri)^2))
}

save.image('bodyfat_step.Rdata')

mean(rmse_ref)
mean(rmse_data)

