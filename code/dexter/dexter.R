library(tidyverse)
library(rstan)
library(dimreduce)
library(rstanarm)

load('dexter.RData')

## Original dataset:
## - 9947 features (of which 2562 always zero) that represent frequencies of
##   occurences of word stems in text
## - target: binary variable indexing which article (statistical unit) is about
##   corporate acquisitions
##
## The frequency of appearance of words in text is known to follow Zipf's law.
## Such law has been used to input 10053 features drawn at random.


dexter <- x
target <- y

n <- nrow(dexter)
p <- ncol(dexter)

hist(cor(as.matrix(x),y))

true <- cor(as.matrix(x),y)

sub <- sample(1:n,n/20)

sample <- cor(as.matrix(x[sub,]),y[sub])

qplot(y=true,x=sample) +
    geom_abline(intercept=0,slope=1,linetype='dashed')

## Define reference model

x.sample <- as.matrix(x[sub,])
y.sample <- y[sub]

## Stan models
nc <- 50
dr <- spca(x.sample,y.sample, nctot=nc, screenthresh = NULL, window=p, sup.only = F, verbose = F)
z <- predict(dr, x.sample)

## fit the model
fit <- stan_glm('y~.',data=data.frame(y=y.sample,z=z), family=binomial(link=logit),
                prior = normal(), prior_intercept=normal(),
                chains=2)


bayesplot::mcmc_intervals(as.matrix(fit))
## In-sample predictions

yfit <- apply(posterior_predict(fit),2,mean)

ref <- cor(x.sample,yfit)

qplot(x=sample,y=true) +
    geom_point(data=data.frame(ref=ref,true=true),aes(x=ref,y=true),color='red') +
    geom_abline(intercept=0,slope=1,linetype='dashed')

## Out-of-sample predictions

yfit <- numeric(n/20)
for(i in 1:(n/20)){
  nc <- 10
  dr <- spca(x.sample[-i,],y.sample[-i], nctot=nc, screenthresh = NULL, window = p, sup.only = F, verbose = F)
  z <- predict(dr, x.sample[-i,])
  fit <- stan_glm('y~.', family = binomial(link="logit"), data = data.frame(y=y.sample[-i],X=z),
                    prior = normal(0,5), prior_intercept = normal(0,5), chains=2)
  colMns <- colMeans(x.sample[-i,])
  colSds <- apply(x.sample[-i,],2,sd)
  z_new <- t(apply(rbind(x.sample[i,]),1,function(x){ifelse(colSds==0,x,(x-colMns)/colSds)}))%*%dr$w
  yfit[i] <- apply(posterior_predict(fit, newdata=data.frame(y=y.sample[i],X=z_new)),2,mean)
}

length(y.sample)


#### Bodyfat check
load('../bodyfat_type1_varyP_p5000.Rdata')

names(rhs_ref)

hist(x=apply(extract(rhs_ref, pars='lambda')$lambda,2,mean),breaks=50)
