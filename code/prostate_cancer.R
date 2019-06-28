library(tidyverse)
library(corrplot)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(loo)
library(dimreduce)
library(locfdr)
SEED=14124869

# file preview shows a header row
diabetes <- read.csv("prostmat.csv", header = TRUE)

data <- diabetes %>% 
  rownames_to_column() %>%
  gather(target,value,-rowname) %>%
  spread(rowname,value) %>%
  mutate_at('target',~as.numeric(grepl('cancer',.)))

order <- as.numeric(colnames(data)[-1])
  
n <- nrow(data)
p <- ncol(data)-1

x <- data[,-1]
y <- data$target

nc <- 10

dr1 <- spca(as.matrix(data[,-1]),data$target, nctot=1, screenthresh = NULL, window = p, sup.only = F, verbose = F)
aa <- prcomp(as.matrix(data[,-1]))
ggplot(data.frame(spc1=aa$x[,1],spc2=aa$x[,2],label=as.factor(y)),aes(x=spc1,y=spc2,color=label)) + geom_point()


dr <- spca(as.matrix(data[,-1]),data$target, nctot=nc, screenthresh = NULL, window = p, sup.only = T, verbose = F)
z <- predict(dr, x)

# fit the model
fit <- stan_glm('y~.', family = binomial(link="logit"), data = data.frame(y=y,X=z),
                prior = normal(), prior_intercept = normal(), chains=4)

folds <- sample(1:5,n,replace=TRUE)
acc <- numeric(5)
for(k in 1:5){
  sel <- which(folds==k)
  x_k <- x[sel,]
  y_k <- y[sel]
  
  nc <- 10
  dr_k <- spca(as.matrix(x_k),y_k, nctot=nc, screenthresh = NULL, window = p, sup.only = T, verbose = F)
  z_k <- predict(dr_k, x_k)
  
  fit_k <- stan_glm('y~.', family = binomial(link="logit"), data = data.frame(y=y_k,X=z_k),
                  prior = normal(0,5), prior_intercept = normal(0,5), chains=2)

  colMns <- colMeans(x_k)
  colSds <- apply(x_k,2,sd)
  z_new <- t(apply(x[-sel,],1,function(x){(x-colMns)/colSds}))%*%dr_k$w
  
  yfit_k <- apply(posterior_predict(fit_k, newdata=data.frame(y=y[-sel],X=z_new)),2,mean)

  acc[k] <- sum(round(yfit_k)==y[-sel])/length(y[-sel])
}


yfit <- numeric(n)
for(i in 1:n){
  nc <- 10
  dr <- spca(as.matrix(x[-i,]),y[-i], nctot=nc, screenthresh = 0.7, window = p, sup.only = T, verbose = F)
  z <- predict(dr, x[-i,])
  
  fit <- stan_glm('y~.', family = binomial(link="logit"), data = data.frame(y=y[-i],X=z),
                    prior = normal(0,5), prior_intercept = normal(0,5), chains=2)
  
  colMns <- colMeans(x[-i,])
  colSds <- apply(x[-i,],2,sd)
  z_new <- t(apply(x[i,],1,function(x){(x-colMns)/colSds}))%*%dr$w
  
  yfit[i] <- apply(posterior_predict(fit, newdata=data.frame(y=y[i],X=z_new)),2,mean)
}

r_ref <- as.vector(cor(yfit,x))
z_ref <- sqrt(n-3) * 0.5*log((1+r_ref)/(1-r_ref)) # Fisher transformation * (n-3)

r_data <- as.vector(cor(y,x))
z_data <- sqrt(n-3) * 0.5*log((1+r_data)/(1-r_data)) # Fisher transformation * (n-3)



z_fdr <- read_csv('prostz.txt',col_names=FALSE) %>% pull(X1)
lcfdr <- locfdr(z_fdr)
sel <- which(lcfdr$fdr<=0.2)


data_plot <- tibble(z_ref=z_ref,
                    z_data=z_data,
                    sel.lfdr=numeric(p))
data_plot$sel.lfdr[which((order %in% sel ==TRUE))] <- 1

ggplot(data_plot,aes(x=z_data,y=z_ref,color=as.factor(sel.lfdr))) + 
  geom_point(size=0.7,alpha=1) +
  geom_abline(slope=1,intercept=0,linetype='dashed',size=0.4) +
  # scale_x_continuous(limits=c(2.5,6)) +
  # scale_y_continuous(limits=c(2.5,6)) +
  labs(color='Locfdr selection') +
  theme_light()

sel_fdr <- which((order %in% sel ==TRUE))

lcfdr_ref <- locfdr(z_ref)
sel_ref <- which(lcfdr_ref$fdr<=0.2)

lcfdr_data <- locfdr(z_data)
sel_data <- which(lcfdr_data$fdr<=0.2)


sel_fdr %in% sel_ref
sel_fdr %in% sel_data
