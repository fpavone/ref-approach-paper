# library(dplyr)
# library(purrr)
# library(tibble)
library(tidyverse)
library(dimreduce)
library(locfdr)
library(EbayesThresh)
library(rstan)
library(rstanarm)
library(latex2exp)
options(mc.cores = parallel::detectCores())
set.seed(453876)

## Simulation parameters
n <- 50  # number of observations
rho <- 0.3  # correlation level
p <- 1000   # total number of features
k <- 100    # number of relevant features

## Data simulation mechanism
simulate_data <- function(){ # n, rho, p and k from global environment
  f <- rnorm(n) # the true underlying function
  
  data <- tibble(y = f + rnorm(n)) %>%   # observed targe variable y
    bind_cols(as_tibble(sqrt(1-rho)*matrix(rnorm(n*k), ncol=k) + sqrt(rho)*f)) %>%  # set of relevant covariates
    bind_cols(as_tibble(matrix(rnorm(n*(p-k)), ncol=p-k))) %>%  # set of spurious covariates
    set_names(c("y",paste("r.",1:k,sep=""),paste("s.",(k+1):p,sep="")))
  
  return(data)
}

## Stan models
fit_model <- stan_model("model.stan", model_name = "reference_model",
                        auto_write = rstan_options("auto_write"=TRUE))
rhs_model <- stan_model("mht_rhs.stan", model_name = "rhs_model",
                        auto_write = rstan_options("auto_write"=TRUE))

k_RHS_ref <- numeric(p)
k_RHS_data <- numeric(p)

for(i in 1:10){
data <- simulate_data()
y <- dplyr::pull(data, var="y")
x <- dplyr::select(data, -contains("y"))


## Fitting the reference model
nc <- 5
dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F)
z <- predict(dr, x)

model_data <- list(y = y,
                   X = z,
                   p = nc,
                   n = n,
                   s_max = dr$sdev[1])
fit <- sampling(fit_model, data = model_data, 
                chains = 1, iter = 2000, seed = 45342)

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

mufit <- predfun(x)
yfit <- apply(mufit,1,function(x){ # predict with the posterior predictive mean
  mean(rnorm(1000,x,sigma))
})

## Computing sample correlation and z-values using reference model's predictions and observed data
r_ref <- as.vector(cor(yfit,x))
z_ref <- sqrt(n-3) * 0.5*log((1+r_ref)/(1-r_ref)) # Fisher transformation * (n-3)

r_data <- as.vector(cor(y,x))
z_data <- sqrt(n-3) * 0.5*log((1+r_data)/(1-r_data)) # Fisher transformation * (n-3)

true <- sqrt(n-3)*0.5*log((1+sqrt(rho))/(1-sqrt(rho)))

# ## Horseshoe plus prior
# rhs_plus_ref <- stan("mht_rhs_plus.stan", data = list(n=p, y=z_ref, sigma=1), 
#                      control = list(adapt_delta = 0.95),
#                      chain = 1, seed = 143246)


## Regularized horseshoe prior
data_stan <- list()
p0 <- 1 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0)
rhs <- hs(global_scale=tau0)
data_stan$n <- p
data_stan$scale_global <- p0/(p-p0) # p0/(p-p0) * 1/sqrt(p)
data_stan$nu_global <- 1
data_stan$nu_local <- 1
data_stan$slab_scale <- rhs$slab_scale
data_stan$slab_df <- 4
data_stan$sigma <- 1 # Features are already transformed to have scale = 1 

rhs_ref <- sampling(rhs_model, data = c(data_stan,list(y = z_ref)), 
                    chain = 1, iter = 2000, seed = 38925)  # No problems
rhs_data <- sampling(rhs_model, data =  c(data_stan,list(y = z_data)),
                     chain = 1, iter = 2000, seed = 38925)  # No problems

lambda_RHS_ref <- extract(rhs_ref,pars=paste("lambda[",1:p,"]",sep=""))
lambda_RHS_data <- extract(rhs_data,pars=paste("lambda[",1:p,"]",sep=""))

tau_RHS_ref <- extract(rhs_ref,par="tau")
tau_RHS_data <- extract(rhs_data,par="tau")


k_RHS_ref <- k_RHS_ref + (lambda_RHS_ref %>% 
  map(~.*tau_RHS_ref$tau) %>%
  map(~.^2) %>%
  map(~. + 1) %>%
  map(~.^(-1)) %>%
  map_dbl(~mean(.)))

k_RHS_data <- k_RHS_data + (lambda_RHS_data %>% 
  map(~.*tau_RHS_data$tau) %>%
  map(~.^2) %>%
  map(~. + 1) %>%
  map(~.^(-1)) %>%
  map_dbl(~mean(.)))
}

plot_data <- data.frame(k_ref=k_RHS_ref/10,
                        k_data=k_RHS_data/10,
                        label=c(rep("r",k),rep("s",p-k)))

ggplot(plot_data, aes(x=k_data,y=k_ref,color=label)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",size=0.5) + 
  geom_point(size=0.5) + 
  theme_light()


## Posterior intervals: done with the last iteration

theta_ref <- extract(rhs_ref,pars=paste("theta[",1:p,"]",sep=""))
theta_data <- extract(rhs_data,pars=paste("theta[",1:p,"]",sep=""))

int_ref <- theta_ref %>% 
  map(~c(mean=mean(.),sd=sd(.))) %>%
  transpose() %>%
  map(~flatten_dbl(.)) %>%
  map_dfc(~.) %>%
  add_column(label=c(rep("r",k),rep("s",p-k)), approach=rep('ref',p))
  
int_data <- theta_data %>% 
  map(~c(mean=mean(.),sd=sd(.))) %>%
  transpose() %>%
  map(~flatten_dbl(.)) %>%
  map_dfc(~.) %>%
  add_column(label=c(rep("r",k),rep("s",p-k)), approach=rep('data',p))

plot_data <- int_ref %>% bind_rows(int_data)

post_int <- ggplot(plot_data, aes(x=rep(1:p,2),y=mean,color=label)) + 
          geom_hline(yintercept=true, linetype="dashed") +
          geom_pointrange(aes(ymin=mean-sd,ymax=mean+sd),size=0.05) + 
          facet_grid(~approach) + 
          labs(x="",y=TeX("$\\theta_{j}$")) +
          guides(color=FALSE) +
          theme_light() + theme(axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.line.x=element_blank()) 

ggsave("../paper/graphics/post_int.pdf",post_int,width=10,height=3)
