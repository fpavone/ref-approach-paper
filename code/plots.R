##########################################
### Script for the plots of the paper ####
##########################################
library(tidyverse)
theme_set(theme_light() +
          theme(strip.background = element_rect(color='black',fill='white'),
                strip.text.x = element_text(color='black'),
                strip.text.y = element_text(color='black')))
library(latex2exp)
source("getStability.R")

p <- 1000
k <- 100

avg.sensitivity <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)<=k)/k})
  return(mean(tmp, na.rm = TRUE))
}

sd.sensitivity <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)<=k)/k})
  return(sd(tmp, na.rm = TRUE))
}

avg.fdr <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)>k)/sum(x)})
  return(mean(tmp, na.rm = TRUE))
}

sd.fdr <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)>k)/sum(x)})
  return(sd(tmp, na.rm = TRUE))
}

data.plot <- tibble(n = numeric(),
                    rho = numeric(),
                    method = character(),
                    approach = character(),
                    sensitivity = numeric(),
                    sensitivity.sd = numeric(),
                    fdr = numeric(),
                    fdr.sd = numeric(),
                    stab.low = numeric(),
                    stab.mean = numeric(),
                    stab.up = numeric())

for(nn in c(50,70,100)){
  for(rr in c(0.3,0.5)){
    load(paste("ref_approach_n",nn,"rho",rr,".Rdata",sep=""))

    data.plot <- data.plot %>%
      # Credibility intervals inclusion probabilities with regularized horseshoe prior
      bind_rows(tibble(n = rep(n,2),
              rho = rep(rho,2),
              method = rep("ci.90",2),
              approach = c("ref","data"),
              sensitivity = c(avg.sensitivity(ci90_X_ref),avg.sensitivity(ci90_X_data)),
              sensitivity.sd = c(sd.sensitivity(ci90_X_ref),sd.sensitivity(ci90_X_data)),
              fdr = c(avg.fdr(ci90_X_ref),avg.fdr(ci90_X_data)),
              fdr.sd = c(sd.fdr(ci90_X_ref),sd.fdr(ci90_X_data)),
              stab.low = c(getStability(ci90_X_ref)$lower,getStability(ci90_X_data)$lower),
              stab.mean = c(getStability(ci90_X_ref)$stability,getStability(ci90_X_data)$stability),
              stab.up = c(getStability(ci90_X_ref)$upper,getStability(ci90_X_data)$upper))) %>%
      # Control of the local false discovery rate
      bind_rows(tibble(n = rep(n,2),
              rho = rep(rho,2),
              method = rep("loc.fdr",2),
              approach = c("ref","data"),
              sensitivity = c(avg.sensitivity(lfdr_X_ref),avg.sensitivity(lfdr_X_data)),
              sensitivity.sd = c(sd.sensitivity(lfdr_X_ref),sd.sensitivity(lfdr_X_data)),
              fdr = c(avg.fdr(lfdr_X_ref),avg.fdr(lfdr_X_data)),
              fdr.sd = c(sd.fdr(lfdr_X_ref),sd.fdr(lfdr_X_data)),
              stab.low = c(getStability(lfdr_X_ref)$lower,getStability(lfdr_X_data)$lower),
              stab.mean = c(getStability(lfdr_X_ref)$stability,getStability(lfdr_X_data)$stability),
              stab.up = c(getStability(lfdr_X_ref)$upper,getStability(lfdr_X_data)$upper))) %>%
      # Empirical Bayes median thresholding
      bind_rows(tibble(n = rep(n,2),
              rho = rep(rho,2),
              method = rep("EB.med",2),
              approach = c("ref","data"),
              sensitivity = c(avg.sensitivity(ebmt_X_ref),avg.sensitivity(ebmt_X_data)),
              sensitivity.sd = c(sd.sensitivity(ebmt_X_ref),sd.sensitivity(ebmt_X_data)),
              fdr = c(avg.fdr(ebmt_X_ref),avg.fdr(ebmt_X_data)),
              fdr.sd = c(sd.fdr(ebmt_X_ref),sd.fdr(ebmt_X_data)),
              stab.low = c(getStability(ebmt_X_ref)$lower,getStability(ebmt_X_data)$lower),
              stab.mean = c(getStability(ebmt_X_ref)$stability,getStability(ebmt_X_data)$stability),
              stab.up = c(getStability(ebmt_X_ref)$upper,getStability(ebmt_X_data)$upper)))

  }
}


facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
                         rho=function(x){paste("rho=",x,sep="")})

## Sensitivity vs False discovery rate plot
plot1 <- ggplot(data.plot,aes(x=fdr,y=sensitivity,col=method)) +
  facet_grid(rho~n, labeller=facet.labels) +
  geom_point(aes(shape=approach),size=2.5) +
 # geom_errorbar(aes(ymin=sensitivity-sensitivity.sd, ymax=sensitivity+sensitivity.sd)) +
 # geom_errorbarh(aes(xmin=fdr-fdr.sd, xmax=fdr+fdr.sd)) +
  scale_shape_manual(values = c(16,15)) +
  geom_line(aes(col=method)) +
  labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")
#  theme_light()

ggsave("../paper/graphics/sensitivity_vs_fdr.pdf",plot1,width=10,height=3)

## Stability plot
plot2 <- ggplot(data.plot,aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(rho~n, labeller=facet.labels) +
  geom_point(size=2.5) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))
#  theme_light()

ggsave("../paper/graphics/stability.pdf",plot2,width=10,height=3)



###########################################
########### FULL BAYES PLOTS ##############

table_results <- tibble(error = numeric(),
                        error_type = character(),
                        prior = character(),
                        approach = character(),
                        n = numeric(),
                        rho = numeric())

k_plot_RHS <- tibble(k.value = numeric(),
                 n = numeric(),
                 rho = numeric(),
                 label = character(),
                 approach = character())

k_plot_DL <- tibble(k.value = numeric(),
                     n = numeric(),
                     rho = numeric(),
                     label = character(),
                     approach = character())

for(nn in c(50,70,100)){
  for(rr in c(0.3,0.5)){
    load(paste("fullBayes_n",nn,"rho",rr,".Rdata",sep=""))

    # Correcting SSE and SESE to be independent of n
    #results$error <- results$error/(n-3)

    table_results <- table_results %>%
      bind_rows(results)

    k_plot_RHS <- k_plot_RHS %>%
      bind_rows((k_result_RHS %>%
                  add_column(id = rep(1:p,2))))

    k_plot_DL <- k_plot_DL %>%
      bind_rows((k_result_DL %>%
                   add_column(id = rep(1:p,2))))
  }
}

facet.labels <- labeller(error_type = function(x){x},
                         rho=function(x){paste("rho=",x,sep="")})

plot_SESE_SSE <- table_results %>%
  group_by(error_type,prior,approach,n,rho) %>%
  summarise(avg.error=mean(error), sd.error = sd(error)) %>%
  ggplot(aes(x=n, y=avg.error, color=approach)) +
  geom_line(aes(linetype=prior)) +
  geom_errorbar(aes(ymin=avg.error-sd.error,ymax=avg.error+sd.error), width = 1) +
  geom_point() +
  scale_x_continuous(breaks=c(50,70,100)) +
  facet_grid(error_type~rho, labeller=facet.labels) +
  labs(x="Number of observations", y="", linetype="Prior", color="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))
 # theme_light()

ggsave("../paper/graphics/SESE_SSE.pdf",plot_SESE_SSE,width=10,height=3)



facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
                         rho=function(x){paste("rho=",x,sep="")})

plot_k_RHS <- k_plot_RHS %>%
  spread(key=approach,value=k.value) %>%
  ggplot(aes(x=data,y=ref,color=label)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",size=0.5) +
  geom_point(size=0.5) +
  facet_grid(rho~n, labeller=facet.labels) +
  labs(x="kappa.data",y="kappa.ref") +
  guides(color=FALSE)
#  theme_light()

ggsave("../paper/graphics/k_RHS.pdf",plot_k_RHS,width=10,height=3)

plot_k_DL <- k_plot_DL %>%
  spread(key=approach,value=k.value) %>%
  ggplot(aes(x=data,y=ref,color=label)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",size=0.5) +
  geom_point(size=0.5) +
  facet_grid(rho~n, labeller=facet.labels) +
  labs(x="kappa.data",y="kappa.ref") +
  guides(color=FALSE)
#  theme_light()

ggsave("../paper/graphics/k_DL.pdf",plot_k_DL,width=10,height=3)




## Posterior intervals: done with the last iteration
load(paste("fullBayes_n50rho0.3.Rdata",sep=""))

theta_ref <- rstan::extract(rhs_ref,pars=paste("theta[",1:p,"]",sep=""))
theta_data <- rstan::extract(rhs_data,pars=paste("theta[",1:p,"]",sep=""))

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
    labs(x="",y="") +
    ##labs(x="",y=TeX("$\\theta_{j}$")) +
    guides(color=FALSE) +
  #  theme_light() +
    theme(axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.line.x=element_blank())

ggsave("../paper/graphics/post_int.pdf",post_int,width=10,height=3)



#####################################################
###### BOOTSTRAP INCLUSION PROBABILITES PLOT ########

load("bodyfat_notebook.RData") # originally bodyfat_bootstrap.Rdata
ordered <- boot_inclusion %>%
    arrange(desc(projpred))
inc_prob <- boot_inclusion %>%
  gather(key="method", value="freq", projpred,steplm) %>%
  ggplot(aes(x=variable,y=freq,fill=method)) +
  geom_bar(stat="identity",position=position_dodge()) +
  labs(x="",y="",fill="") +
  #guides(fill=FALSE) +
  scale_x_discrete(limits=ordered$variable) +
  scale_y_continuous(labels = scales::number_format(suffix="%")) +
  scale_fill_manual(values=c('steplm'="#819FF7",'projpred'="#FAAC58"))
#  theme_light()

ggsave("../paper/graphics/inc_prob.pdf",inc_prob,width=10,height=3)


# Model selection frequencies
models_step <- as_tibble(boot_step) %>%
    group_by_at(pred) %>%
    count() %>%
    ungroup() %>%
    arrange(desc(n))

for(i in 1:min(10,nrow(models_step))){
    print(paste(paste(pred[models_step[i,pred]==1],collapse=' + '),
                models_step$n[i]/bootnum*100))
}

models_proj <- as_tibble(boot_proj) %>%
    group_by_at(pred) %>%
    count() %>%
    ungroup() %>%
    arrange(desc(n))

for(i in 1:min(10,nrow(models_proj))){
    print(paste(paste(pred[models_proj[i,pred]==1],collapse=' + '),
                models_proj$n[i]/bootnum*100))
}


#####################################################
################ CORRELATION PLOT  ##################
library(rstan)
library(dimreduce)

## Simulation parameters
n <- 70  # number of observations
rho <- 0.3  # correlation level
p <- 1000   # total number of features
k <- 100    # number of relevant features

## Data simulation mechanism
f <- rnorm(n) # the true underlying function
y <- f + rnorm(n)
x <- as_tibble(sqrt(1-rho)*matrix(rnorm(n*k), ncol=k) + sqrt(rho)*f) %>%  # set of relevant covariates
    bind_cols(as_tibble(matrix(rnorm(n*(p-k)), ncol=p-k))) %>%  # set of spurious covariates
    set_names(c(paste("r.",1:k,sep=""),paste("s.",(k+1):p,sep="")))

nc <- 5
dr <- spca(x,y, nctot=nc, screenthresh = 0.6, window=p, sup.only = T, verbose = F)
z <- predict(dr, x)

model_data <- list(y = y,
                   X = z,
                   p = nc,
                   n = n,
                   s_max = dr$sdev[1])
fit <- stan("model.stan", data = model_data,
            chains=1, seed=45342, save_dso=F)

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

cor_plot <- tibble(corXf = abs(cor(x,f)),
                   corXYfit = abs(cor(x,yfit)),
                   corXY = abs(cor(x,y)),
                   label = c(rep('r',k),rep('s',p-k))) %>%
  ggplot(aes(x=corXY,y=corXYfit,color=label)) +
  geom_point(size=0.3) +
  labs(x=TeX('$|Cor(x,y)|$'),y=TeX('$|Cor(x,\\hat{y})|$')) +
  guides(color=FALSE) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1))
#  theme_light()


cor_plot_final <- ggExtra::ggMarginal(cor_plot,
                                      groupColour = TRUE,
                                      groupFill = TRUE,
                                     # size = 8,
                                      type = 'density',
                                      xparams = list(size=0.5),
                                      yparams = list(size=0.5))
ggsave("../paper/graphics/correlation.pdf",cor_plot_final, width=4, height=4)


##################################################
########## BODYFAT VARIABLE SELECTION ############

k <- 13

avg.sensitivity <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)<=k)/k})
  return(mean(tmp, na.rm = TRUE))
}

avg.fdr <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)>k)/sum(x)})
  return(mean(tmp, na.rm = TRUE))
}

data.plot <- tibble(n = numeric(),
                    rho = numeric(),
                    method = character(),
                    approach = character(),
                    recall = numeric(),
                    fdr = numeric(),
                    stab.low = numeric(),
                    stab.mean = numeric(),
                    stab.up = numeric())

for(nn in c(50,70,100,251)){
  load(paste("bodyfat_type1_varyN_n",nn,".Rdata",sep=""))

  data.plot <- data.plot %>%
    # Credibility intervals inclusion probabilities with regularized horseshoe prior
    bind_rows(tibble(n = rep(n,2),
                     method = rep("ci.90",2),
                     approach = c("ref","data"),
                     sensitivity = c(avg.sensitivity(ci90_X_ref),avg.sensitivity(ci90_X_data)),
                     fdr = c(avg.fdr(ci90_X_ref),avg.fdr(ci90_X_data)),
                     stab.low = c(getStability(ci90_X_ref)$lower,getStability(ci90_X_data)$lower),
                     stab.mean = c(getStability(ci90_X_ref)$stability,getStability(ci90_X_data)$stability),
                     stab.up = c(getStability(ci90_X_ref)$upper,getStability(ci90_X_data)$upper))) %>%
    # Control of the local false discovery rate
    bind_rows(tibble(n = rep(n,2),
                     method = rep("loc.fdr",2),
                     approach = c("ref","data"),
                     sensitivity = c(avg.sensitivity(lfdr_X_ref),avg.sensitivity(lfdr_X_data)),
                     fdr = c(avg.fdr(lfdr_X_ref),avg.fdr(lfdr_X_data)),
                     stab.low = c(getStability(lfdr_X_ref)$lower,getStability(lfdr_X_data)$lower),
                     stab.mean = c(getStability(lfdr_X_ref)$stability,getStability(lfdr_X_data)$stability),
                     stab.up = c(getStability(lfdr_X_ref)$upper,getStability(lfdr_X_data)$upper))) %>%
    # Empirical Bayes median thresholding
    bind_rows(tibble(n = rep(n,2),
                     method = rep("EB.med",2),
                     approach = c("ref","data"),
                     sensitivity = c(avg.sensitivity(ebmt_X_ref),avg.sensitivity(ebmt_X_data)),
                     fdr = c(avg.fdr(ebmt_X_ref),avg.fdr(ebmt_X_data)),
                     stab.low = c(getStability(ebmt_X_ref)$lower,getStability(ebmt_X_data)$lower),
                     stab.mean = c(getStability(ebmt_X_ref)$stability,getStability(ebmt_X_data)$stability),
                     stab.up = c(getStability(ebmt_X_ref)$upper,getStability(ebmt_X_data)$upper)))
}



facet.labels <- labeller(n = function(x){paste("n=",x,sep="")})

## Sensitivity vs False discovery rate plot
plot1 <- ggplot(data.plot,aes(x=fdr,y=sensitivity,col=method)) +
    facet_grid(~n, labeller=facet.labels) +
    ##geom_abline(intercept=0, slope=1, linetype='dashed') +
    scale_x_continuous(limits=c(0,0.6)) +
    scale_y_continuous(limits=c(0.3,0.9)) +
    geom_point(aes(shape=approach),size=2.5) +
    scale_shape_manual(values = c(16,15)) +
    geom_line(aes(col=method)) +
    labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")
   # theme_light() +
   # theme(strip.background = element_rect(color='black',fill='white'),
  #        strip.text.x = element_text(color='black'),
   #       strip.text.y = element_text(color='black'))


ggsave("../paper/graphics/bodyfat_sensitivity_vs_fdr.pdf",plot1,width=10,height=2.1)

## Stability plot
plot2 <- ggplot(data.plot,aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(~n, labeller=facet.labels) +
  geom_point(size=2.5) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))
#  theme_light()

ggsave("../paper/graphics/bodyfat_stability.pdf",plot2,width=10,height=2)


############################################################
########## BODYFAT VARIABLE SELECTION VARYING P ############

k <- 13

avg.sensitivity <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)<=k)/k})
  return(mean(tmp, na.rm = TRUE))
}

avg.fdr <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)>k)/sum(x)})
  return(mean(tmp, na.rm = TRUE))
}

data.plot <- tibble(n = numeric(),
                    rho = numeric(),
                    method = character(),
                    approach = character(),
                    recall = numeric(),
                    fdr = numeric(),
                    stab.low = numeric(),
                    stab.mean = numeric(),
                    stab.up = numeric())

for(pp in c(100,500,1000,2000,5000)){
  load(paste("bodyfat_type1_varyP_p",pp,".Rdata",sep=""))

  data.plot <- data.plot %>%
    # Credibility intervals inclusion probabilities with regularized horseshoe prior
    bind_rows(tibble(p = rep(p,2),
                     method = rep("ci.90",2),
                     approach = c("ref","data"),
                     sensitivity = c(avg.sensitivity(ci90_X_ref),avg.sensitivity(ci90_X_data)),
                     fdr = c(avg.fdr(ci90_X_ref),avg.fdr(ci90_X_data)),
                     stab.low = c(getStability(ci90_X_ref)$lower,getStability(ci90_X_data)$lower),
                     stab.mean = c(getStability(ci90_X_ref)$stability,getStability(ci90_X_data)$stability),
                     stab.up = c(getStability(ci90_X_ref)$upper,getStability(ci90_X_data)$upper))) %>%
    # Control of the local false discovery rate
    bind_rows(tibble(p = rep(p,2),
                     method = rep("loc.fdr",2),
                     approach = c("ref","data"),
                     sensitivity = c(avg.sensitivity(lfdr_X_ref),avg.sensitivity(lfdr_X_data)),
                     fdr = c(avg.fdr(lfdr_X_ref),avg.fdr(lfdr_X_data)),
                     stab.low = c(getStability(lfdr_X_ref)$lower,getStability(lfdr_X_data)$lower),
                     stab.mean = c(getStability(lfdr_X_ref)$stability,getStability(lfdr_X_data)$stability),
                     stab.up = c(getStability(lfdr_X_ref)$upper,getStability(lfdr_X_data)$upper))) %>%
    # Empirical Bayes median thresholding
    bind_rows(tibble(p = rep(p,2),
                     method = rep("EB.med",2),
                     approach = c("ref","data"),
                     sensitivity = c(avg.sensitivity(ebmt_X_ref),avg.sensitivity(ebmt_X_data)),
                     fdr = c(avg.fdr(ebmt_X_ref),avg.fdr(ebmt_X_data)),
                     stab.low = c(getStability(ebmt_X_ref)$lower,getStability(ebmt_X_data)$lower),
                     stab.mean = c(getStability(ebmt_X_ref)$stability,getStability(ebmt_X_data)$stability),
                     stab.up = c(getStability(ebmt_X_ref)$upper,getStability(ebmt_X_data)$upper)))
}



facet.labels <- labeller(p = function(x){paste("p=",x,sep="")})

## Sensitivity vs False discovery rate plot
plot1 <- ggplot(data.plot,aes(x=fdr,y=sensitivity,col=method)) +
  facet_grid(~p, labeller=facet.labels) +
  #geom_abline(intercept=0, slope=1, linetype='dashed') +
  #scale_x_continuous(limits=c(0,0.6)) +
  #scale_y_continuous(limits=c(0.3,0.9)) +
  geom_point(aes(shape=approach),size=2.5) +
  scale_shape_manual(values = c(16,15)) +
  geom_line(aes(col=method)) +
  labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")
#  theme_light()

ggsave("../paper/graphics/bodyfat_varyingP_sensitivity_vs_fdr.pdf",plot1,width=12,height=2.1)

## Stability plot
plot2 <- ggplot(data.plot,aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(~p, labeller=facet.labels) +
  geom_point(size=2.5) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))
#  theme_light()

ggsave("../paper/graphics/bodyfat_varyingP_stability.pdf",plot2,width=12,height=2)

#################################################
###### BODYFAT STEPWISE w/wout REF APPROACH #####

load('bodyfat_step.Rdata')

data.plot <- data.frame(stat=c(noisy_ref,noisy_data,rmse_ref,rmse_data),
                        method=rep(c(rep('ref',100),rep('data',100)),2),
                        type=c(rep('noisy',200),rep('rmse',200)))

facet.labels <- labeller(type=function(x){ifelse(x=='noisy','Noisy features selected','RMSE')})

plot1 <- ggplot(data.plot,aes(x=stat,fill=method,color=method)) +
    facet_grid(~type, scales='free', labeller=facet.labels) +
    geom_histogram(position = "identity", alpha=0.7, bins=30) +
    labs(x='',y='',fill='Approach') +
    guides(color=FALSE) +
    scale_fill_manual(values=c("#819FF7","#FAAC58")) +
    scale_color_manual(values=c("#819FF7","#FAAC58"))
#    theme_light()

ggsave("../paper/graphics/bodyfat_step_refvsdata.pdf",plot1,width=10,height=2)
