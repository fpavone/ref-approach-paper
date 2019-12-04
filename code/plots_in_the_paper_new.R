##########################################
### Script for the plots of the paper ####
##########################################
library(tidyverse)
library(ggrepel)
theme_set(theme_light() +
          theme(strip.background = element_rect(color='black',fill='white'),
                strip.text.x = element_text(color='black'),
                strip.text.y = element_text(color='black')))
shapes <- c('ref'=8, 'data'=20)
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

colors <- c(
  'projpred.iter' = '#2E2E2E',
  'loc.fdr' = '#619CFF',
  'ci.90' = '#F8766D',
  'EB.med' = '#00BA38'
)

## If you want to run the script to create the data object
## set TRUE, otherwise set FALSE to load the data.
saveMode <- FALSE

######################################################
###### COMPLETE SELECTION (OLD VERSION) ##############
######################################################

if(saveMode){
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
                ## Credibility intervals inclusion probabilities
                ## with regularized horseshoe prior
                bind_rows(tibble(n = rep(n,2),
                                 rho = rep(rho,2),
                                 method = rep("ci.90",2),
                                 approach = c("ref","data"),
                                 sensitivity = c(avg.sensitivity(ci90_X_ref),
                                                 avg.sensitivity(ci90_X_data)),
                                 sensitivity.sd = c(sd.sensitivity(ci90_X_ref),
                                                    sd.sensitivity(ci90_X_data)),
                                 fdr = c(avg.fdr(ci90_X_ref),
                                         avg.fdr(ci90_X_data)),
                                 fdr.sd = c(sd.fdr(ci90_X_ref),
                                            sd.fdr(ci90_X_data)),
                                 stab.low = c(getStability(ci90_X_ref)$lower,
                                              getStability(ci90_X_data)$lower),
                                 stab.mean = c(getStability(ci90_X_ref)$stability,
                                               getStability(ci90_X_data)$stability),
                                 stab.up = c(getStability(ci90_X_ref)$upper,
                                             getStability(ci90_X_data)$upper))) %>%
                ## Control of the local false discovery rate
                bind_rows(tibble(n = rep(n,2),
                                 rho = rep(rho,2),
                                 method = rep("loc.fdr",2),
                                 approach = c("ref","data"),
                                 sensitivity = c(avg.sensitivity(lfdr_X_ref),
                                                 avg.sensitivity(lfdr_X_data)),
                                 sensitivity.sd = c(sd.sensitivity(lfdr_X_ref),
                                                    sd.sensitivity(lfdr_X_data)),
                                 fdr = c(avg.fdr(lfdr_X_ref),
                                         avg.fdr(lfdr_X_data)),
                                 fdr.sd = c(sd.fdr(lfdr_X_ref),
                                            sd.fdr(lfdr_X_data)),
                                 stab.low = c(getStability(lfdr_X_ref)$lower,
                                              getStability(lfdr_X_data)$lower),
                                 stab.mean = c(getStability(lfdr_X_ref)$stability,
                                               getStability(lfdr_X_data)$stability),
                                 stab.up = c(getStability(lfdr_X_ref)$upper,
                                             getStability(lfdr_X_data)$upper))) %>%
                ## Empirical Bayes median thresholding
                bind_rows(tibble(n = rep(n,2),
                                 rho = rep(rho,2),
                                 method = rep("EB.med",2),
                                 approach = c("ref","data"),
                                 sensitivity = c(avg.sensitivity(ebmt_X_ref),
                                                 avg.sensitivity(ebmt_X_data)),
                                 sensitivity.sd = c(sd.sensitivity(ebmt_X_ref),
                                                    sd.sensitivity(ebmt_X_data)),
                                 fdr = c(avg.fdr(ebmt_X_ref),
                                         avg.fdr(ebmt_X_data)),
                                 fdr.sd = c(sd.fdr(ebmt_X_ref),
                                            sd.fdr(ebmt_X_data)),
                                 stab.low = c(getStability(ebmt_X_ref)$lower,
                                              getStability(ebmt_X_data)$lower),
                                 stab.mean = c(getStability(ebmt_X_ref)$stability,
                                               getStability(ebmt_X_data)$stability),
                                 stab.up = c(getStability(ebmt_X_ref)$upper,
                                             getStability(ebmt_X_data)$upper)))
        }
    }
    save(data.plot, file='complete_selection_old_plot.RData')
} else
    load('complete_selection_old_plot.RData')

data.plot <- data.plot %>%
  mutate(method = factor(method, levels = names(colors)))

facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
                         rho=function(x){paste("rho=",x,sep="")})


## Sensitivity vs False discovery rate plot
plot1 <- ggplot(data.plot,aes(x=fdr,y=sensitivity,col=method)) +
  facet_grid(rho~n, labeller=facet.labels) +
  geom_point(aes(shape=approach),size=2.5) +
 # geom_errorbar(aes(ymin=sensitivity-sensitivity.sd, ymax=sensitivity+sensitivity.sd)) +
 # geom_errorbarh(aes(xmin=fdr-fdr.sd, xmax=fdr+fdr.sd)) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  geom_line(aes(col=method)) +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 2)
  ) +
  labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")
#  theme_light()

ggsave("graphics/sensitivity_vs_fdr.pdf",plot1,width=10,height=3)

## Stability plot
plot2 <- ggplot(data.plot,aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(rho~n, labeller=facet.labels) +
  geom_point(size=2.5) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))
#  theme_light()

ggsave("graphics/stability.pdf",plot2,width=10,height=3)






#####################################################
###### BOOTSTRAP INCLUSION PROBABILITES PLOT ########
#####################################################

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

ggsave("graphics/inc_prob.pdf",inc_prob,width=10,height=3)


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
#####################################################

library(rstan)
# remotes::install_github("jpiironen/dimreduce")
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
  scale_y_continuous(limits = c(0,1)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12)
  )
#  theme_light()

cor_plot


cor_plot_final <- ggExtra::ggMarginal(cor_plot,
                                      groupColour = TRUE,
                                      groupFill = TRUE,
                                     # size = 8,
                                      type = 'density',
                                      xparams = list(size=0.5),
                                      yparams = list(size=0.5))
ggsave("graphics/correlation.pdf",cor_plot_final, width=4, height=4)


##################################################
########## BODYFAT VARIABLE SELECTION ############
##################################################
k <- 13

avg.sensitivity <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)<=k)/k})
  return(mean(tmp, na.rm = TRUE))
}

avg.fdr <- function(X){
  tmp <- apply(X,1,function(x){sum(which(x==1)>k)/sum(x)})
  return(mean(tmp, na.rm = TRUE))
}

if(saveMode){
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
            ## Credibility intervals inclusion probabilities
            ## with regularized horseshoe prior
            bind_rows(tibble(n = rep(n,2),
                             method = rep("ci.90",2),
                             approach = c("ref","data"),
                             sensitivity = c(avg.sensitivity(ci90_X_ref),
                                             avg.sensitivity(ci90_X_data)),
                             fdr = c(avg.fdr(ci90_X_ref),
                                     avg.fdr(ci90_X_data)),
                             stab.low = c(getStability(ci90_X_ref)$lower,
                                          getStability(ci90_X_data)$lower),
                             stab.mean = c(getStability(ci90_X_ref)$stability,
                                           getStability(ci90_X_data)$stability),
                             stab.up = c(getStability(ci90_X_ref)$upper,
                                         getStability(ci90_X_data)$upper))) %>%
            ## Control of the local false discovery rate
            bind_rows(tibble(n = rep(n,2),
                             method = rep("loc.fdr",2),
                             approach = c("ref","data"),
                             sensitivity = c(avg.sensitivity(lfdr_X_ref),
                                             avg.sensitivity(lfdr_X_data)),
                             fdr = c(avg.fdr(lfdr_X_ref),
                                     avg.fdr(lfdr_X_data)),
                             stab.low = c(getStability(lfdr_X_ref)$lower,
                                          getStability(lfdr_X_data)$lower),
                             stab.mean = c(getStability(lfdr_X_ref)$stability,
                                           getStability(lfdr_X_data)$stability),
                             stab.up = c(getStability(lfdr_X_ref)$upper,
                                         getStability(lfdr_X_data)$upper))) %>%
            ## Empirical Bayes median thresholding
            bind_rows(tibble(n = rep(n,2),
                             method = rep("EB.med",2),
                             approach = c("ref","data"),
                             sensitivity = c(avg.sensitivity(ebmt_X_ref),
                                             avg.sensitivity(ebmt_X_data)),
                             fdr = c(avg.fdr(ebmt_X_ref),avg.fdr(ebmt_X_data)),
                             stab.low = c(getStability(ebmt_X_ref)$lower,
                                          getStability(ebmt_X_data)$lower),
                             stab.mean = c(getStability(ebmt_X_ref)$stability,
                                           getStability(ebmt_X_data)$stability),
                             stab.up = c(getStability(ebmt_X_ref)$upper,
                                         getStability(ebmt_X_data)$upper)))
        ## Iterative projection
        load(paste("bodyfat_type1_varyN_iteratedproj_n",nn,".Rdata",sep=""))
        data.plot <- data.plot %>%
            bind_rows(tibble(n = n,
                             method = "projpred.iter",
                             approach = "ref",
                             sensitivity =  avg.sensitivity(projpred_X),
                             fdr = avg.fdr(projpred_X),
                             stab.low = getStability(projpred_X)$lower,
                             stab.mean = getStability(projpred_X)$stability,
                             stab.up = getStability(projpred_X)$upper))
    }
    save(data.plot, file='bodyfat_complete_selection_plot.RData')
} else
    load('bodyfat_complete_selection_plot.RData')

## colors <- c('loc.fdr' = '#819FF7',
##             'ci.90' = '#F78181',
##             'EB.med' = '#298A08')

facet.labels <- labeller(n = function(x){paste("n=",x,sep="")})

## Sensitivity vs False discovery rate plot
plot1 <- data.plot %>%
    mutate(method = factor(method, levels = names(colors))) %>%
    ggplot(aes(x=fdr,y=sensitivity,col=method)) +
    facet_grid(~n, labeller=facet.labels) +
    ##geom_abline(intercept=0, slope=1, linetype='dashed') +
    scale_x_continuous(limits=c(0,0.6)) +
    scale_y_continuous(limits=c(0.3,0.99)) +
    geom_point(aes(shape=approach), size=2) +
    scale_shape_manual(values = shapes) +
    scale_color_manual(values = colors) +
    geom_line(aes(col=method)) +
    guides(
      shape = guide_legend(order = 1),
      colour = guide_legend(order = 2)
    ) +
    labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")
   # theme_light() +
   # theme(strip.background = element_rect(color='black',fill='white'),
  #        strip.text.x = element_text(color='black'),
   #       strip.text.y = element_text(color='black'))


ggsave("graphics/bodyfat_sensitivity_vs_fdr.pdf",plot1,width=10,height=2.1)

## Stability plot
plot2 <- data.plot %>%
  mutate(method = factor(method, levels = rev(names(colors)))) %>%
  ggplot(aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(~n, labeller=facet.labels) +
  geom_point(size=3, alpha = 0.8) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))
#  theme_light()

ggsave("graphics/bodyfat_stability.pdf",plot2,width=10,height=2)



#################################################
###### BODYFAT STEPWISE w/wout REF APPROACH #####
#################################################


# TODO: add new data to the plot
# load('bodyfat_step.Rdata')
#
# data.plot <- data.frame(stat=c(noisy_ref,noisy_data,rmse_ref,rmse_data),
#                         method=rep(c(rep('ref',100),rep('data',100)),2),
#                         type=c(rep('noisy',200),rep('rmse',200)))
#
# facet.labels <- labeller(type=function(x){ifelse(x=='noisy','Noisy features selected','RMSE')})
#
# plot1 <- ggplot(data.plot,aes(x=stat,fill=method,color=method)) +
#     facet_grid(~type, scales='free', labeller=facet.labels) +
#     geom_histogram(position = "identity", alpha=0.7, bins=30) +
#     labs(x='',y='',fill='Approach') +
#     guides(color=FALSE) +
#     scale_fill_manual(values=c("#819FF7","#FAAC58")) +
#     scale_color_manual(values=c("#819FF7","#FAAC58")) +
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank())
# #    theme_light()
#
# ggsave("graphics/bodyfat_step_refvsdata.pdf",plot1,width=10,height=2)



##++++++++++++++++++++++++++++++++++++++++++++++++#
##++++++++++++++++++++++++++++++++++++++++++++++++#
####### MINIMAL SUBSET SELECTION PARALLEL #########
##++++++++++++++++++++++++++++++++++++++++++++++++#
##++++++++++++++++++++++++++++++++++++++++++++++++#
library(entropy)

if(saveMode){
    data.plot <- tibble(n = numeric(),
                        rho = numeric(),
                        method = character(),
                        fdr = numeric(),
                        rmse = numeric(),
                        entropy = numeric())

    intervals <- matrix(1:100,nrow=25,ncol=4)

    for(nn in c(80,100,150)){
        for(rr in c(0.3,0.5)){
            X.projpred.tot <- numeric()
            X.step.data.tot <- numeric()
            X.step.ref.tot <- numeric()
            X.bayes.step.data.tot <- numeric()
            X.bayes.step.ref.tot <- numeric()
            rmse.projpred.tot <- numeric()
            rmse.step.data.tot <- numeric()
            rmse.step.ref.tot <- numeric()
            rmse.bayes.step.data.tot <- numeric()
            rmse.bayes.step.ref.tot <- numeric()
            for(ii in seq(1,25)){
                load(paste('minimal_subset_parallel/minimal_subset_final_n',
                           nn,'_rho',rr,'_inter',ii,'.RData',sep=''))
                which.intervals <- intervals[ii,]
                X.projpred.tot <- rbind(X.projpred.tot,
                                        X.projpred[which.intervals,])
                X.step.data.tot <- rbind(X.step.data.tot,
                                         X.step.data[which.intervals,])
                X.step.ref.tot <- rbind(X.step.ref.tot,
                                        X.step.ref[which.intervals,])
                X.bayes.step.data.tot <- rbind(X.bayes.step.data.tot,
                                               X.bayes.step.data[which.intervals,])
                X.bayes.step.ref.tot <- rbind(X.bayes.step.ref.tot,
                                              X.bayes.step.ref[which.intervals,])
                rmse.projpred.tot <- c(rmse.projpred.tot,
                                       rmse.projpred[which.intervals])
                rmse.step.data.tot <- c(rmse.step.data.tot,
                                        rmse.step.data[which.intervals])
                rmse.step.ref.tot <- c(rmse.step.ref.tot,
                                       rmse.step.ref[which.intervals])
                rmse.bayes.step.data.tot <- c(rmse.bayes.step.data.tot,
                                              rmse.bayes.step.data[which.intervals])
                rmse.bayes.step.ref.tot <- c(rmse.bayes.step.ref.tot,
                                             rmse.bayes.step.ref[which.intervals])
            }
            data.plot <- data.plot %>%
                bind_rows(tibble(n = rep(n,5),
                                 rho = rep(rho,5),
                                 approach = c('ref',
                                              'data',
                                              'ref',
                                              'data',
                                              'ref'),
                                 method = c('projpred',
                                            'step.lm',
                                            'step.lm',
                                            'step.bayes',
                                            'step.bayes'),
                                 fdr = c(avg.fdr(X.projpred.tot),
                                         avg.fdr(X.step.data.tot),
                                         avg.fdr(X.step.ref.tot),
                                         avg.fdr(X.bayes.step.data.tot),
                                         avg.fdr(X.bayes.step.ref.tot)),
                                 rmse = c(mean(rmse.projpred.tot),
                                          mean(rmse.step.data.tot),
                                          mean(rmse.step.ref.tot),
                                          mean(rmse.bayes.step.data.tot),
                                          mean(rmse.bayes.step.ref.tot)),
                                 entropy = c(entropy(y=apply(X.projpred.tot,2,sum)),
                                             entropy(apply(X.step.data.tot,2,sum)),
                                             entropy(apply(X.step.ref.tot,2,sum)),
                                             entropy(apply(X.bayes.step.data.tot,2,
                                                           sum)),
                                             entropy(apply(X.bayes.step.ref.tot,2,
                                                           sum)))))
        }
    }
    save(data.plot,file='minimal_subset_selection_parallel_plot.RData')
} else
    load('minimal_subset_selection_parallel_plot.RData')

facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
                         rho=function(x){paste("rho=",x,sep="")})

data.plot <- data.plot %>%
  mutate(method = factor(method, levels = c("projpred", "step.bayes", "step.lm")))

plot1 <- ggplot(data.plot, aes(x=fdr,y=rmse, color=method)) +
    facet_grid(rho~n, labeller=facet.labels) +
    geom_point(aes(shape=approach),size=3) +
    scale_shape_manual(values = shapes) +
    geom_line(data=filter(data.plot,method=='step.lm'),aes(x=fdr,y=rmse)) +
    geom_line(data=filter(data.plot,method=='step.bayes'),aes(x=fdr,y=rmse)) +
    guides(
      shape = guide_legend(order = 1),
      colour = guide_legend(order = 2)
    ) +
    labs(x='False discovery rate', y='RMSE',color='Method', shape='Approach')

ggsave('graphics/rmse_vs_fdr_parallel.pdf',plot1,width=10,height=3)

plot2 <- data.plot %>%
    mutate(method = factor(method, levels = c("step.lm", "step.bayes", "projpred"))) %>%
    ggplot(aes(x=entropy,y=method, color=approach)) +
    facet_grid(rho~n, labeller = facet.labels) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values=c("#819FF7","#FAAC58")) +
    labs(x='Entropy',y='',color='Approach')

ggsave('graphics/entropy_parallel.pdf',plot2,width=10,height=3)


################################################################
####### COMPLETE SELECTION WITH ITERATED PROJ ##################
################################################################

if(saveMode){
    data.plot <- tibble(n = numeric(),
                        rho = numeric(),
                        alpha = numeric(),
                        method = character(),
                        approach = character(),
                        sensitivity = numeric(),
                        ## sensitivity.sd = numeric(),
                        fdr = numeric(),
                        ## fdr.sd = numeric(),
                        stab.low = numeric(),
                        stab.mean = numeric(),
                        stab.up = numeric())

    for(nn in c(50,70,100)){
        for(rr in c(0.3,0.5)){
            load(paste('iterated_n',
                       nn,'_rho',rr,'_3alphas_suponlyT.Rdata',sep = ''))
            aa <- 2 ## correspond to alpha = 0.16, default value ## for(aa in 1:length(alpha_vec)){
            X_ref_test <- X_proj_test_list[[aa]]
            X_ref_cv <- X_proj_cv_list[[aa]]
            X_data_test <- X_lasso_test_list[[aa]]
            X_data_cv <- X_lasso_cv_list[[aa]]
            ## Stability
            stab_ref_test <- getStability(X_ref_test)
            stab_ref_cv <- getStability(X_ref_cv)
            stab_data_test <- getStability(X_data_test)
            stab_data_cv <- getStability(X_data_cv)
            ## Save the data
            data.plot <- data.plot %>%
                bind_rows(tibble(n = rep(n,4),
                                 rho = rep(rho,4),
                                 alpha = rep(alpha_vec[aa],4),
                                 method = c('test.iter',
                                            'test.iter',
                                            'projpred.iter',
                                            'projpred.iter'),
                                 approach = rep(c('ref','data'),2),
                                 sensitivity = c(avg.sensitivity(X_ref_test),
                                                 avg.sensitivity(X_data_test),
                                                 avg.sensitivity(X_ref_cv),
                                                 avg.sensitivity(X_data_cv)),
                                 fdr = c(avg.fdr(X_ref_test),
                                         avg.fdr(X_data_test),
                                         avg.fdr(X_ref_cv),
                                         avg.fdr(X_data_cv)),
                                 stab.low = c(stab_ref_test$lower,
                                              stab_data_test$lower,
                                              stab_ref_cv$lower,
                                              stab_data_cv$lower),
                                 stab.mean = c(stab_ref_test$stability,
                                               stab_data_test$stability,
                                               stab_ref_cv$stability,
                                               stab_data_cv$stability),
                                 stab.up = c(stab_ref_test$upper,
                                             stab_data_test$upper,
                                             stab_ref_cv$upper,
                                             stab_data_cv$upper)))
            ## }
            load(paste("ref_approach_n",
                       nn,"rho",rr,".Rdata",sep=""))
            ## Save the data
            data.plot <- data.plot %>%
                ## Credibility intervals inclusion probabilities
                ## with regularized horseshoe prior
                bind_rows(tibble(n = rep(n,2),
                                 rho = rep(rho,2),
                                 alpha = rep(NA,2),
                                 method = rep("ci.90",2),
                                 approach = c("ref","data"),
                                 sensitivity = c(avg.sensitivity(ci90_X_ref),
                                                 avg.sensitivity(ci90_X_data)),
                                 fdr = c(avg.fdr(ci90_X_ref),
                                         avg.fdr(ci90_X_data)),
                                 stab.low = c(getStability(ci90_X_ref)$lower,
                                              getStability(ci90_X_data)$lower),
                                 stab.mean = c(getStability(ci90_X_ref)$stability,
                                               getStability(ci90_X_data)$stability),
                                 stab.up = c(getStability(ci90_X_ref)$upper,
                                             getStability(ci90_X_data)$upper))) %>%
                ## Control of the local false discovery rate
                bind_rows(tibble(n = rep(n,2),
                                 rho = rep(rho,2),
                                 alpha = rep(NA,2),
                                 method = rep("loc.fdr",2),
                                 approach = c("ref","data"),
                                 sensitivity = c(avg.sensitivity(lfdr_X_ref),
                                                 avg.sensitivity(lfdr_X_data)),
                                 fdr = c(avg.fdr(lfdr_X_ref),
                                         avg.fdr(lfdr_X_data)),
                                 stab.low = c(getStability(lfdr_X_ref)$lower,
                                              getStability(lfdr_X_data)$lower),
                                 stab.mean = c(getStability(lfdr_X_ref)$stability,
                                               getStability(lfdr_X_data)$stability),
                                 stab.up = c(getStability(lfdr_X_ref)$upper,
                                             getStability(lfdr_X_data)$upper))) %>%
                ## Empirical Bayes median thresholding
                bind_rows(tibble(n = rep(n,2),
                                 rho = rep(rho,2),
                                 alpha = rep(NA,2),
                                 method = rep("EB.med",2),
                                 approach = c("ref","data"),
                                 sensitivity = c(avg.sensitivity(ebmt_X_ref),
                                                 avg.sensitivity(ebmt_X_data)),
                                 fdr = c(avg.fdr(ebmt_X_ref),
                                         avg.fdr(ebmt_X_data)),
                                 stab.low = c(getStability(ebmt_X_ref)$lower,
                                              getStability(ebmt_X_data)$lower),
                                 stab.mean = c(getStability(ebmt_X_ref)$stability,
                                               getStability(ebmt_X_data)$stability),
                                 stab.up = c(getStability(ebmt_X_ref)$upper,
                                             getStability(ebmt_X_data)$upper)))
        }
    }
    save(data.plot, file='complete_selection_plot.RData')
} else
    load('complete_selection_plot.RData')

facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
                         rho=function(x){paste("rho=",x,sep="")})

# colors <- c('projpred.iter' = '#2E2E2E',
#             'loc.fdr' = '#819FF7',
#             'ci.90' = '#F78181',
#             'EB.med' = '#298A08')

## Sensitivity vs False discovery rate plot
data.plot.clean <- data.plot %>%
  filter(method != 'test.iter') %>%
  filter(!(method=='projpred.iter' & approach=='data')) %>%
  mutate(method = factor(method, levels = names(colors)))


plot1 <- ggplot(data.plot.clean,aes(x=fdr,y=sensitivity,col=method)) +
    facet_grid(rho~n, labeller=facet.labels) +
    geom_point(aes(shape=approach),size=2) +
    ## geom_label_repel(data=filter(data.plot.clean,!is.na(alpha)),
    ##                 aes(label=alpha),
    ##                 segment.size = 0.2,
    ##                 size = 3,
    ##                 alpha = 0.6) +
    scale_shape_manual(values = shapes) +
    scale_colour_manual(values = colors) +
    guides(
      shape = guide_legend(order = 1),
      colour = guide_legend(order = 2)
    ) +
    geom_line(data=filter(data.plot.clean,method!='projpred.iter'),aes(col=method)) +
    guides(label=FALSE) +
    labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")

ggsave(filename='graphics/sensitivity_vs_fdr_iterated.pdf',plot1,width=10,height=3)

plot2 <- data.plot.clean %>%
  mutate(method = factor(method, levels = rev(names(colors)))) %>%
  ggplot(aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(rho~n, labeller=facet.labels) +
  geom_point(size = 3, alpha = 0.8) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))

ggsave(filename='graphics/stability_iterated.pdf',plot2, width=10,height=3)


