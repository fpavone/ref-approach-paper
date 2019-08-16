library(tidyverse)
library(ggrepel)
theme_set(theme_light() +
          theme(strip.background = element_rect(color='black',fill='white'),
                strip.text.x = element_text(color='black'),
                strip.text.y = element_text(color='black')))
library(latex2exp)
source("getStability.R")

## load('iterated_n50_rho0.5.Rdata')

alpha_vec <- c(0.05,0.1,0.16,0.2,0.25)

## Sensitivity vs FDR, stability for test and alpha = 0.16

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
                 nn,'_rho',rr,'_3alphas.Rdata',sep = ''))
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
                               method = c('iter.test',
                                          'iter.test',
                                          'iter.cv',
                                          'iter.cv'),
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

facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
                         rho=function(x){paste("rho=",x,sep="")})

## Sensitivity vs False discovery rate plot
data.plot.clean <- data.plot %>%
    filter(method!='iter.test') %>%
    filter(!(method=='iter.cv' & approach=='data'))

colors <- c('iter.cv' = '#2E2E2E',
            'loc.fdr' = '#819FF7',
            'EB.med' = '#298A08',
            'ci.90' = '#F78181')

plot1 <- ggplot(data.plot.clean,aes(x=fdr,y=sensitivity,col=method)) +
    facet_grid(rho~n, labeller=facet.labels) +
    geom_point(aes(shape=approach),size=2) +
    ## geom_label_repel(data=filter(data.plot.clean,!is.na(alpha)),
    ##                 aes(label=alpha),
    ##                 segment.size = 0.2,
    ##                 size = 3,
    ##                 alpha = 0.6) +
    scale_shape_manual(values = c(16,15)) +
    scale_colour_manual(values = colors) +
    geom_line(data=filter(data.plot.clean,method!='iter.cv'),aes(col=method)) +
    guides(label=FALSE) +
    labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")

ggsave(filename='../paper/graphics/sensitivity_vs_fdr_iterated.pdf',plot1,width=10,height=3)

plot2 <- ggplot(data.plot.clean,aes(y=stab.mean,x=method,col=approach)) +
  facet_grid(rho~n, labeller=facet.labels) +
  geom_point(size=2.5) +
  geom_linerange(aes(ymin=stab.low,ymax=stab.up)) +
  coord_flip() +
  labs(x="",y="Stability", col="Approach") +
  scale_color_manual(values=c("#819FF7","#FAAC58"))

ggsave(filename='../paper/graphics/stability_iterated.pdf',plot2, width=10,height=3)



## #### Test vs K10 ####


## data.plot <- tibble(n = numeric(),
##                     rho = numeric(),
##                     alpha = numeric(),
##                     method = character(),
##                     approach = character(),
##                     sensitivity = numeric(),
##                     ## sensitivity.sd = numeric(),
##                     fdr = numeric(),
##                     ## fdr.sd = numeric(),
##                     stab.low = numeric(),
##                     stab.mean = numeric(),
##                     stab.up = numeric())


## for(nn in c(50,70,100)){
##   for(rr in c(0.3,0.5)){
##       load(paste('iterated_cv_n',nn,'_rho',rr,'.Rdata',sep = ''))

##       for(j in 1:5){
##           X_ref.cv <- X_proj_list[[j]]
##           X_ref.test <- X_proj_test_list[[j]]
##           X_data.cv <- X_lasso_list[[j]]
##           X_data.test <- X_lasso_test_list[[j]]

##           stab_ref.cv <- getStability(X_ref.cv)
##           stab_ref.test <- getStability(X_ref.test)
##           stab_data.cv <- getStability(X_data.cv)
##           stab_data.test <- getStability(X_data.test)

##           data.plot <- data.plot %>%
##               bind_rows(tibble(n = rep(n,4),
##                                rho = rep(rho,4),
##                                alpha = rep(alpha_vec[j],4),
##                                method = c(rep('test',2),
##                                           rep('cv',2)),
##                                approach = rep(c('ref','data'),2),
##                                sensitivity = c(avg.sensitivity(X_ref.test),
##                                                avg.sensitivity(X_data.test),
##                                                avg.sensitivity(X_ref.cv),
##                                                avg.sensitivity(X_data.cv)),
##                                fdr = c(avg.fdr(X_ref.test),
##                                        avg.fdr(X_data.test),
##                                        avg.fdr(X_ref.cv),
##                                        avg.fdr(X_data.cv)),
##                                stab.low = c(stab_ref.test$lower,
##                                             stab_data.test$lower,
##                                             stab_ref.cv$lower,
##                                             stab_data.cv$lower),
##                                stab.mean = c(stab_ref.test$stability,
##                                              stab_data.test$stability,
##                                              stab_ref.cv$stability,
##                                              stab_data.cv$stability),
##                                stab.up = c(stab_ref.test$upper,
##                                            stab_data.test$upper,
##                                            stab_ref.cv$upper,
##                                            stab_data.cv$upper)))
##       }
##   }
## }

## facet.labels <- labeller(n = function(x){paste("n=",x,sep="")},
##                          rho=function(x){paste("rho=",x,sep="")})

## ## Sensitivity vs False discovery rate plot
## data.plot %>%
##     filter(alpha==0.16) %>%
##     ggplot(aes(x=fdr,y=sensitivity,col=method)) +
##     facet_grid(rho~n, labeller=facet.labels) +
##     geom_point(aes(shape=approach),size=2.5) +
##     ## geom_errorbar(aes(ymin=sensitivity-sensitivity.sd, ymax=sensitivity+sensitivity.sd)) +
##     ## geom_errorbarh(aes(xmin=fdr-fdr.sd, xmax=fdr+fdr.sd)) +
##     scale_shape_manual(values = c(16,15)) +
##     geom_line(aes(col=method)) +
##     labs(x="False discovery rate",y="Sensitivity", shape="Approach", col="Method")


