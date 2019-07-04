library(tidyverse)
library(rstanarm)
library(projpred)
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

p0 <- 2 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
rhs_prior <- hs(global_scale=tau0)

## Analysis with original covariates
bootnum <- 100
boot_proj <- boot_step <- matrix(0, ncol = length(pred), nrow = bootnum,
                                 dimnames = list(NULL, pred))
for (i in 1:bootnum) {
  print(i)
  set.seed(5437854+i)
  data_id <- sample(1:dim(df)[1], replace = T)
  fitb <- stan_glm(formula, data = df[data_id, ],
                   prior=rhs_prior, QR=TRUE, seed=i, refresh=0)
  bcvvs <- cv_varsel(fitb, method='forward', cv_method='LOO', nloo=n,
                     verbose = FALSE)
  nv <- suggest_size(bcvvs,alpha=0.1)
  boot_proj[i, bcvvs$vind[1:nv]] <- 1

  steplm <- step(lm(formula, data = df[data_id, ]), direction = 'backward')
  boot_step[i, names(steplm$coefficients)[-1]] <- 1
}

boot_inclusion <- tibble(variable = pred,
                         projpred = colMeans(boot_proj)*100,
                         steplm = colMeans(boot_step)*100)

## Save.image('bodyfat_notebook.Rdata')


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
