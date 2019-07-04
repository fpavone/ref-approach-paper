library(tidyverse)
library(rstanarm)
library(projpred)
set.seed(453872238)

## 10-fold cross-validation RMSE for projpred vs step.lm with and without noisy covariates

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

K <- 10

## Original covariates
sel <- step(lm(formula, data = df,  x=T,y=T),
                direction = "backward",
                trace = 0)

p0 <- 2 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
rhs_prior <- hs(global_scale=tau0)
fitrhs <- stan_glm(formula, data = df, prior=rhs_prior, QR=TRUE,
                   seed=1513306866, refresh=0)
fitrhs_cvvs <- cv_varsel(fitrhs, method = 'forward', cv_method = 'LOO',
                          nloo=n, verbose = FALSE)
nv <- suggest_size(fitrhs_cvvs, alpha=0.1)
print(fitrhs_cvvs$vind[1:nv])


perm <- sample.int(n)
idx <- ceiling(seq(from = 1, to = n, length.out = K + 1))
bin <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
lmmuss <- list()
lmvsnlmuss <- list()
lmvsnlnvss <- list()
muss <- list()
vsmuss <- list()
vsnvss <- list()
for (k in 1:K) {
    message("Fitting model ", k, " out of ", K)
    omitted <- which(bin == k)
    ## Stepwise regression
    lmfit_k <- lm(formula, data = df[-omitted,, drop=FALSE],  x=T,y=T)
    lmmuss[[k]] <- predict(lmfit_k, newdata = df[omitted, , drop = FALSE])
    selnl_k <- step(lm(formula, data = df[-omitted,, drop=FALSE],  x=T,y=T),
                    direction = "backward",
                    trace = 0)
    lmvsnlmuss[[k]] <- predict(selnl_k, newdata = df[omitted, , drop = FALSE])
    lmvsnlnvss[[k]] <- length(coef(selnl_k))-1
    ## Reference model + projpred
    fit_k <- update(
        object = fitrhs,
        data = df[-omitted,, drop=FALSE],
        weights = NULL,
        refresh = 0
    )
    fit_cvvs_k <- cv_varsel(fit_k, method='forward', cv_method='LOO',
                            nloo=nrow(df[-omitted,, drop=FALSE]))
    nvk <- suggest_size(fit_cvvs_k,alpha=0.1)
    vsnvss[[k]] <- nvk
    proj_k <- project(fit_cvvs_k, nv = nvk, ns = 4000)
    muss[[k]] <- colMeans(posterior_linpred(fit_k,
                                            newdata = df[omitted, , drop = FALSE]))
    vsmuss[[k]] <- colMeans(proj_linpred(proj_k, xnew = df[omitted, , drop = FALSE]))
}
## Stepwise regression
lmmus<-unlist(lmmuss)[order(as.integer(names(unlist(lmmuss))))]
lmvsnlmus<-unlist(lmvsnlmuss)[order(as.integer(names(unlist(lmvsnlmuss))))]
lmvsnlnvs <- unlist(lmvsnlnvss)
rmse_lmfull <- sqrt(mean((df$siri-lmmus)^2))
rmse_lmselnl <- sqrt(mean((df$siri-lmvsnlmus)^2))
## Reference model + projpred
mus<-unlist(muss)[order(as.integer(names(unlist(muss))))]
vsmus<-unlist(vsmuss)[order(as.integer(names(unlist(vsmuss))))]
vsnvs <- unlist(vsnvss)
rmse_full <- sqrt(mean((df$siri-mus)^2))
rmse_proj <- sqrt(mean((df$siri-vsmus)^2))


###### Noisy covariates #######
p <- 100
noise <- array(rnorm((p-13)*n), c(n,p-13))
dfr<-cbind(df,noise=noise)
formula2<-paste(formula,"+",paste(colnames(dfr[,15:(p+1)]), collapse = "+"))

sel2 <- step(lm(formula2, data = dfr,  x=T,y=T),
                direction = "backward",
                trace = 0)

p0 <- 5 # prior guess for the number of relevant variables
tau0 <- p0/(p-p0) * 1/sqrt(n)
hs_prior <- hs(global_scale=tau0)
fitrhs2 <- stan_glm(formula2, data = dfr, prior = hs_prior, QR = TRUE,
                    seed=34523, refresh=0)
fitrhs2_cvvs <- cv_varsel(fitrhs2, method = 'forward', cv_method = 'LOO',
                          nloo=n, verbose = FALSE)
nv2 <- suggest_size(fitrhs2_cvvs, alpha=0.1)
print(fitrhs2_cvvs$vind[1:nv2])

perm <- sample.int(n)
idx <- ceiling(seq(from = 1, to = n, length.out = K + 1))
bin <- .bincode(perm, breaks = idx, right = FALSE, include.lowest = TRUE)
lmmuss2 <- list()
lmvsnlmuss2 <- list()
lmvsnlnvss2 <- list()
muss2 <- list()
vsmuss2 <- list()
vsnvss2 <- list()
for (k in 1:K) {
    message("Fitting model ", k, " out of ", K)
    omitted <- which(bin == k)
    ## Stepwise regression
    lmfit_k2 <- lm(formula2, data = dfr[-omitted,, drop=FALSE],  x=T,y=T)
    lmmuss2[[k]] <- predict(lmfit_k2, newdata = dfr[omitted, , drop = FALSE])
    selnl_k2 <- step(lm(formula2, data = dfr[-omitted,, drop=FALSE],  x=T,y=T),
                    direction = "backward",
                    trace = 0)
    lmvsnlmuss2[[k]] <- predict(selnl_k2, newdata = dfr[omitted, , drop = FALSE])
    lmvsnlnvss2[[k]] <- length(coef(selnl_k2))-1
    ## Reference model + projpred
    fit_k2 <- update(
        object = fitrhs2,
        data = dfr[-omitted,, drop=FALSE],
        weights = NULL,
        refresh = 0
    )
    fit_cvvs_k2 <- cv_varsel(fit_k2, method='forward', cv_method='LOO',
                            nloo=nrow(dfr[-omitted,, drop=FALSE]))
    nvk2 <- suggest_size(fit_cvvs_k2,alpha=0.1)
    vsnvss2[[k]] <- nvk2
    proj_k2 <- project(fit_cvvs_k2, nv = nvk2, ns = 4000)
    muss2[[k]] <- colMeans(posterior_linpred(fit_k2,
                                            newdata = dfr[omitted, , drop = FALSE]))
    vsmuss2[[k]] <- colMeans(proj_linpred(proj_k2, xnew = dfr[omitted, , drop = FALSE]))
}
## Stepwise regression
lmmus2 <- unlist(lmmuss2)[order(as.integer(names(unlist(lmmuss2))))]
lmvsnlmus2 <- unlist(lmvsnlmuss2)[order(as.integer(names(unlist(lmvsnlmuss2))))]
lmvsnlnvs2 <- unlist(lmvsnlnvss2)
rmse_lmfull2 <- sqrt(mean((dfr$siri-lmmus2)^2))
rmse_lmselnl2 <- sqrt(mean((dfr$siri-lmvsnlmus2)^2))
## Reference model + projpred
mus2 <- unlist(muss2)[order(as.integer(names(unlist(muss2))))]
vsmus2 <- unlist(vsmuss2)[order(as.integer(names(unlist(vsmuss2))))]
vsnvs2 <- unlist(vsnvss2)
rmse_full2 <- sqrt(mean((dfr$siri-mus2)^2))
rmse_proj2 <- sqrt(mean((dfr$siri-vsmus2)^2))

save.image( "bodyfat_notebook_RMSE.RData")



