# Function to fit the reference model properly, choosing the
# SPCA threhsold by cross-validation maximazing the lpd

cv_threshold <- function(x,y,nc = 5, K = 5, sup.only.flag = FALSE){
  n <- nrow(x)
  p <- ncol(x)
  
  folds <- sample(1:K,n,replace=TRUE)
  
  thresholds <- seq(0.2,0.9,by=0.05)
  logdens_fit <- numeric(length(thresholds))
  i <- 0
  for(thresh in thresholds){
    i <- i + 1
    for(k in 1:K){
      sel <- which(folds==k)
      x_k <- x[sel,]
      y_k <- y[sel]
      
      dr_k <- spca(as.matrix(x_k),y_k, nctot=nc, screenthresh = thresh, 
                   window = p, sup.only = sup.only.flag, verbose = F)
      z_k <- predict(dr_k, x_k)
      
      model_data_k <- list(y=y_k,X=z_k,p=ncol(z_k),n=length(y_k),s_max=dr_k$sdev[1])
      fit_k <- sampling(fit_model, data = model_data_k, 
                      chains = 2, iter = 2000, seed = 45342)  
      
      draws <- as.matrix(fit_k) # posterior draws
      sigma <- draws[,'sigma'] # noise std
      beta <- draws[,2:(ncol(z_k)+1)] # regression coefficients
      alpha <- draws[,'intercept'] # intercept
      
      predfun <- function(zt){
        colMns <- colMeans(x_k)
        colSds <- apply(x_k,2,sd)
        zt <- t(apply(zt,1,function(x){(x-colMns)/colSds}))%*%dr_k$w
        return(t( beta %*% t(zt) + alpha ))   # n x s
      }
     
      mufit <- predfun(x[-sel,]) 
      logdens_fit[i] <- logdens_fit[i] + sum(dnorm(y[-sel],mean=mufit,sd=sigma))/length(y[-sel])
    }
  }
  
  return(list(optimal=thresholds[which.max(logdens_fit)],all=logdens_fit))
}


res <- cv_threshold(x,y)
