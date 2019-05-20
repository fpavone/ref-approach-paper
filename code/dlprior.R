### Dirichlet-Laplace prior for variable selection by ###
###  Bhattacharya, Pati, Pillai and Dunson (2015)     ###
### "Dirichlet–Laplace Priors for Optimal Shrinkage"  ###

# R-code by Maryclare Griffin  https://github.com/maryclare/dlprior
# Code for paper is given at https://github.com/debdeeptamu/Dirichlet_Laplace/blob/master/DL.m

library(GIGrvg)
library(statmod)

sym.sq.root <- function(A) {
  A.eig <- eigen((A + t(A))/2)
  crossprod(t(A.eig$vectors), tcrossprod(diag(sqrt(ifelse(A.eig$values > 0, A.eig$values, 0)),
                                              nrow = nrow(A), ncol = ncol(A)), A.eig$vectors))
}

ei.inv <- function(A) {
  A.eig <- eigen(A)
  crossprod(t(A.eig$vectors), tcrossprod(diag(ifelse(A.eig$values > 0, 1/A.eig$values, 0)), A.eig$vectors))
}

dl.beta <- function(XtX, Xty, sig.sq, psi, lambda) {
  
  if (is.matrix(XtX)) {
    
    p <- nrow(XtX)
    diagonals <- exp(log(lambda) + log(psi)/2)
    phiphit <- tcrossprod(diagonals)
    A <- (XtX/sig.sq)*(phiphit) + diag(p)
    A.inv <- ei.inv(A)
    mean <- crossprod(A.inv*phiphit, Xty/sig.sq)
    # cat("Diagonals\n")
    # print(summary(psi)) # sqrt(psi)
    # cat("A.inv\n")
    # print(summary(c(abs(A.inv))))
    return(crossprod(t(tcrossprod(diagonals, rep(1, p))*sym.sq.root(A.inv)), rnorm(p)) + mean)
  } else if (is.vector(XtX)) {
    p <- length(XtX)
    
    A <- XtX/sig.sq + 1/(exp(2*log(lambda) + log(psi)))
    A.inv <- 1/(A)
    mean <- A.inv*Xty/sig.sq
    var <- A.inv
    
    return(sqrt(var)*rnorm(p) + mean)
  }
}

dl.psi <- function(beta, lambda) {
  p <- length(beta)
  
  psi.tilde <- statmod::rinvgauss(p, mean = lambda/abs(beta), shape = rep(1, p))
  return(1/psi.tilde)
  # return(1/GIGrvg::rgig(p, lambda = -1/2, chi = 1, psi = beta^2/lambda^2))
  # return(GIGrvg::rgig(p, lambda = 1/2, psi = 1, chi = beta^2/lambda^2)) # Gives same results as above
  
  # I compared statmod::rinvgauss, an inverse gaussian generator from Wikipedia used in Debdeep's code and
  # GIGrvg::rgig
  # - the Wiki generator and GIGrvg::rgig generators were all numerically unstable/performed poorly for very large mu, lambda/abs(beta)
}

dl.tau <- function(beta, phi, a) {
  p <- length(beta)
  return(GIGrvg::rgig(1, lambda = p*a - p, psi = 1, chi = 2*sum(abs(beta/(phi)))))
}

dl.lambda <- function(beta, a) {
  p <- length(beta)
  return(GIGrvg::rgig(p, lambda = 1 - a, chi = 1, psi = 2*abs(beta)))
}

# Uses slice sampling algorithm of Damien, Wakefield and Walker (1999)
dl.phi <- function(beta, a, s) {
  # - Draw slice variable
  u <- runif(p, 0, exp(-1/(2*s)))
  # - Convert slice variable to lower bound on s
  lb <- 1/(2*log(1/u))
  # - Find the gamma cdf value that corresponds to this lower bound
  Flb <- pgamma(lb, shape = 1 - a, rate = abs(beta))
  # - Use inverse CDF method to draw gamma random variable
  uu <- pmin(runif(p, Flb, 1), 1-(1e-16))
  s <- qgamma(uu,shape = 1-a,rate = abs(beta))
  # - Convert back to t and phi
  t <- 1/s
  phi <- t/sum(t)
  phi[phi <= (1e-20)] <- (1e-20)
  return(list("phi" = phi, "s" = s))
}

dl.sampler <- function(y, X, a, sig.sq, num.samp = 10000,
                       burn.in = 0, thin = 1,
                       print.iter = FALSE, lambdapar = TRUE) {
  
  p <- ncol(X)
  n <- nrow(X)
  
  XtX <- crossprod(X)
  if (max(abs(XtX[lower.tri(XtX, diag = FALSE)])) == 0) {
    diagX <- TRUE
    XtX <- diag(crossprod(X))
  }
  Xty <- crossprod(X, y)
  
  betas <- psis <- matrix(nrow = num.samp, ncol = p)
  if (!lambdapar) {
    phis <- matrix(nrow = num.samp, ncol = p)
    taus <- numeric(num.samp)
  } else {
    lambdas <- matrix(nrow = num.samp, ncol = p)
  }
  
  # Starting values
  psi <- rep(1, p)
  if (!lambdapar) {
    t <- s <- rep(1, p)
    phi <- t/sum(t)
    tau <- 1
    lambda <- phi*tau
  } else {
    lambda <- rep(1, p)
  }
  
  
  for (i in 1:(num.samp*thin + burn.in)) {
    
    if (print.iter) {cat("i = ", i, "\n")}
    
    beta <- dl.beta(XtX = XtX, Xty = Xty, sig.sq = sig.sq, psi = psi,
                    lambda = lambda)
    psi <- dl.psi(beta = beta, lambda = lambda)
    
    if (!lambdapar) {
      tau <- dl.tau(beta = beta, phi = phi, a = a)
      samp.phi <- dl.phi(beta = beta, a = a, s = s)
      phi <- samp.phi$phi
      s <- samp.phi$s
      lambda <- phi*tau
    } else {
      lambda <- dl.lambda(beta = beta, a = a)
    }
    
    if (i > burn.in & (i - burn.in)%%thin == 0) {
      betas[(i - burn.in)/thin, ] <- beta
      psis[(i - burn.in)/thin, ] <- psi
      if (!lambdapar) {
        phis[(i - burn.in)/thin, ] <- phi
        taus[(i - burn.in)/thin] <- tau
      } else {
        lambdas[(i - burn.in)/thin, ] <- lambda
      }
    }
  }
  
  if (!lambdapar) {
    res <- list("betas" = betas, "psis" = psis, "phis" = phis, "taus" = taus)
  } else {
    res <- list("betas" = betas, "psis" = psis, "lambdas" = lambdas)
  }
  
  return(res)
  
}