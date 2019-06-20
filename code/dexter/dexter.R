library(tidyverse)

load('dexter.RData')

n <- nrow(x)
p <- ncol(x)

xx <- as.matrix(x)

constant <- numeric()
for(i in 1:p){
  if(sum(xx[,i])==0 | sum(xx[,i])==n) constant <- c(constant,i)
}

xxx <- xx[,-constant]

true <- cor(xxx,y)

sub <- sample(1:n,n/10,replace = FALSE)


