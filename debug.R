
rm(list=ls())
library(corpcor)
library(flare)

rho <- 2
lambda <- rho/10

############################
# Macro
############################

lasso <- function(Y,X,lambda){
  slimi <- flare::slim(X=X,Y=Y,lambda=lambda,rho=1,method='lasso',verbose=FALSE)
  Xsel <- (slimi$beta!=0)
  X <- X[,Xsel,drop=F]
  list(Xsel=Xsel,lm=lm(Y~X))
}
lassot <- function(Y,X,lambda,i=1){
  X <- cbind(Y[,-i],X)
  Y <- Y[,i,drop=F]
  lasso(Y,X,lambda)
}

qpca <- function(A,lambda=0){
  A.svd <- svd(A)
  d <- A.svd$d-lambda*A.svd$d[1]
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
}

ginv<-function(A){
  A_svd<-fast.svd(A)
  if(length(A_svd$d)==1){
    A_inv<-A_svd$v%*%as.matrix(1/A_svd$d)%*%t(A_svd$u)
  }else{
    A_inv<-A_svd$v%*%diag(1/A_svd$d)%*%t(A_svd$u)
  }
  return(A_inv)
}

pn <- function(x,p=2){
  (sum(abs(x)^p))^(1/p)
}

positive <- function(x){
  x * (x>0)
}

apply_random <- function(x,s=1){
  out <- x + rnorm(x,0,s)
  scale(out)[,]
}

############################
# Sample data
############################

#Iris data YX

Y <- scale(data.matrix(iris[,1:4]))
colnames(Y) <- paste0('y',1:4)
Y <- cbind(Y,y5 = runif(nrow(Y)),y6=runif(nrow(Y)))

X <- scale(qpca(sapply(1:4,function(i){lm(Y[,i]~Y[,-i])$residual}),lambda=.2)$X)
colnames(X) <- paste0('x',1:ncol(X))
X <- cbind(X,x4=runif(nrow(X)), x5=runif(nrow(X)), x6=runif(nrow(X)))

#Ramdom data YX

# Y <- cbind(y1=runif(100,10,20))
# Y <- scale(cbind(Y,y2=Y[,1]*20+24+runif(100,-2,2)))
# X <- scale(matrix(rnorm(300),100,3,dimnames=list(NULL,c('x1','x2','x3'))))

raw <- lapply(1:3,function(i){
  set.seed(i)
  Y <- apply_random(Y,0.1)
  X <- apply_random(X,1)
  list(Y=Y,X=X)
})

############################
# Initialization and Interation 1
############################

#Global Config
YX <- raw
# rho <- 1
# lambda <- rho/20
L <- length(YX)
J <- ncol(raw[[1]]$Y)
K <- ncol(raw[[1]]$X)

#for i
i <- 1

#init: DELTAi(0), Zi(0), U0
D0 <- sapply(1:L,function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(int=1,Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  u0 <- 0
  # di0 <- (ginv(t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wi + rho * diag(ncol(Wi)))%*%
            # t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%yi)
  di0 <- (ginv(t(Wi)%*%Wi%*%ginv(t(Wi)%*%Wi)%*%t(Wi)%*%Wi + rho * diag(ncol(Wi)))%*%
            t(Wi)%*%Wi%*%ginv(t(Wi)%*%Wi)%*%t(Wi)%*%yi)
  di0
})
Z0 <- D0
U0 <- D0-Z0

#interation: DELTAi(m)

D1 <- sapply(1:L,function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(int=1,Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  d1 <- (ginv(t(Wi)%*%Wi%*%ginv(t(Wi)%*%Wi)%*%t(Wi)%*%Wi + rho * diag(ncol(Wi)))%*%
            t(Wi)%*%Wi%*%ginv(t(Wi)%*%Wi)%*%t(Wi)%*%yi + 
           rho * (Z0[,l] - U0[,l]))
  d1
})

#interation: Z(m)

# Z1 <- positive(1 - lambda/rho * sqrt(L) / sqrt(rowMeans(D0 - U0)^2)) * (D0-U0)
Z1 <- positive(1-lambda/rho * sqrt(L) / abs(D0-U0)) * (D0-U0)
U1 <- U0 + D1 - Z1
# print(paste(pn(D1-D0),pn(Z1-Z0),pn(U1-U0)))

#interation: m+1
itn <- 0
while(itn<1000){
  # print(itn)
  itn <- itn+1
  D0 <- D1; Z0 <- Z1; U0 <- U1
  D1 <- sapply(1:L,function(l){
    Yl <- YX[[l]]$Y
    Xl <- YX[[l]]$X
    Wi <- cbind(int=1,Yl[,-i,drop=F],Xl)
    yi <- Yl[,i,drop=F]
    d1 <- (ginv(t(Wi)%*%Wi%*%ginv(t(Wi)%*%Wi)%*%t(Wi)%*%Wi + rho * diag(ncol(Wi)))%*%
             t(Wi)%*%Wi%*%ginv(t(Wi)%*%Wi)%*%t(Wi)%*%yi + 
             rho * (Z0[,l] - U0[,l]))
    d1
  })
  # Z1 <- positive(1 - lambda/rho * sqrt(L) / sqrt(rowMeans(D0 - U0)^2)) * (D0-U0)
  Z1 <- positive(1-lambda/rho * sqrt(L) / abs(D0-U0)) * (D0-U0)
  Z1[is.na(Z1)] <- 0
  U1 <- U0 + D1 - Z1
  # print(paste(pn(D1-D0),pn(Z1-Z0),pn(U1-U0)))
}

#output
((Z1!=0)+0)[-1,,drop=F]
sapply(raw,function(x){
  y <- x$Y[,1,drop=F]
  x <- cbind(x$Y[,-1],x$X)
  lasso(y,x,lambda=0.08)[[1]]
}) + 0
