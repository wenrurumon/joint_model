
rm(list=ls())
library(corpcor)
library(flare)

############################
# Macro
############################

lasso <- function(Y,X,lambda){
  slimi <- flare::slim(X=X,Y=Y,lambda=lambda,rho=1,method='lasso',verbose=FALSE)
  Xsel <- (slimi$beta!=0)
  X <- X[,Xsel,drop=F]
  list(Xsel=Xsel,lm=lm(Y~X))
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

############################
# Sample data
############################

Y <- scale(data.matrix(iris[,1:4]))
colnames(Y) <- paste0('y',1:4)
X <- scale(qpca(sapply(1:ncol(Y),function(i){lm(Y[,i]~Y[,-i])$residual}),lambda=.2)$X)
colnames(X) <- paste0('x',1:3)
rawi <- cbind(Y,X)

raw <- lapply(1:3,function(i){
  set.seed(i)
  rawi <- rawi + rnorm(length(rawi),0,.8)
  list(Y=scale(rawi)[,1:4,drop=F],X=scale(rawi)[,-1:-4,drop=F])
})

############################
# Initialization and Interation 1
############################

#Global Config
YX <- raw
lambda <- .1
rho <- 1
L <- length(YX)
J <- ncol(raw[[1]]$Y)
K <- ncol(raw[[1]]$X)

#for i
i <- 1

#init: DELTAi(0), Zi(0), U0
D0 <- lapply(1:L,function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  u0 <- 0
  di0 <- (ginv(t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wi + rho)%*%
            t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%yi)
  di0
})
Z0 <- D0
U0 <- lapply(1:L,function(l){
  D0[[l]]-Z0[[l]]
})

#interation: DELTAi(m)

D1 <- lapply(1:L,function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  d1 <- (ginv(t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wi + rho)%*%
           t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%yi + 
           rho * (Z0[[l]] - U0[[l]]))
  d1
})

#interation: Z(m+1)

Z1 <- lapply(1:L,function(l){
  Zr1 <- positive(1 - lambda * sqrt(L) / (rho * pn(D1[[l]][1:(J-1)] - U0[[l]][1:(J-1)]))) * (D1[[l]] - U0[[l]])
  Zr1 <- Zr1[1:(J-1)]
  Zb1 <- positive(1 - lambda * sqrt(L) / (rho * pn(D1[[l]][-(1:(J-1))] - U0[[l]][-(1:(J-1))]))) * (D1[[l]] - U0[[l]])
  Zb1 <- Zb1[-(1:(J-1))]
  Z1 <- as.matrix(c(Zr1,Zb1))
  Z1
})

#interation: U(m+1)

U1 <- lapply(1:L,function(l){
  U0[[l]] + D1[[l]] - Z1[[l]]
})

#check convergence

D.check <- pn(do.call(cbind,D1) - do.call(cbind,D0))
Z.check <- pn(do.call(cbind,Z1) - do.call(cbind,Z0))
U.check <- pn(do.call(cbind,U1) - do.call(cbind,U0))
if(any(c(D.check,Z.check,U.check) < 0.0001)){
  print('converge')
} else {
  print(c(D.check,Z.check,U.check))
  D0 <- D1; Z0 <- Z1; U0 <- U1
}

############################
# Interation 1+
############################

D1 <- lapply(1:L,function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  d1 <- (ginv(t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wi + rho)%*%
           t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%yi + 
           rho * (Z0[[l]] - U0[[l]]))
  d1
})
Z1 <- lapply(1:L,function(l){
  Zr1 <- positive(1 - lambda * sqrt(L) / (rho * pn(D1[[l]][1:(J-1)] - U0[[l]][1:(J-1)]))) * (D1[[l]] - U0[[l]])
  Zr1 <- Zr1[1:(J-1)]
  Zb1 <- positive(1 - lambda * sqrt(L) / (rho * pn(D1[[l]][-(1:(J-1))] - U0[[l]][-(1:(J-1))]))) * (D1[[l]] - U0[[l]])
  Zb1 <- Zb1[-(1:(J-1))]
  Z1 <- as.matrix(c(Zr1,Zb1))
  Z1
})
U1 <- lapply(1:L,function(l){
  U0[[l]] + D1[[l]] - Z1[[l]]
})
D.check <- pn(do.call(cbind,D1) - do.call(cbind,D0))
Z.check <- pn(do.call(cbind,Z1) - do.call(cbind,Z0))
U.check <- pn(do.call(cbind,U1) - do.call(cbind,U0))
if(any(c(D.check,Z.check,U.check) < 0.0001)){
  print('converge')
} else {
  print(c(D.check,Z.check,U.check))
  D0 <- D1; Z0 <- Z1; U0 <- U1
}