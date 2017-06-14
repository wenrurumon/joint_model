
rm(list=ls())
library(corpcor)
library(flare)

############################
# Macro
############################

lasso <- function(Y,X=NULL,lambda=0.3){
  if(is.list(Y)){
    Y <- do.call(cbind,Y)
  }
  if(is.null(X)){
    X <- Y[,-1]
    Y <- Y[,1,drop=F]
  }
  slimi <- flare::slim(X=X,Y=Y,lambda=lambda,rho=1,method='lasso',verbose=FALSE)
  Xsel <- (slimi$beta!=0)
  X <- X[,Xsel,drop=F]
  list(Xsel=Xsel,lm=lm(Y~X-1))
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
  (sum((x)^p))^(1/p)
}

positive <- function(x){
  x * (x>0)
}

dummy <- function(x,error=1){
  e <- rnorm(length(x),0,error*sd(x))
  x + e
}

#calculation

calc_d <- function(Wli,Xl,Yli,rho=0,Z=0,U=0){
  ginv(t(Wli)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wli+rho*diag(ncol(Wli))) %*%
    (t(Wli)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Yli+rho*(Z-U))
}

############################
# Sample data
############################

Y <- scale(matrix(rnorm(500),100,5,dimnames=list(NULL,paste0('y',1:5))))
X <- scale(matrix(rnorm(300),100,3,dimnames=list(NULL,paste0('x',1:3))))
Y[,1] <- Y[,3] + Y[,5] + X[,1]
Y[,2] <- Y[,3] + Y[,4] + X[,2]
Y <- scale(Y)
X <- scale(X)

raw <- lapply(1:3,function(l){
  list(Y=dummy(Y,0),X=dummy(X,0))
})

#############################
# Config
#############################

l <- 1
i <- 1
rho <- 1
l1 <- 0.1
l2 <- 0.4

#init

L <- length(raw)
M <- ncol(raw[[1]]$Y)
K <- ncol(raw[[1]]$X)

d0 <- sapply(1:L,function(l){
  Yl <- raw[[l]]$Y
  Xl <- raw[[l]]$X
  Yli <-Yl[,i,drop=F]
  Wli <- cbind(Yl[,-i],Xl)
  calc_d(Wli,Xl,Yli,rho)
})
z0 <- d0
u0 <- z0-d0

#itn
itn <- 0
itnmax <- 20
while(TRUE){
  itn <- itn+1
  if(itn>itnmax){break}
  print(itn)
  d1 <- sapply(1:L,function(l){
    Yl <- raw[[l]]$Y
    Xl <- raw[[l]]$X
    Yli <-Yl[,i,drop=F]
    Wli <- cbind(Yl[,-i],Xl)
    d1 <- calc_d(Wli,Wli,Yli,rho,Z=z0[,l],U=u0[,l])
    d1
  })
  # print(d1)
  # print(u0)
  z1 <- rep(0,nrow(d1))
  z1[1:(M-1)] <- positive(1 - l1 * sqrt(L)/apply(d1-u0,1,pn))[1:M-1]
  z1[-1:(1-M)] <- positive(1 - l2 * sqrt(L)/apply(d1-u0,1,pn))[-1:(1-M)]
  z1 <- ifelse(is.na(z1),0,z1) * (d1-u0)
  u1 <- u0 + d1 - z1
  print((z1!=0)+0)
  print(pn(d1-d0))
  d0 <- d1; z0 <- z1; u0 <- u1
}




