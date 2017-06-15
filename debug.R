
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
  if(sum(abs(Z),abs(U))==0){
      Z <- 0
    } else {
      Z <- Z-U
    }
  ginv(t(Wli)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wli+rho*diag(ncol(Wli))) %*%
    (t(Wli)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Yli+rho*Z)
}

############################
# Sample data
############################

Y <- sapply(1:5,function(x) rnorm(100))
colnames(Y) <- paste0('y',1:5)
X <- sapply(1:3,function(x) rnorm(100))
colnames(X) <- paste0('x',1:3)
Y[,1] <- Y[,3] + Y[,5] + X[,1]
Y[,2] <- Y[,3] + Y[,4] + X[,2]
Y <- scale(Y)
X <- scale(X)

raw <- lapply(1:3,function(l){
  list(Y=dummy(Y,0.1),X=dummy(X,0.1))
})

#############################
# Config
#############################

l <- 1
i <- 1
rho <- .3

#init

L <- length(raw)
M <- ncol(raw[[1]]$Y)
K <- ncol(raw[[1]]$X)

d0 <- sapply(1:L,function(l){
  Yl <- raw[[l]]$Y
  Xl <- raw[[l]]$X
  Yli <-Yl[,i,drop=F]
  # Yl <- Xl %*% ginv(t(Xl)%*%Xl) %*% t(Xl) %*% Yl
  Wli <- cbind(Yl[,-i],Xl)
  calc_d(Wli,Xl,Yli,rho)
})
# z0 <- positive(1-rho*sqrt(L)/apply(d0,1,pn)) * d0
z0 <- d0
u0 <- d0 - z0
print(d0)
print(z0)

#itn
d1 <- sapply(1:L,function(l){
  Yl <- raw[[l]]$Y
  Xl <- raw[[l]]$X
  Yli <-Yl[,i,drop=F]
  Wli <- cbind(Yl[,-i],Xl)
  d1 <- calc_d(Wli,Xl,Yli,rho,Z=z0[,l],U=u0[,l])
  d1
})
z1 <- positive(1-sqrt(L)/rho/apply(d0,1,pn))
z1 <- ifelse(is.na(z1),0,z1) * (d1-u0)
u1 <- u0 + d1 - z1
print((z1!=0)+0)
print(pn(d1-d0))

#itn+1
z0 <- z1; d0 <- d1; u0 <- u1
d1 <- sapply(1:L,function(l){
  Yl <- raw[[l]]$Y
  Xl <- raw[[l]]$X
  Yli <-Yl[,i,drop=F]
  Wli <- cbind(Yl[,-i],Xl)
  d1 <- calc_d(Wli,Xl,Yli,rho,Z=z0[,l],U=u0[,l])
  d1
})
z1 <- positive(1-sqrt(L)/rho/apply(d1-u0,1,pn))
z1 <- ifelse(is.na(z1),0,z1) * (d1-u0)
u1 <- u0 + d1 - z1
print((z1!=0)+0)
print(pn(d1-d0))



