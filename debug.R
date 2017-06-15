
rm(list=ls())

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

selbycor <- function(yi,w,nx){
  w.cor <- t(abs(cor(yi,w)))
  w.cor >= quantile(w.cor,1-nx/ncol(w))
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

#2 stage ols

calc_d <- function(wi,x,yi,rho=0,Z=0,U=0){
  if(sum(abs(Z),abs(U))==0){
    Z <- 0
  } else {
    Z <- Z-U
  }
  out <- ginv(t(wi)%*%x%*%ginv(t(x)%*%x)%*%t(x)%*%wi+rho*diag(ncol(wi))) %*%
    (t(wi)%*%x%*%ginv(t(x)%*%x)%*%t(x)%*%yi+rho*Z)
  rownames(out) <- colnames(wi)
  out
}

#sparse coef given Y, X, i within l
equali <- function(Y,X,i,lambda1=0.5,rho=1,wsel=NULL,Z=0,U=0){
  yi <- Y[,i,drop=F]
  x <- X
  w <- cbind(Y[,-i],x)
  if(is.null(wsel)){
    wsel <- selbycor(yi,w,ncol(x)*lambda1)  
  }
  wi <- w[,wsel,drop=F]
  di <- calc_d(wi,x,yi,rho,Z,U)
  wsel[wsel] <- di
  wsel
}

#sparse coef given Y,X, i cross L
maini <- function(raw,i=1,rho=1,lambda1=.7,lambda2=.1,a=.3,itnmax=100){
  #config_out
  # L <- length(raw)
  # M <- ncol(Y)
  # rho <- 1
  # lambda1 <- .7
  # lambda2 <- .1
  # a <- 0.3
  # i <- 1
  # itnmax <- 100
  #config_in
  L <- length(raw)
  M <- ncol(raw[[1]]$Y)
  J <- ncol(raw[[1]]$X)
  out <- matrix(0,M+J,L,dimnames=list(c(colnames(raw[[1]]$Y),colnames(raw[[1]]$X)),1:L))
  #init0
  d0 <- sapply(1:L,function(l){
    Yl <- raw[[l]]$Y
    Xl <- raw[[l]]$X
    equali(Yl,Xl,i,lambda1,NULL,rho=rho)
  })
  z0 <- d0
  wsel <- rowSums(z0!=0)>0
  u0 <- d0-z0
  d1 <- d0
  #initm
  itn <- 0
  while(TRUE){
    itn <- itn+1
    if(itn>itnmax){break}
    rho <- rho * a
    lambda2 <- lambda2 * a
    d1 <- sapply(1:L,function(l){
      Yl <- raw[[l]]$Y
      Xl <- raw[[l]]$X
      equali(Y=Yl,X=Xl,i=i,lambda1=NULL,rho=rho,wsel=wsel,Z=z0[wsel,l],U=u0[wsel,l])
    })
    if(pn(d1-d0)<=1e-8){break}
    z1 <- positive(1-sqrt(L)*lambda2/rho/apply(d1-u0,1,pn)) * (d1-u0)
    z1[is.na(z1)] <- 0
    wsel <- rowSums(z1!=0)>0
    if(sum(wsel)==0){
      d1[,] <- 0
      break
    }
    u1 <- u0 + d1 - z1
    d0 <- d1; z0 <- z1; u0 <- u1
  }
  #result
  out[-i,] <- d1
  out
}

#full network 
main <- function(raw,rho=1,lambda1=.7,lambda2=.1,a=.3,itnmax=100){
  out <- lapply(1:ncol(raw[[1]]$Y),maini,raw=raw)
  abind(out,along=0)
}

############################
# Sample data
############################

Y <- sapply(1:5,function(x) rnorm(1000))
colnames(Y) <- paste0('y',1:5)
X <- sapply(1:10,function(x) rnorm(1000))
colnames(X) <- paste0('x',1:10)
Y[,1] <- 3 * Y[,3] + 5 * Y[,5] + 2 * X[,2] - 2 * X[,3]
Y[,2] <- 1 * Y[,3] + 3 * X[,5] - 5 * X[,9]
Y <- scale(Y)
X <- scale(X)

raw <- lapply(1:3,function(l){
  Y <- dummy(Y,.1)
  X <- dummy(X,.3)
  list(Y=Y,X=X)
})

############################
# Main
############################

test <- main(raw)
