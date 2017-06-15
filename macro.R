
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

#Calculation in model

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

#equali(Yl,Xl,i,NULL,wsel,Z=z0[wsel,l],U=u0[wsel,l])
equali <- function(Y,X,i,lambda1=0.5,wsel=NULL,Z=0,U=0){
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
