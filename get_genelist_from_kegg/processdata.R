
rm(list=ls())

####################################

setwd('E:\\uthealth\\joint_model\\sepcompare')
load('adnidata.rda')
load('rushdata.rda')
load('genelist.rda')

qpca <- function(A,rank=0,ifscale=TRUE){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}

####################################

getdata <- function(p,d){
  gsel <- lapply(p,function(g){
    which(d$map%in%g)
  })
  p <- p[sapply(gsel,length)>0]
  out <- sapply(p,function(g){
    x <- d$data[,which(d$map%in%g),drop=F]
    qpca(x)$X[,1,drop=F]
  })
  colnames(out) <- sapply(p,function(x){x[1]})
  return(out)
}

pathinrush <- lapply(genelist,getdata,d=rushdata)
pathinadni <- lapply(genelist,getdata,d=adnidata)
save(pathinrush,pathinadni,genelist,file='datainpathways.rda')

####################################
