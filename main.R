############################
# Algorithm
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
  D0 <- D1; Z0 <- Z1; U0 <- U1
}
