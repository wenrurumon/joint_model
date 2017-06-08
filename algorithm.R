
############################
# Algorithm
############################

#Global Config
YX <- raw
rho <- 0.1
L <- length(YX)

#for i
i <- 1

#init: DELTAi(0), Zi(0), U0
d0 <- lapply(1:length(YX),function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  u0 <- 0
  di0 <- (ginv(t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wi + rho)%*%
                      t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%yi)
  di0
})
Z0 <- d0
u0 <- 0

#interation: DELTAi(m)

d1 <- lapply(1:length(YX),function(l){
  Yl <- YX[[l]]$Y
  Xl <- YX[[l]]$X
  Wi <- cbind(Yl[,-i,drop=F],Xl)
  yi <- Yl[,i,drop=F]
  d1 <- (ginv(t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%Wi + rho)%*%
      t(Wi)%*%Xl%*%ginv(t(Xl)%*%Xl)%*%t(Xl)%*%yi + 
      rho * (Z0[[l]] - u0))
  d1
})

#interation: Z(m+1), U(m+1)
