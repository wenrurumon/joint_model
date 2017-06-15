############################
# Sample data
############################

Y <- sapply(1:5,function(x) rnorm(1000))
colnames(Y) <- paste0('y',1:5)
X <- sapply(1:10,function(x) rnorm(1000))
colnames(X) <- paste0('x',1:10)
Y[,1] <- Y[,3] + Y[,5] + X[,2] + X[,3]
Y[,2] <- Y[,3] + Y[,4] + X[,5] + X[,9]
Y <- scale(Y)
X <- scale(X)

raw <- lapply(1:3,function(l){
  Y <- dummy(Y,.1)
  X <- dummy(X,.3)
  list(Y=Y,X=X)
})
