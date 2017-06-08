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
  rawi <- rawi + rnorm(length(rawi),0,0.3)
  list(Y=scale(rawi)[,1:4,drop=F],X=scale(rawi)[,-1:-4,drop=F])
})
