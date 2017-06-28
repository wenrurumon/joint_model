
rm(list=ls())

#Nodes for AD KEGG
html <- readLines('http://www.genome.jp/kegg-bin/show_pathway?hsa05010')
nodes <- grep('area shape=rect\tcoords',html,value=T)
nodes <- lapply(strsplit(nodes,'\"'),function(x)x[[4]])
nodes2 <- lapply(nodes,function(x){
  x <- strsplit(x,', ')[[1]]
  do.call(rbind,strsplit(x,'\\(|: |)'))[,2]
})

#data A
setwd('C:\\Users\\zhu2\\Documents\\joint_model')
load('Arush.rda')
load('Yrs.rda')
load('geneinpath.rda')

#Process
A <- A.rush[match(rownames(Yrs),rownames(A.rush)),colnames(A.rush)%in%unlist(nodes2)]
nodes2 <- nodes2[sapply(nodes2,function(x){any(x%in%colnames(A))})]
A <- sapply(nodes2,function(x){
  rowMeans(scale(A[,colnames(A)%in%x,drop=F]))
})
colnames(A) <- sapply(nodes2,function(x) x[[1]])
A <- scale([,match(unique(colnames(A)),colnames(A))])

As <- lapply(1:ncol(Yrs),function(i){
  A * Yrs[,i]
})

raw <- lapply(As,function(x){
  Y <- x;
  # X <- matrix(rnorm(length(x)),ncol=ncol(x),nrow=nrow(x))
  # colnames(X) <- paste0("X",1:ncol(x))
  X <- matrix(0,ncol=0,nrow=nrow(x))
  list(Y=Y,X=X)
})

plotnet <- function(x,mode='undirected'){
  plot(igraph::graph_from_adjacency_matrix(x,mode),
     edge.arrow.size=.2,
     vertex.size=3,
     vertex.label.cex=1,
     edge.width=.1)
}

#############################################################
#Expression Network
#############################################################

setwd('C:\\Users\\zhu2\\Documents\\signaling\\codes\\')
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")
setwd('C:\\Users\\zhu2\\Documents\\getpathway')
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)
setwd("C:/Users/zhu2/Documents/joint_model")

sem <- sparse_2sem(Y=A,lambda=0.4)
plotnet(sem[[1]]>0.8)
sem_cnif <- CNIF(data=A,init.adj=sem[[1]]>0.8,max_parent=3)
sum(sem_cnif)
plotnet(sem_cnif,'directed')

#############################################################
#Joint Model
#############################################################

# devtools::install_github("wenrurumon/mysrc/jointnet",force=T)
library(jointnet)
library(corpcor)
library(abind)
sem_deconv <- (jointnet(raw,rho=1,lambda1=0.3,lambda2=0.028,a=0.3,itnmax=100)[,,1]>0)+0
rownames(sem_deconv) <- colnames(sem_deconv)
sum((sem_deconv+t(sem_deconv))>0)/2
Apool <- abind(As,along=0)
Apool <- apply(Apool,3,as.vector)
sem_dcnif <- CNIF(data=Apool,init.adj=sem_cnif,max_parent=3)
plotnet(sem_dcnif,'directed')

#############################################################
#Plot
#############################################################

save(sem_cnif,sem_dcnif,file='result.rda')
par(mfrow=c(1,2))
plotnet(sem_cnif,'directed');plotnet(sem_dcnif,'directed')
par(mfrow=c(1,1))

