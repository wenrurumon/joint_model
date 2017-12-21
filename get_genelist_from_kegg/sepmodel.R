
rm(list=ls())

####################################

setwd('e:\\uthealth\\signaling\\codes\\')
library(flare)
library(grplasso)
library(data.table)
library(dplyr)
library(igraph)
source('sparse_2sem_final.R')
source('local_cnif_macro.R')
source('CNIF.R')
sourceCpp("score_function_regression.cpp")
sourceCpp("simple_cycle.cpp")
sourceCpp("initial_sem.cpp")

plotnet <- function(x,mode='undirected'){
  diag(x) <- 0
  plot(graph_from_adjacency_matrix(t(x),mode=mode),
       edge.arrow.size=.1,
       vertex.size=3,
       vertex.label.cex=1,
       edge.width=.1)
}
plot2 <- function(x1,x2){
  par(mfrow=c(1,2))
  plotnet(x1[[1]],'directed')
  plotnet(x2[[1]],'directed')
  par(mfrow=c(1,1))
}

####################################

setwd('E:\\uthealth\\joint_model\\sepcompare')
load('datainpathways.rda')

####################################

test <- function(p){
  print(p<<-p)
  rush <- pathinrush[[p]]
  adni <- pathinadni[[p]]
  gene <- colnames(rush)[colnames(rush)%in%colnames(adni)]
  rush <- rush[,match(gene,colnames(rush)),drop=F]
  adni <- adni[,match(gene,colnames(adni)),drop=F]
  #SEM
  sem_rush <- sparse_2sem(rush,lambda=0.12)
  sem_adni <- sparse_2sem(adni,lambda=0.12)
  #CNIF
  sem_rush <- CNIF(rush,init.adj=sem_rush[[1]],max_parent=2)
  sem_adni <- CNIF(adni,init.adj=sem_adni[[1]],max_parent=2)
  #SEM
  m.rush <- sparse_2sem(rush,Y.fixed=sem_rush)
  m.adni <- sparse_2sem(adni,Y.fixed=sem_adni)
  #G
  g.rush <- graph_from_adjacency_matrix(t(m.rush[[1]]),mode='directed')
  g.adni <- graph_from_adjacency_matrix(t(m.adni[[1]]),mode='directed')
  par(mfrow=c(1,2));plot(g.rush);plot(g.adni);par(mfrow=c(1,1))
  #check
  x <- m.adni$eq_scorenet
  x <- lapply(1:nrow(x),function(i){
    x <- x[i,c(2,4)]
    shortest_paths(g.rush,from=x[2],to=x[1])$vpath
  })
  print(xadni <- c(length(x),sum(sapply(lapply(x,unlist),length)>0)/length(x)))
  x <- m.rush$eq_scorenet
  x <- lapply(1:nrow(x),function(i){
    x <- x[i,c(2,4)]
    shortest_paths(g.adni,from=x[2],to=x[1])$vpath
  })
  print(xrush <- c(length(x),sum(sapply(lapply(x,unlist),length)>0)/length(x)))
}
trytest <- function(p){try(test(p))}

rlt <- lapply(1:length(genelist),trytest)

####################################
