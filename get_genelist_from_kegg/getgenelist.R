rm(list=ls())

library(rvest)
library(dplyr)

getgenelist <- function(url){
  if(!grepl('http://www.kegg.jp',url)){url <- paste0('http://www.kegg.jp',url)}
  kegg <- readLines(url)
  kegg <- kegg[(grep('mapdata',kegg)[2]+1):(grep('</map>',kegg)-1)]
  kegg <- kegg[!grepl('onmouseover',kegg)]
  kegg <- substr(kegg,regexpr('title',kegg)+7,nchar(kegg))
  kegg <- substr(kegg,1,regexpr('\"',kegg)-1)
  kegg <- strsplit(gsub(',','',kegg),' ')
  genelist2 <- lapply(kegg,function(x){
    x <- x[1:(length(x)/2)*2]
    substr(x,2,nchar(x)-1)
  })
  genelist2 <- unique(genelist2)
}

test <- apply(pathmap,1,function(x){
  getgenelist(x[2])
})
