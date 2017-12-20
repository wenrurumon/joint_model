
rm(list=ls())

library(rvest)
library(dplyr)

url <- 'http://www.genome.jp/kegg/pathway.html'
html <- read_html(url)
title <- html %>% html_nodes('dd a') %>% html_text()
href <- html %>% html_nodes('dd a') %>% html_attr('href')
pathmap <- data.frame(title=title,href=href)

#pathmap <- filter(pathmap,grepl('hsa',href))
