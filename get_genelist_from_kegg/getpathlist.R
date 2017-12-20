library(rvest)
url <- 'http://www.genome.jp/kegg/pathway.html'
html <- read_html(url)
title <- html %>% html_nodes('dd a') %>% html_text()
href <- html %>% html_nodes('dd a') %>% html_attr('href')
pathmap <- cbind(title=title,href=href)
