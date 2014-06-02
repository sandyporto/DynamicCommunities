library(igraph)

setwd("C:/Users/Sandy/Dropbox/2014_Sandy/DynamicCommunities/DynamicCommunities/")

pasta = "C:/Users/Sandy/Dropbox/2014_Sandy/DynamicCommunities/DynamicCommunities/"



criarGrafoInicial <- function(){
  arquivo = paste(pasta,"main.exe",sep="")
  
  system(arquivo)
  
  arquivo = paste(pasta,"network.dat",sep="")
  rede = as.matrix(read.table(arquivo))
  arquivo = paste(pasta,"community.dat",sep="")
  comus = as.matrix(read.table(arquivo))
  
  G = graph.edgelist(rede,directed=F)
  G = simplify(G)
  V(G)$p = 0
  V(G)$p = comus[,2]
  
  return(G)
}


grafoInicial = criarGrafoInicial()