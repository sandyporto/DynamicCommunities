library(igraph)
library('testthat')

setwd("C:/Users/Sandy/Dropbox/2014_Sandy/DynamicCommunities/DynamicCommunities/")

pasta = "C:/Users/Sandy/Dropbox/2014_Sandy/DynamicCommunities/DynamicCommunities/"


nvertices = 300
avgdegree = 20
maxdegree = 40
mixing = 0.05
minsize = 40
maxsize = 80


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

born <- function(g, nmin = minsize, nmax = maxsize, mi = mixing)


##################################################
#Funções Auxiliares
##################################################
  
calculaMixing <- function(g){
  nc = length(unique(V(g)$p))
  temp = 0
  for (i in 1:nc){
    vout = length(E(g)[V(g)[V(g)$p==i] %--% V(g)[V(g)$p!=i]])
    vtotal = length(E(g)[V(g)[V(g)$p==i] %--% V(g)])
    temp = temp +vout/vtotal
  }
  
  return(temp/nc)
}

menorComunidade <- function(g){
  nc = length(unique(V(g)$p))
  temp = 0
  tamanho = vcount(g)
  for (i in 1:nc){
    temp = length(V(g)$p[V(g)$p==i])
    if (temp<tamanho){
      tamanho = temp
    }
  }
  return(tamanho)
}

maiorComunidade <- function(g){
  nc = length(unique(V(g)$p))
  temp = 0
  tamanho = 0
  for (i in 1:nc){
    temp = length(V(g)$p[V(g)$p==i])
    if (temp>tamanho){
      tamanho = temp
    }
  }
  return(tamanho)
}  
  

test_file("testDynamicalCommunities.R")

