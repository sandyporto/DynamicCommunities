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

born <- function(g, nmin = minsize, nmax = maxsize, mi = mixing){
  taminicial = vcount(g)
  tamcomu = sample(nmin:nmax,1)
  idcomu = max(V(g)$p)+1
  
  for (i in 1:tamcomu){
    g = add.vertices(g,1)
    V(g)[vcount(g)]$p = idcomu
    if (i==1){
      g = add.edges(g,c(vcount(g),sample(1:(vcount(g)-1),1)))
    }else{
      grau = sample(2:(i+1),1)
      
      for (j in 1:grau){
        conexao = sample(c("in","out"),1,replace=F,c(1-mi,mi))
        
        if (conexao=="out"){
          espaco = as.vector(V(g)[V(g)$p!=idcomu])
          v1 = vcount(g)
          v2 = sample(espaco,1)
          
          g = add.edges(g,c(v1,v2))
          g = simplify(g)
          
        }
        if (conexao == "in"){
          espaco = as.vector(V(g)[V(g)$p==idcomu])
          v1 = vcount(g)
          v2 = sample(espaco,1)
          
          g = add.edges(g, c(v1,v2))
          g = simplify(g)
          
        }
      }
    }
  }
  
  
}


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

