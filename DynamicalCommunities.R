
##################################################
#Funções Principais
##################################################


criarGrafoInicial <- function(p){
  arquivo = paste(p,"main.exe",sep="")
  
  system(arquivo)
  
  arquivo = paste(p,"network.dat",sep="")
  rede = as.matrix(read.table(arquivo))
  arquivo = paste(p,"community.dat",sep="")
  comus = as.matrix(read.table(arquivo))
  
  G = graph.edgelist(rede,directed=F)
  G = simplify(G)
  V(G)$p = 0
  V(G)$p = comus[,2]
  
  if (maiorComunidade(G)>maxsize){
    G = criarGrafoInicial()
  }
  
  return(G)
}

born <- function(g, nmin = minsize, nmax = maxsize, dmax = maxdegree, mi = mixing){
  taminicial = vcount(g)
  tamcomu = sample(nmin:nmax,1)
  if (taminicial != 0){
    idcomu = max(V(g)$p)+1
  }else{
    idcomu = 1
  }
 
  for (i in 1:tamcomu){
    g = add.vertices(g,1)
    V(g)[vcount(g)]$p = idcomu
    if (i==1){
      if (taminicial != 0){
        g = add.edges(g,c(vcount(g),sample(1:(vcount(g)-1),1)))
      }
    }else{
      auxgrau = min((i+1),maxdegree)
      grau = sample(2:auxgrau,1)
      
      for (j in 1:grau){
        conexao = sample(c("in","out"),1,replace=F,c(1-mi,mi))
        
        if(taminicial == 0){
          conexao = "in"
        }
        
        if (conexao=="out"){
          espaco = as.vector(V(g)[V(g)$p!=idcomu])
          v1 = vcount(g)
          v2 = sample(espaco,1)
          
          while(degree(g,v2) == maxdegree){
            v2 = sample(espaco,1)
          }
          
          g = add.edges(g,c(v1,v2))
          g = simplify(g)
          
        }
        if (conexao == "in"){
          espaco = as.vector(V(g)[V(g)$p==idcomu])
          v1 = vcount(g)
          v2 = sample(espaco,1)
          
          while(degree(g,v2) == maxdegree){
            v2 = sample(espaco,1)
          }
          
          g = add.edges(g, c(v1,v2))
          g = simplify(g)
          
        }
      }
    }
  }

  dens = densidadeComunidade(g)
  narestasout = length(E(g)[V(g)[V(g)$p==idcomu] %--% V(g)[V(g)$p!=idcomu]])
  comumaxavgdegree = (((maxdegree*tamcomu)/2) - narestasout)/(tamcomu*(tamcomu-1)/2)
  
  d = min(mean(dens),comumaxavgdegree)
  d = max(d,(avgdegree-1)/tamcomu)
  aux = T

  while(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])) < (d) && aux){
    espaco = as.vector(V(g)[V(g)$p==idcomu])
    espaco = espaco[degree(g,espaco) < maxdegree]
    
    if (length(espaco)>=2){
      v1 = sample(espaco,1)
      v2 = sample(espaco,1)
      
      g = add.edges(g,c(v1,v2))
      g = simplify(g)
    }else{
      aux = F
    }
    
  }

  if (calculaMixing(g) < (mi-toleranciamixing)){
    g = corrigeMixing(g,mi)
  }

  return(g)
}

extinction <- function(g, comu = 0){
  if(vcount(g)==0){
    return(g)
  }
  
  if (comu==0){
    idcomu = sample(V(g)$p,1)
  }else{
    idcomu=comu
  }
  
  aux = c(1:vcount(g))
  aux = aux[V(g)$p==idcomu]
  
  while(length(aux)>1){
    g = delete.vertices(g,sample(aux,1))
    aux = c(1:vcount(g))
    aux = aux[V(g)$p==idcomu]
  }
  
  g = delete.vertices(g,aux)
  
  return(g)
}


##################################################
#Funções Auxiliares
##################################################

corrigeMixing <- function(g,mi){
  
  aux = T
  while(aux){
    v1 = sample(1:vcount(g),1)
    while (degree(g,v1) == maxdegree){
      v1 = sample(1:vcount(g),1)
    }
    vout = length(E(g)[v1 %--% V(g)[V(g)$p!=V(g)[v1]$p]])
    vtotal = degree(g,v1)
    if (vout/vtotal < (mi-toleranciamixing)){
      v2 = sample(1:vcount(g),1)
      while (degree(g,v2) == maxdegree){
        v2 = sample(1:vcount(g),1)
      }
      vout = length(E(g)[v2 %--% V(g)[V(g)$p!=V(g)[v2]$p]])
      vtotal = degree(g,v2)
      if (vout/vtotal < (mi-toleranciamixing)){
        g = add.edges(g, c(v1,v2))
        g = simplify(g)
      }     
    }
    
    if (calculaMixing(g) > (mi-toleranciamixing)){
      aux = F
    }
  }
  return(g)
}

densidadeComunidade <- function(g){
  dens = c()
  nc = unique(V(g)$p)
  for (i in 1:length(nc)){
    #print(graph.density(induced.subgraph(g,V(g)[V(g)$p==nc[i]])))
    dens = c(dens,graph.density(induced.subgraph(g,V(g)[V(g)$p==nc[i]])))
    
  }
  return(dens)
}
  
calculaMixing <- function(g){
  temp = 0
  for (i in 1:vcount(g)){
    vout = length(E(g)[i %--% V(g)[V(g)$p!=V(g)[i]$p]])
    vtotal = degree(g,i)
    temp = temp +vout/vtotal
  }
  
  if(temp == 0){
    return(mixing)
  }
  
  return(temp/vcount(g))
}

menorComunidade <- function(g){
  nc = unique(V(g)$p)
  temp = 0
  tamanho = vcount(g)
  for (i in nc){
    temp = length(V(g)$p[V(g)$p==i])
    if (temp<tamanho){
      tamanho = temp
    }
  }
  return(tamanho)
}

maiorComunidade <- function(g){
  nc = unique(V(g)$p)
  temp = 0
  tamanho = 0
  for (i in nc){
    temp = length(V(g)$p[V(g)$p==i])
    if (temp>tamanho){
      tamanho = temp
    }
  }
  return(tamanho)
}  

ilustrarGrafo <- function(g, nome = "", pasta = ""){
  nc = unique(V(g)$p)
  cores = rainbow(max(nc))
  for (i in nc){
    V(g)[V(g)$p == i]$color = cores[i]
  }
  
  arquivo = paste("g",toString(length(nc)),".png",sep="")
  
  if (nome != ""){
    if (pasta != ""){
      arquivo = paste(pasta,nome,".png",sep="")
    }else{
      arquivo = paste(nome,".png",sep="")
    }
  }else{
    if (pasta!=""){
      arquivo = paste(pasta,arquivo,sep="")
    }
  }
  png(filename=arquivo)
  plot.igraph(g)
  dev.off()
}
  
#################################################
#Testando Funções
#################################################

resultadosTestes = test_file("testDynamicalCommunities.R", reporter = "summary")

ntestes = length(resultadosTestes)
for(i in 1:ntestes){
  nresults = length(resultadosTestes[[i]]$results)
  for (j in 1:nresults){
    if (!resultadosTestes[[i]]$results[[j]]$passed){
      st = resultadosTestes[[i]]$results[[j]]$failure_msg
      st = str_replace(st,"\n",", ")
      aux = c(seed,st)
      write(aux,arquivoErro,ncolumns=2,append=T,sep="\t")
    }
  }
}

