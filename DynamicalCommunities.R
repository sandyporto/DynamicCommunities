
##################################################
#Funções Principais
##################################################


criarGrafoInicial <- function(p){
  arquivo = paste(p,"main.exe",sep="")
  setwd(p)
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
    G = criarGrafoInicial(p)
  }
  setwd(pasta)
  return(G)
}

born <- function(g, nmin = minsize, nmax = maxsize, mi=mixing){
  taminicial = vcount(g)
  tamcomu = sample(nmin:nmax,1)
  if (taminicial != 0){
    idcomu = max(V(g)$p)+1
  }else{
    idcomu = 1
  }
  if(msgDebug){
    cat("\nTamanho da comunidade:",tamcomu)
    cat("\nVertice sendo adicionado: ")
  }
  for (i in 1:tamcomu){
    g = add.vertices(g,1)
    V(g)[vcount(g)]$p = idcomu
    if(msgDebug){
      cat(i,"")
    }
    
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
  
  if(msgDebug){
    cat("\nDens:",dens)
    cat("\nMeanDens:",mean(dens))
    cat("\nCMAD:",comumaxavgdegree)
    cat("\navg/tam:",avgdegree/tamcomu)
  }
  
  d = min(mean(dens),comumaxavgdegree)
  d = max(d,(avgdegree)/tamcomu)
  flag = T
  if(msgDebug){
    cat("\nValor d:",d)
    cat("\nGrau Médio: ")
  }
  x = 0
  while(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])) < (d) && flag){
    espaco = as.vector(V(g)[V(g)$p==idcomu])
    espaco = espaco[degree(g,espaco) < maxdegree]
    if(msgDebug){
      if ( x != round(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])),2)){
        x = round(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])),2)
        cat(x,"")
      }
    }
    
    if (length(espaco)>=2){
      v1 = sample(espaco,1)
      v2 = sample(espaco,1)
      
      g = add.edges(g,c(v1,v2))
      g = simplify(g)
    }else{
      if(msgDebug){
        cat("\nEspaco pequeno:",espaco,"\n")
      }
      
      flag = F
    }
    
  }

  if (calculaMixing(g) < (mi)){
    if(msgDebug){
      cat("\nCorrigindo")
    }
    
    g = corrigeMixing(g,idcomu)
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
  
  if (calculaMixing(g) < (mixing)){
    if(msgDebug){
      cat("\nCorrigindo")
    }
    g = corrigeMixing(g)
  }
  
  return(g)
}

growth <- function(g, comu = 0, nmax = maxsize, mi = mixing){
  if(vcount(g)==0){
    return(g)
  }
  if (comu==0){
    idcomu = sample(V(g)$p,1)
  }else{
    idcomu = comu
  }
  tamcomuinicial = length(V(g)[V(g)$p==idcomu])
  while(tamcomuinicial>=nmax){
    if (comu==0){
      idcomu = sample(V(g)$p,1)
    }else{
      return(g)
    }
    tamcomuinicial = length(V(g)[V(g)$p==idcomu])
  }
  tamcomufinal = sample(tamcomuinicial:nmax,1)
  novosvertices = tamcomufinal-tamcomuinicial
  
  if (msgDebug){
    cat("\nVertices a adicionar:",novosvertices)
    cat("\nVertice sendo adicionado: ")
  }
  
  for (i in 1:novosvertices){
    if(msgDebug){
      cat(i,"")
    }
    
    g = add.vertices(g,1)
    V(g)[vcount(g)]$p = idcomu
    auxgrau = min((tamcomuinicial+i),maxdegree)
    grau = sample(2:auxgrau,1)
    for (j in 1:grau){
      conexao = sample(c("in","out"),1,replace=F,c(1-mi,mi))
      
      if(length(unique(V(g)$p)) == 1){
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
  tamcomu = length(V(g)[V(g)$p==idcomu])
  dens = densidadeComunidade(g)
  narestasout = length(E(g)[V(g)[V(g)$p==idcomu] %--% V(g)[V(g)$p!=idcomu]])
  comumaxavgdegree = (((maxdegree*tamcomu)/2) - narestasout)/(tamcomu*(tamcomu-1)/2)
  
  if(msgDebug){
    cat("\nDens:",dens)
    cat("\nMeanDens:",mean(dens))
    cat("\nCMAD:",comumaxavgdegree)
    cat("\navg/tam:",avgdegree/tamcomu)
  }
  
  d = min(mean(dens),comumaxavgdegree)
  d = max(d,(avgdegree)/tamcomu)
  flag = T
  if(msgDebug){
    cat("\nValor d:",d)
    cat("\nGrau Médio: ")
  }
  x=0
  while(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])) < (d) && flag){
    espaco = as.vector(V(g)[V(g)$p==idcomu])
    espaco = espaco[degree(g,espaco) < maxdegree]
    if(msgDebug){
      if ( x != round(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])),2)){
        x = round(graph.density(induced.subgraph(g,V(g)[V(g)$p==idcomu])),2)
        cat(x,"")
      }
    }
    
    if (length(espaco)>=2){
      v1 = sample(espaco,1)
      v2 = sample(espaco,1)
      
      g = add.edges(g,c(v1,v2))
      g = simplify(g)
    }else{
      if(msgDebug){
        cat("\nEspaco pequeno:",espaco,"\n")
      }
      
      flag = F
    }
    
  }
  
  if (calculaMixing(g) < (mi)){
    if(msgDebug){
      cat("\nCorrigindo")
    }
    
    g = corrigeMixing(g,idcomu)
  }
  
  return(g)
}


##################################################
#Funções Auxiliares
##################################################

corrigeMixing <- function(g,comu=0){
  if(msgDebug){
    cat("\nMixing: ")
  }
  
  flag = T
  x=0
  while(flag){
    if (comu!=0){
      idcomu=comu
    }else{
      idcomu = sample(V(g)$p,1)
    }
    aux = as.vector(V(g)[V(g)$p==idcomu])
    v1 = sample(aux,1)
    vizinhos = as.vector(neighbors(g,v1))
    vizinhos = vizinhos[V(g)[vizinhos]$p==idcomu]
    while (length(vizinhos)<=2){
      v1 = sample(aux,1)
      vizinhos = as.vector(neighbors(g,v1))
      vizinhos = vizinhos[V(g)[vizinhos]$p==idcomu]
    }
    v2 = sample(vizinhos,1)
    vizinhos2 = as.vector(neighbors(g,v2))
    vizinhos2 = vizinhos2[V(g)[vizinhos2]$p==idcomu]
    while(length(vizinhos2)<=2){
      v2 = sample(vizinhos,1)
      vizinhos2 = as.vector(neighbors(g,v2))
      vizinhos2 = vizinhos2[V(g)[vizinhos2]$p==idcomu]
    }
    
    aux = as.vector(V(g)[V(g)$p!=idcomu])
    v3 = sample(aux,1)
    while(degree(g,v3)>= maxdegree){
      v3 = sample(aux,1)
    }
    v4 = sample(aux,1)
    while(degree(g,v4)>= maxdegree){
      v4 = sample(aux,1)
    }
    
    aresta = get.edge.ids(g,c(v1,v2))
    g = delete.edges(g,aresta)
    g = add.edges(g,c(v1,v3))
    g = add.edges(g,c(v2,v4))
    g = simplify(g)
    
    if(msgDebug){
      if (x != round(calculaMixing(g),2)){
        x = round(calculaMixing(g),2)
        cat(x,"")
      }
    }
    
    if(calculaMixing(g)>=(mixing)){
      flag = F
    }
  }
  
  return(g)
}

densidadeComunidade <- function(g){
  dens = c()
  nc = unique(V(g)$p)
  for (i in 1:length(nc)){
    dens = c(dens,graph.density(induced.subgraph(g,V(g)[V(g)$p==nc[i]])))
    
  }
  return(dens)
}
  
calculaMixing <- function(g){
  nc = unique(V(g)$p)
  if (length(nc)>2){
    temp = 0
    for (i in nc){
      vout = length(E(g)[V(g)[V(g)$p==i] %--% V(g)[V(g)$p!=i]])
      vtotal = sum(degree(g)[V(g)$p==i])
      temp = temp +vout/vtotal
    }
    return(temp/length(nc))
  }else{
    return(mixing)
  }
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

if (testeCompleto){
  resultadosTestes2 = test_file("testDynamicalCommunities2.R", reporter = "summary")
  ntestes = length(resultadosTestes2)
  for(i in 1:ntestes){
    nresults = length(resultadosTestes2[[i]]$results)
    for (j in 1:nresults){
      if (!resultadosTestes2[[i]]$results[[j]]$passed){
        st = resultadosTestes2[[i]]$results[[j]]$failure_msg
        st = str_replace(st,"\n",", ")
        aux = c(seed,st)
        write(aux,arquivoErro,ncolumns=2,append=T,sep="\t")
      }
    }
  }
}

