
##################################################
#Funções Principais
##################################################


criarGrafoInicial <- function(p){
  arquivo = paste(p,"main.exe",sep="")
  setwd(p)
  if(msgDebug){
    cat("\nCriando Grafo Inicial")
  }
  system(arquivo,show.output.on.console = F)
  
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
    cat("\nTamanho da nova comunidade:",tamcomu)
    cat("\nVertice sendo adicionado: ")
  }
  for (i in 1:tamcomu){
    g = add.vertices(g,1)
    V(g)[vcount(g)]$p = idcomu
    if(msgDebug){
      cat(vcount(g),"")
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

  dc = densidadeComunidade(g)
  nc = sort(unique(V(g)$p))
  aux = which(nc==idcomu)
  densidademedia = mean(c(dc[-aux],edge_density(g)))
  densidademedia = max(densidademedia,edge_density(g))
  
  if (dc[aux] < densidademedia){
    g = corrigeDensidadeComu(g,idcomu)
  }
  
  if(mean(degree(g)) < (avgdegree - toleranciagrau)){
    g = corrigeGrauUp(g)
  }
  
  if(mean(degree(g)) > (avgdegree + toleranciagrau)){
    g = corrigeGrauDown(g)
  }

  if (calculaMixing(g) < (mi-toleranciamixing)){
    g = corrigeMixingUp(g,idcomu)
  }
  
  if (calculaMixing(g) > (mi+toleranciamixing)){
    g = corrigeMixingDown(g,idcomu)
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
  
  espaco = V(g)[V(g)$p==idcomu]
  if(msgDebug){
    cat("\nVertices a serem excluidos:",length(espaco))
    cat("\nVertice sendo excluído:")
  }
  
  while(length(espaco)>1){
    v1 = sample(espaco,1)
    if(msgDebug){
      cat("",v1)
    }
    g = delete.vertices(g,v1)
    espaco = V(g)[V(g)$p==idcomu]
  }
  if(msgDebug){
    cat("",espaco)
  }
 
  g = delete.vertices(g,espaco)
  
  if(mean(degree(g)) < (avgdegree - toleranciagrau)){
    g = corrigeGrauUp(g)
  }
  
  if(mean(degree(g)) > (avgdegree + toleranciagrau)){
    g = corrigeGrauDown(g)
  }
  
  if (calculaMixing(g) < (mixing-toleranciamixing)){
    g = corrigeMixingUp(g)
  }
  
  if (calculaMixing(g) > (mixing+toleranciamixing)){
    g = corrigeMixingDown(g)
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
    g = add.vertices(g,1)
    V(g)[vcount(g)]$p = idcomu
    if(msgDebug){
      cat(vcount(g),"")
    }
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
  dc = densidadeComunidade(g)
  nc = sort(unique(V(g)$p))
  aux = which(nc==idcomu)
  densidademedia = mean(c(dc[-aux],edge_density(g)))
  densidademedia = max(densidademedia,edge_density(g))
  
  if (dc[aux] < densidademedia){
    g = corrigeDensidadeComu(g,idcomu)
  }
  
  if(mean(degree(g)) < (avgdegree - toleranciagrau)){
    g = corrigeGrauUp(g)
  }
  
  if(mean(degree(g)) > (avgdegree + toleranciagrau)){
    g = corrigeGrauDown(g)
  }
  
  if (calculaMixing(g) < (mi-toleranciamixing)){
    g = corrigeMixingUp(g,idcomu)
  }
  
  if (calculaMixing(g) > (mi+toleranciamixing)){
    g = corrigeMixingDown(g,idcomu)
  }
  
  return(g)
}

contraction <- function(g,comu=0,nmin=minsize, mi = mixing){
  if(vcount(g)==0){
    return(g)
  }
  if (comu==0){
    idcomu = sample(V(g)$p,1)
  }else{
    idcomu = comu
  }
  tamcomuinicial = length(V(g)[V(g)$p==idcomu])
  while(tamcomuinicial<=nmin){
    if (comu==0){
      idcomu = sample(V(g)$p,1)
    }else{
      return(g)
    }
    tamcomuinicial = length(V(g)[V(g)$p==idcomu])
  }
  tamcomufinal = sample(nmin:(tamcomuinicial-1),1)
  velhosvertices = tamcomuinicial-tamcomufinal
  
  if (msgDebug){
    cat("\nVertices a excluir:",velhosvertices)
    cat("\nVertice sendo excluído: ")
  }
  
  for ( i in 1:velhosvertices){
    espaco = as.vector(V(g)[V(g)$p==idcomu])
    v1 = sample(espaco,1)
    if(msgDebug){
      cat("",v1)
    }
    g = delete.vertices(g,v1)
    
    
  }
  
  dc = densidadeComunidade(g)
  nc = sort(unique(V(g)$p))
  aux = which(nc==idcomu)
  densidademedia = mean(c(dc[-aux],edge_density(g)))
  densidademedia = max(densidademedia,edge_density(g))
  
  if (dc[aux] < densidademedia){
    g = corrigeDensidadeComu(g,idcomu)
  }
  
  if(mean(degree(g)) < (avgdegree - toleranciagrau)){
    g = corrigeGrauUp(g)
  }
  
  if(mean(degree(g)) > (avgdegree + toleranciagrau)){
    g = corrigeGrauDown(g)
  }
  
  if (calculaMixing(g) < (mi-toleranciamixing)){
    g = corrigeMixingUp(g,idcomu)
  }
  
  if (calculaMixing(g) > (mi+toleranciamixing)){
    g = corrigeMixingDown(g,idcomu)
  }
  
  return(g)
}

##################################################
#Funções Auxiliares
##################################################

corrigeMixingDown <-function(g,comu=0){
  if(msgDebug){
    cat("\nMixing Down: ")
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
    vizinhos = vizinhos[V(g)[vizinhos]$p!=idcomu]
    while (length(vizinhos)<=2){
      v1 = sample(aux,1)
      vizinhos = as.vector(neighbors(g,v1))
      vizinhos = vizinhos[V(g)[vizinhos]$p!=idcomu]
    }
    v2 = sample(vizinhos,1)
    while(degree(g,v2)<=2){
      v2 = sample(vizinhos,1)
    }
    
    aux = as.vector(V(g)[V(g)$p==idcomu])
    v3 = sample(aux,1)
    vizinhos2 = as.vector(neighbors(g,v3))
    vizinhos2 = vizinhos2[V(g)[vizinhos2]$p!=idcomu]
    while (length(vizinhos2)<=2){
      v3 = sample(aux,1)
      vizinhos2 = as.vector(neighbors(g,v3))
      vizinhos2 = vizinhos2[V(g)[vizinhos2]$p!=idcomu]
    }
    v4 = sample(vizinhos2,1)
    while(degree(g,v4)<=2){
      v4 = sample(vizinhos2,1)
    }
    
    aresta1 = get.edge.ids(g,c(v1,v2))
    aresta2 = get.edge.ids(g,c(v3,v4))
    g = delete.edges(g,aresta1)
    g = delete.edges(g,aresta2)
    g = add.edges(g,c(v1,v3))
    g = simplify(g)
    
    if(msgDebug){
      if (x != round(calculaMixing(g),2)){
        x = round(calculaMixing(g),2)
        cat(x,"")
      }
    }
    
    if(calculaMixing(g)<=(mixing)){
      flag = F
    }
    
  }
}

corrigeDensidadeComu <- function(g,idcomu){
  dc = densidadeComunidade(g)
  nc = sort(unique(V(g)$p))
  aux = which(nc==idcomu)
  flag = T
  densidademedia = mean(c(dc[-aux],edge_density(g)))
  densidademedia = max(densidademedia,edge_density(g))
  tamcomu = length(V(g)[V(g)$p==idcomu])
  graumaximo = min(tamcomu,maxdegree)
  if(msgDebug){
    cat("\nCorrigindo Densidade Comu:",idcomu)
    cat("\nDensidade da Rede:",edge_density(g))
    cat("\nDensidade Média:",densidademedia)
    cat("\nTamanho da Comu:",tamcomu)
    cat("\nGrau Maximo:",graumaximo)
    cat("\nDensidade Comu:")
  }
  x = 0
  while(dc[aux] < densidademedia && flag){
    espaco = V(g)[V(g)$p==idcomu]
    espaco = espaco[degree(g,espaco)<graumaximo]
    
    if(length(espaco)>1){
      v1 = sample(espaco,1)
      while(degree(g,v1)>=maxdegree){
        v1 = sample(espaco,1)
      }
      v2 = sample(espaco,1)
      while(degree(g,v2)>=maxdegree){
        v2 = sample(espaco,1)
      }
      
      g = add.edges(g,c(v1,v2))
      g = simplify(g)
    }else{
      if(msgDebug){
        cat("\nDensidade máxima atingida!")
      }
      flag = F
    }
    dc = densidadeComunidade(g)
    if(msgDebug){
      if (round(dc[aux],2) != x){
        x = round(dc[aux],2)
        cat("",x)
      }
    }
  }
  return(g)
}

corrigeMixingUp <- function(g,comu=0){
  if(msgDebug){
    cat("\nMixing Up: ")
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

corrigeGrauUp <- function(g){
  if(msgDebug){
    cat("\nGrau Up:")
  }
  x = 0
  while(mean(degree(g)) < min(avgdegree,vcount(g))){
    conexao = sample(c("in","out"),1,replace=F,prob=c(1-mixing,mixing))
    espaco = V(g)[degree(g,V(g))<maxdegree]
    v1 = sample(espaco,1)
    idcomu = V(g)[v1]$p
    if (conexao == "in"){
      espaco = V(g)[V(g)$p==idcomu]
    }else{
      if (conexao == "out"){
        espaco = V(g)[V(g)$p!=idcomu]
      }
    }
    espaco = espaco[degree(g,espaco)<maxdegree]
    v2 = sample(espaco,1)
    
    g = add.edges(g,c(v1,v2))
    g = simplify(g)
    if(msgDebug){
      if(x != round(mean(degree(g)))){
        x = round(mean(degree(g)))
        cat("",x)
      }
    }
  }
  return(g)
}

corrigeGrauDown <- function(g){
  if(msgDebug){
    cat("\nGrau Down:")
  }
  x = 0
  while(mean(degree(g)) > (avgdegree)){
    espaco = as.vector(V(g)[degree(g,V(g))<maxdegree])
    v1 = sample(espaco,1)
    idcomu = V(g)[v1]$p
    vizinhos = as.vector(neighbors(g,v1))
    if(sample(c(F,T),1)){
      vizinhos = vizinhos[V(g)[vizinhos]$p!=idcomu]
    }else{
      vizinhos = vizinhos[V(g)[vizinhos]$p==idcomu]
    }
    vizinhos = vizinhos[degree(g,vizinhos)>2]
    while(length(vizinhos)<2){
      v1 = sample(espaco,1)
      idcomu = V(g)[v1]$p
      vizinhos = as.vector(neighbors(g,v1))
      if(sample(c(F,T),1)){
        vizinhos = vizinhos[V(g)[vizinhos]$p!=idcomu]
      }else{
        vizinhos = vizinhos[V(g)[vizinhos]$p==idcomu]
      }
      vizinhos = vizinhos[degree(g,vizinhos)>2]
    }
    v2 = sample(vizinhos,1)
    
    aresta = get.edge.ids(g,c(v1,v2))
    g = delete.edges(g,aresta)
    g = simplify(g)
    if(msgDebug){
      if(x != round(mean(degree(g)))){
        x = round(mean(degree(g)))
        cat("",x)
      }
    }
  }
  return(g)
}

densidadeComunidade <- function(g){
  dens = c()
  nc = sort(unique(V(g)$p))
  for (i in 1:length(nc)){
    dens = c(dens,graph.density(induced.subgraph(g,V(g)[V(g)$p==nc[i]])))
  }
  return(dens)
}
  
calculaMixing <- function(g){
  nc = unique(V(g)$p)
  if (length(nc)>=2){
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
  if(vcount(g)==0){
    return(minsize)
  }
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
  if(vcount(g)==0){
    return(maxsize)
  }
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

