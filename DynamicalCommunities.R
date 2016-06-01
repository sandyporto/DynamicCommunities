
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
  V(G)$mi = 0
  G = atualizaMixing(G)
  
  if (maiorComunidade(G)>maxsize 
      || mean(degree(G))>(avgdegree+toleranciagrau)
      || mean(degree(G))<(avgdegree-toleranciagrau)){
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
    g = adicionaVertice(g,idcomu)
    if(msgDebug){
      cat(vcount(g),"")
    }
    
    if (i==1){
      if (taminicial != 0){
        vlist = c(vcount(g),sample(1:(vcount(g)-1),1))
        g = adicionaAresta(g,vlist)
      }
    }else{
      auxgrau = min((i+1),avgdegree)
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
          
          vlist = c(v1,v2)
          g = adicionaAresta(g,vlist)
          
        }
        if (conexao == "in"){
          espaco = as.vector(V(g)[V(g)$p==idcomu])
          v1 = vcount(g)
          v2 = sample(espaco,1)
          
          while(degree(g,v2) == maxdegree){
            v2 = sample(espaco,1)
          }
          
          vlist = c(v1,v2)
          g = adicionaAresta(g,vlist)
          
        }
      }
    }
    g = corrigeRede(g)
  }
  
  g = corrigeRede(g,idcomu)
  
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
    g = deletaVertice(g,v1)
    espaco = V(g)[V(g)$p==idcomu]
    
  }
  if(msgDebug){
    cat("",espaco)
  }
 
  g = deletaVertice(g,espaco)
  
  if(vcount(g)==0){
    return(g)
  }
  
  g = corrigeRede(g)
  
  return(g)
  
}

growth <- function(g, comu = 0, nmax = maxsize, mi = mixing){
  if(vcount(g)==0){
    return(g)
  }
  if (comu==0){
    listidcomu = sample(unique(V(g)$p),replace=F)
  }else{
    listidcomu = comu
  }
  tamcomuinicial = 0
  for(idcomu in listidcomu){
    if(length(V(g)[V(g)$p==idcomu])>=nmax){next()}
    tamcomuinicial = length(V(g)[V(g)$p==idcomu])
    break()
  }
  if(tamcomuinicial==0){return(g)}
  if(tamcomuinicial<nmax){
    tamcomufinal = sample((tamcomuinicial):nmax,1)
  }else{
    tamcomufinal = tamcomuinicial+1
  }
  
  novosvertices = tamcomufinal-tamcomuinicial
  if(novosvertices==0){return(g)}
  
  if (msgDebug){
    cat("\nVertices a adicionar:",novosvertices)
    cat("\nVertice sendo adicionado: ")
  }
  
  for (i in 1:novosvertices){
    g = adicionaVertice(g,idcomu)
    if(msgDebug){
      cat(vcount(g),"")
    }
    auxgrau = min((tamcomuinicial+i),avgdegree)
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
        
        while(degree(g,v2) >= maxdegree){
          v2 = sample(espaco,1)
        }
        
        vlist = c(v1,v2)
        g = adicionaAresta(g,vlist)
        
      }
      if (conexao == "in"){
        espaco = as.vector(V(g)[V(g)$p==idcomu])
        v1 = vcount(g)
        v2 = sample(espaco,1)
        
        while(degree(g,v2) >= maxdegree){
          v2 = sample(espaco,1)
        }
        
        vlist = c(v1,v2)
        g = adicionaAresta(g,vlist)
        
      }
    }
    g = corrigeRede(g)
  }
  
  g = corrigeRede(g,idcomu)
  
  return(g)
}

contraction <- function(g,comu=0,nmin=minsize, mi = mixing){
  if(vcount(g)==0){
    return(g)
  }
  if (comu==0){
    listidcomu = sample(unique(V(g)$p),replace=F)
  }else{
    listidcomu = comu
  }
  tamcomuinicial = 0
  for(idcomu in listidcomu){
    if(length(V(g)[V(g)$p==idcomu])<=nmin){next()}
    tamcomuinicial = length(V(g)[V(g)$p==idcomu])
    break()
  }
  
  if(tamcomuinicial==0){return(g)}
  if(tamcomuinicial>=(nmin)){
    tamcomufinal = sample(nmin:(tamcomuinicial),1)
  }else{
    tamcomufinal = tamcomuinicial-1
  }
  
  velhosvertices = tamcomuinicial-tamcomufinal
  if(velhosvertices==0){return(g)}
  
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
    g = deletaVertice(g,v1)
    g = corrigeRede(g)
  }
  
  g = corrigeRede(g,idcomu)
  
  return(g)
}

##################################################
#Funções Auxiliares
##################################################

corrigeRede <- function(g,idcomu=0){
  if(idcomu!=0){
    dc = densidadeComunidade(g)
    nc = sort(unique(V(g)$p))
    aux = which(nc==idcomu)
    
    if(length(dc)>=2){
      densidademedia = mean(dc[-aux])
    }else{
      densidademedia = edge_density(g)
    }
    densidademedia = max(densidademedia,edge_density(g))
    
    tamcomu = (length(V(g)[V(g)$p==idcomu]))
    graumaximo = min((tamcomu-1),maxdegree)
    if(dc[aux]<densidademedia){
      g = corrigeDensidadeComu(g,idcomu)
    }
  }
  
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

corrigeMixingDown <-function(g){
  if(msgDebug){
    cat("\nMixing Down: ")
  }
  
  x=0
  times = 0
  maxtimes = vcount(g)
  while(calculaMixing(g)>mixing && times<maxtimes){
    idcomuv1 = sample(V(g)$p,1)
    espaco = as.vector(V(g)[V(g)$p==idcomuv1])
    espaco = espaco[degree(g,espaco)>2]
    
    v2ok=F
    if(length(espaco)<2){next()}
    for (v1 in sample(espaco)){
      vizinhos1 = as.vector(neighbors(g,v1))
      vizinhos1 = vizinhos1[V(g)[vizinhos1]$p!=idcomuv1]
      vizinhos1 = vizinhos1[degree(g,vizinhos1)>2]
      if(length(vizinhos1)<2){next()}
      v2 = sample(vizinhos1,1)
      v2ok=T
      break()
    }
    if(!v2ok){times=times+1;next()}
    idcomuv2 = V(g)[v2]$p
    vizinhos2 = as.vector(neighbors(g,v2))
    
    espaco = espaco[espaco!=v1]
    vizinhos1 = as.vector(neighbors(g,v1))
    espaco = espaco[!espaco %in% vizinhos1]
    if(length(espaco)<2){next()}
    v4ok=F
    for(v3 in sample(espaco)){
      vizinhos3 = as.vector(neighbors(g,v3))
      vizinhos3 = vizinhos3[V(g)[vizinhos3]$p==idcomuv2]
      vizinhos3 = vizinhos3[!vizinhos3 %in% vizinhos2]
      vizinhos3 = vizinhos3[degree(g,vizinhos3)>2]
      
      if(length(vizinhos3)<2){next()}
      v4 = sample(vizinhos3,1)
      v4ok=T
      break()
    }
    if(!v4ok){times=times+1;next()}
    
    vlist = c(v1,v2)
    g = deletaAresta(g,vlist)
    vlist = c(v3,v4)
    g = deletaAresta(g,vlist)
    vlist = c(v1,v3)
    g = adicionaAresta(g,vlist)
    vlist = c(v2,v4)
    g = adicionaAresta(g,vlist)
    
    if(msgDebug){
      if (x != round(calculaMixing(g),2)){
        x = round(calculaMixing(g),2)
        cat(x,"")
      }
    }
  }
  if (times==maxtimes){
    if(msgDebug){
      cat("\nNão é mais possível diminuir o mixing!")
    }
  }
  cat("\n")
  return(g)
}

corrigeDensidadeComu <- function(g,idcomu){
  dc = densidadeComunidade(g)
  nc = sort(unique(V(g)$p))
  aux = which(nc==idcomu)
  
  densidademedia = mean(dc[-aux])
  densidademedia = max(densidademedia,edge_density(g))
  tamcomu = (length(V(g)[V(g)$p==idcomu]))
  graumaximo = min((tamcomu-1),maxdegree)
  
  if(msgDebug){
    cat("\nDensidade Média:",densidademedia)
    cat("\nDensidade Comu",idcomu,":")
  }

  x = 0
  while(dc[aux] < densidademedia){
    espaco = as.vector(V(g)[V(g)$p==idcomu])
    espaco = espaco[degree(g,espaco)<graumaximo]
    #cat("\nEspaco:",espaco)
    if(length(espaco)<2){
      if(msgDebug){
        cat("\nDensidade máxima atingida!")
      }
      break()
    }
    v2 = 0
    for(v1 in sample(espaco,replace=F)){
      vizinhos = as.vector(neighbors(g,v1))
      #cat("\nVizinhos",vizinhos)
      aux2 = match(c(v1,vizinhos),espaco,nomatch = 0)
      aux2 = aux2[aux2>0]
      espaco2 = espaco[-aux2]
      #cat("\nEspaco2:",espaco2)
      if(length(espaco2)<2){next()}
      v2 = sample(espaco2,1)
      break()
    }
    if(v2==0){
      if(msgDebug){
        cat("\nDensidade máxima atingida!")
      }
      break()
    }
    
    #cat("\nV1:",v1,", v2:",v2)
    vlist = c(v1,v2)
    g = adicionaAresta(g,vlist)
    #cat("\necount:",ecount(g))
    dc = densidadeComunidade(g)
    if(msgDebug){
      if (round(dc[aux],2) != x){
        x = round(dc[aux],2)
        cat("",x)
      }
    }
  }
  cat("\n")
  return(g)
}

corrigeMixingUp <- function(g){
  if(msgDebug){
    cat("\nMixing Up: ")
  }
  
  x=0
  times = 0
  maxtimes = vcount(g)
  while(calculaMixing(g)<mixing && times<maxtimes){
    idcomuv1 = sample(V(g)$p,1)
    espaco = as.vector(V(g)[V(g)$p==idcomuv1])
    espaco = espaco[degree(g,espaco)>2]                  
    
    v2ok = F
    if(length(espaco)<2){times=times+1;next()}
    for(v1 in sample(espaco)){
      vizinhos1 = as.vector(neighbors(g,v1))
      vizinhos1 = vizinhos1[V(g)[vizinhos1]$p == idcomuv1]
      vizinhos1 = vizinhos1[degree(g,vizinhos1)>2]
      if(length(vizinhos1)<2){next()}
      for(v2 in sample(vizinhos1)){
        idcomuv2 = V(g)[v2]$p
        vizinhos2 = as.vector(neighbors(g,v2))
        vizinhos2 = vizinhos2[V(g)[vizinhos2]$p==idcomuv2]
        if(length(vizinhos2)<3){next()}
        v2ok=T
        break()
      }
      break()
    }
    if(!v2ok){times=times+1;next()}
    vizinhos2 = as.vector(neighbors(g,v2))
    
    espaco = as.vector(V(g)[V(g)$p!=idcomuv1])
    vizinhos1 = as.vector(neighbors(g,v1))
    espaco = espaco[!espaco %in% vizinhos1]
    if(length(espaco)<2){next()}
    v4ok=F
    for(v3 in sample(espaco)){
      vizinhos3 = as.vector(neighbors(g,v3))
      idcomuv3 = V(g)[v3]$p
      vizinhos3 = vizinhos3[V(g)[vizinhos3]$p==idcomuv3]
      vizinhos3 = vizinhos3[!vizinhos3 %in% vizinhos2]
      vizinhos3 = vizinhos3[degree(g,vizinhos3)>2]
      if(length(vizinhos3)<2){next()}
      for(v4 in sample(vizinhos3)){
        idcomuv4 = V(g)[v4]$p
        vizinhos4 = as.vector(neighbors(g,v4))
        vizinhos4 = vizinhos4[V(g)[vizinhos4]$p==idcomuv4]
        if(length(vizinhos4)<3){next()}
        v4ok=T
        break()
      }
      break()
    }
    if(!v4ok){times=times+1;next()}
    
    vlist = c(v1,v2)
    g = deletaAresta(g,vlist)
    vlist = c(v3,v4)
    g = deletaAresta(g,vlist)
    vlist = c(v1,v3)
    g = adicionaAresta(g,vlist)
    vlist = c(v2,v4)
    g = adicionaAresta(g,vlist)
    
    if(msgDebug){
      if (x != round(calculaMixing(g),2)){
        x = round(calculaMixing(g),2)
        cat(x,"")
      }
    }
  }
  if (times==maxtimes){
    if(msgDebug){
      cat("\nNão é mais possível aumentar o mixing!")
    }
  }
  cat("\n")
  return(g)
}

corrigeGrauUp <- function(g){
  if(msgDebug){
    cat("\nGrau Up:")
  }
  x = 0
  times = 0
  maxtimes = vcount(g)
  while(mean(degree(g)) < min(avgdegree,(vcount(g)-1)) 
        && times<maxtimes){
    espaco = as.vector(V(g)[degree(g,V(g))<maxdegree])
    v2ok=F
    for(v1 in sample(espaco,replace=F)){
      idcomu = V(g)[v1]$p
      conexao = sample(c("in","out"),1,replace=F,prob=c(1-mixing,mixing))
      if(length(unique(V(g)$p))<2){conexao="in"}
      if (conexao == "in"){
        espaco2 = as.vector(V(g)[V(g)$p==idcomu])
      }else{
        if (conexao == "out"){
          espaco2 = as.vector(V(g)[V(g)$p!=idcomu])
        }
      }
      espaco2 = espaco2[degree(g,espaco2)<maxdegree]
      vizinhos = as.vector(neighbors(g,v1))
      espaco2 = espaco2[!espaco2 %in% vizinhos]
      
      if(length(espaco2)<2){next()}
      v2 = sample(espaco2,1)
      v2ok=T
      break()
    }
    if(!v2ok){times=times+1;next()}
    
    vlist = c(v1,v2)
    g = adicionaAresta(g,vlist)
    
    if(msgDebug){
      if(x != round(mean(degree(g)))){
        x = round(mean(degree(g)))
        cat("",x)
      }
    }
  }
  if (times==maxtimes){
    if(msgDebug){
      cat("\nNão é mais possível aumentar o grau!")
    }
  }
  cat("\n")
  return(g)
}

corrigeGrauDown <- function(g){
  if(msgDebug){
    cat("\nGrau Down:")
  }
  x = 0
  times = 0
  maxtimes = vcount(g)
  while(mean(degree(g)) > (avgdegree)
        && times<maxtimes){
    espaco = as.vector(V(g)[degree(g,V(g))<maxdegree])
    v2ok=F
    for(v1 in sample(espaco,replace=F)){
      idcomu = V(g)[v1]$p
      vizinhos = as.vector(neighbors(g,v1))
      if(sample(c(T,F),1,prob=c(1-mixing,mixing))){
        vizinhos = vizinhos[V(g)[vizinhos]$p==idcomu]
      }else{
        vizinhos = vizinhos[V(g)[vizinhos]$p!=idcomu]
      }
      vizinhos = vizinhos[degree(g,vizinhos)>2]
      if(length(vizinhos)<2){next()}
      v2=sample(vizinhos,1)
      v2ok=T
      break()
    }
    if(!v2ok){times=times+1;next()}
    
    vlist = c(v1,v2)
    g = deletaAresta(g,vlist)
    
    if(msgDebug){
      if(x != round(mean(degree(g)))){
        x = round(mean(degree(g)))
        cat("",x)
      }
    }
  }
  if (times==maxtimes){
    if(msgDebug){
      cat("\nNão é mais possível diminuir o grau!")
    }
  }
  cat("\n")
  return(g)
}

densidadeComunidade <- function(g){
  dens = c()
  nc = sort(unique(V(g)$p))
  for (i in 1:length(nc)){
    dens = c(dens,edge_density(induced_subgraph(g,V(g)[V(g)$p==nc[i]],impl="auto")))
  }
  return(dens)
}

calculaMixing <- function(g){
  nc = unique(V(g)$p)
  if(length(nc)>1){
    if(menorComunidade(g)<minsize){
      return(mixing)
    }
    return(sum(V(g)$mi)/vcount(g))
  }else{
    return(mixing)
  }
}

atualizaMixing <- function(g, vlist = V(g)){
  eout <- function(g,v){
    idcomu = V(g)[v]$p
    vizinhos = as.vector(neighbors(g,v))
    vizinhos = vizinhos[V(g)[vizinhos]$p!=idcomu]
    eout = length(vizinhos)
    return(eout)
  }
  for (v in vlist){
    nedgeout = eout(g,v)
    if(degree(g,v)!=0){
      V(g)[v]$mi = nedgeout/degree(g,v)
    }else{
      V(g)[v]$mi = 0
    }
  }
  return(g)
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

adicionaVertice <- function(g,idcomu){
  g = add.vertices(g,1)
  v = vcount(g)
  V(g)[v]$p = idcomu
  V(g)[v]$mi = 0
  gravaGrafo(g)
  return(g)
}

adicionaAresta <- function(g,vlist){
  g = add.edges(g,vlist)
  g = simplify(g)
  g = atualizaMixing(g,vlist)
  gravaGrafo(g)
  return(g)
}

deletaVertice <- function(g,v){
  vizinhos = as.vector(neighbors(g,v))
  g = delete.vertices(g,v)
  if (length(vizinhos)>0){
    aux = which(vizinhos>v)
    vizinhos[aux]=vizinhos[aux]-1
    g = atualizaMixing(g,vizinhos)
  }
  gravaGrafo(g)
  return(g)
}

deletaAresta <- function(g,vlist){
  aresta = get.edge.ids(g,vlist)
  g = delete.edges(g,aresta)
  g = atualizaMixing(g,vlist)
  return(g)
}

gravaGrafo <- function(g){
  if (listaGrafos){
    pathGraph = paste(pastaGrafos,"g",ngraph,".dat",sep="")
    write_graph(g,pathGraph,format = "edgelist")
    ngraph = ngraph+1
  }
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
