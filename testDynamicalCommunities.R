

grafoInicial = criarGrafoInicial(path)
cat("\n")
test_that("Função criarGrafoInicial",{
  
  
  expect_that(is.igraph(grafoInicial), is_true())
  
  expect_that(vcount(grafoInicial),equals(nvertices))
  expect_that(mean(degree(grafoInicial)), equals(avgdegree, tolerance=toleranciagrau, scale=1))
  expect_that(max(degree(grafoInicial)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(grafoInicial)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(grafoInicial), is_more_than(minsize-1))
  expect_that(maiorComunidade(grafoInicial), is_less_than(maxsize+1))
  
})

g = grafoInicial
nfuncoes = sample(4:8,1)

for(i in 1:nfuncoes){
  funcao = sample(c("b","e","g","c"),1,prob=probfuncao)
  if (funcao == "b"){
    if(msgDebug){
      cat("\nFunção",i,": Born","\n")
    }
    nv = vcount(g)
    g = born(g,nmax=round(maxsize/3))
    cat("\n")
    test_that("Função born",{
      idcomu = max(V(g)$p)
      tamcomu = length(V(g)$p[V(g)$p==idcomu])
      
      expect_gt(vcount(g),nv)
      expect_that(mean(degree(g)), equals(avgdegree, tolerance=toleranciagrau, scale=1))
      expect_that(max(degree(g)),is_less_than(maxdegree+1))
      mixingReal = calculaMixing(g)
      expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
      expect_that(menorComunidade(g), is_more_than(minsize-1))
      expect_that(maiorComunidade(g), is_less_than(maxsize+1))
      
      
    })
  }
  
  if(funcao == "e"){
    if(msgDebug){
      cat("\nFunção",i,": Extinction","\n")
    }
    nv = vcount(g)
    if(nv!=0){
      idcomu = sample(V(g)$p,1)
      tamcomu = length(V(g)$p[V(g)$p==idcomu])
      g = extinction(g, comu = idcomu)
      if(vcount(g)!=0){
        cat("\n")
        test_that("Função extinction",{
          
          expect_lt(vcount(g),nv)
          nv2 = length(V(g)$p[V(g)$p==idcomu])
          expect_that(nv2,equals(0))
          expect_that(mean(degree(g)), equals(avgdegree, tolerance=toleranciagrau, scale=1))
          expect_that(max(degree(g)),is_less_than(maxdegree+1))
          mixingReal = calculaMixing(g)
          expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
          expect_that(menorComunidade(g), is_more_than(minsize-1))
          expect_that(maiorComunidade(g), is_less_than(maxsize+1))
          
        })
      }
    }else{
      g = extinction(g)
    }
    
    
  }
  
  if(funcao == "g"){
    if(msgDebug){
      cat("\nFunção",i,": Growth","\n")
    }
    nv = vcount(g)
    if(nv!=0){
      idcomu = sample(V(g)$p,1)
      tamcomu = length(V(g)$p[V(g)$p==idcomu])
      
      g = growth(g,comu=idcomu)
      cat("\n")
      test_that("Função growth",{
        
        expect_gte(vcount(g),nv)
        nv2 = length(V(g)$p[V(g)$p==idcomu])
        expect_gte(nv2,tamcomu)
        expect_that(mean(degree(g)), equals(avgdegree, tolerance=toleranciagrau, scale=1))
        expect_that(max(degree(g)),is_less_than(maxdegree+1))
        mixingReal = calculaMixing(g)
        expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
        expect_that(menorComunidade(g), is_more_than(minsize-1))
        expect_that(maiorComunidade(g), is_less_than(maxsize+1))
      })
    }else{
      g = growth(g)
    }
    
  }
  
  if(funcao == "c"){
    if(msgDebug){
      cat("\nFunção",i,": Contraction","\n")
    }
    nv = vcount(g)
    if(nv!=0){
      idcomu = sample(V(g)$p,1)
      tamcomu = length(V(g)[V(g)$p==idcomu])
      
      g = contraction(g,comu=idcomu)
      cat("\n")
      test_that("Função contraction",{
        
        expect_lte(vcount(g),nv)
        nv2 = length(V(g)$p[V(g)$p==idcomu])
        expect_lte(nv2,tamcomu)
        expect_that(mean(degree(g)), equals(avgdegree, tolerance=toleranciagrau, scale=1))
        expect_that(max(degree(g)),is_less_than(maxdegree+1))
        mixingReal = calculaMixing(g)
        expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
        expect_that(menorComunidade(g), is_more_than(minsize-1))
        expect_that(maiorComunidade(g), is_less_than(maxsize+1))
      })
    }else{
      g = contraction(g)
    }
    
  }
  
  
  if(msgDebug){
    cat("\n_____________________________________________\n")
    cat("Resumo Grafo Atual:\n")
    cat("Número de Vértices:",vcount(g),"\n")
    cat("Número de Comunidades:",length(unique(V(g)$p)),"\n")
    cat("Grau médio:",mean(degree(g)),"\n")
    cat("Grau máximo:",max(degree(g)),"\n")
    cat("Densidades:",densidadeComunidade(g),"\n")
    cat("Mixing:",calculaMixing(g),"\n")
    cat("Menor Comu:",menorComunidade(g),"\n")
    cat("Maior Comu:",maiorComunidade(g),"\n")
    cat("_____________________________________________\n")
  }
}









