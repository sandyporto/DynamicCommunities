library(igraph)


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

test_that("",{
  grafoInicial = criarGrafoInicial()
  
  expect_that(is.igraph(grafoInicial), is_true())
  
  expect_that(vcount(grafoInicial),equals(300))
  avgdegree = 20
  expect_that(mean(degree(grafoInicial)), equals(avgdegree, tolerance=1))
  maxdegree = 40
  expect_that(max(degree(grafoInicial)),equals(maxdegree))
  mixingEsperado = 0.05
  mixingReal = calculaMixing(grafoInicial)
  expect_that(mixingReal, equals(mixingEsperado, tolerance = 0.05))
  expect_that(menorComunidade(grafoInicial), is_more_than(40-1))
  expect_that(maiorComunidade(grafoInicial), is_less_than(80+1))
  
})