library(igraph)


test_that("Grafo vazio para função born",{
  g = make_empty_graph(0,F)
  V(g)$p = 0
  g = born(g)
  nc = length(unique(V(g)$p[V(g)$p!=0]))
  expect_that(nc, equals(1))
  expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(g)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(g)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(g), is_more_than(minsize-1))
  expect_that(maiorComunidade(g), is_less_than(maxsize+1))
})

test_that("Grafo vazio para função extinction",{
  g = make_empty_graph(0,F)
  g = extinction(g)
  nv = vcount(g)
  expect_that(nv, equals(0))
})


grafoInicial = criarGrafoInicial()
#ilustrarGrafo(grafoInicial)


test_that("Função criarGrafoInicial",{
  
  
  expect_that(is.igraph(grafoInicial), is_true())
  
  expect_that(vcount(grafoInicial),equals(nvertices))
  expect_that(mean(degree(grafoInicial)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(grafoInicial)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(grafoInicial)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(grafoInicial), is_more_than(minsize-1))
  expect_that(maiorComunidade(grafoInicial), is_less_than(maxsize+1))
  
})

g = born(grafoInicial)

test_that("Função born",{
  #grafoInicial = criarGrafoInicial()
  
  
  
  idcomu = max(V(g)$p)
  tamcomu = length(V(g)$p[V(g)$p==idcomu])
  
  expect_that(vcount(g),equals(nvertices+tamcomu))
  expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(g)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(g)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(g), is_more_than(minsize-1))
  expect_that(maiorComunidade(g), is_less_than(maxsize+1))
  
  
})

# g = born(g)
# g = born(g)
# 
# test_that("Usando a função born 3 vezes",{
#   idcomu = max(V(g)$p)
#   tamcomu = length(V(g)$p[V(g)$p==idcomu])
#   
#   expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
#   expect_that(max(degree(g)),is_less_than(maxdegree+1))
#   mixingReal = calculaMixing(g)
#   expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
#   expect_that(menorComunidade(g), is_more_than(minsize-1))
#   expect_that(maiorComunidade(g), is_less_than(maxsize+1))
# })
# 
# g = born(g)
# g = born(g)


# test_that("Usando a função born 5 vezes",{
#   idcomu = max(V(g)$p)
#   tamcomu = length(V(g)$p[V(g)$p==idcomu])
#   
#   expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
#   expect_that(max(degree(g)),is_less_than(maxdegree+1))
#   mixingReal = calculaMixing(g)
#   expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
#   expect_that(menorComunidade(g), is_more_than(minsize-1))
#   expect_that(maiorComunidade(g), is_less_than(maxsize+1))
# })
# 
idcomu = sample(V(g)$p,1)
tamcomu = length(V(g)$p[V(g)$p==idcomu])
nvertices = vcount(g)
g = extinction(g, comu = idcomu)

test_that("Função extinction",{
  
  expect_that(vcount(g),equals(nvertices-tamcomu))
  nv = length(V(g)$p[V(g)$p==idcomu])
  expect_that(nv,equals(0))
  expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(g)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(g)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(g), is_more_than(minsize-1))
  expect_that(maiorComunidade(g), is_less_than(maxsize+1))
  
})




#ilustrarGrafo(g)



