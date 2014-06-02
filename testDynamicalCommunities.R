library(igraph)




test_that("Verificando o grafo retornado pela função criarGrafoInicial",{
  grafoInicial = criarGrafoInicial()
  
  expect_that(is.igraph(grafoInicial), is_true())
  
  expect_that(vcount(grafoInicial),equals(nvertices))
  expect_that(mean(degree(grafoInicial)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(grafoInicial)),equals(maxdegree))
  mixingReal = calculaMixing(grafoInicial)
  expect_that(mixingReal, equals(mixing, tolerance = 0.05))
  expect_that(menorComunidade(grafoInicial), is_more_than(minsize-1))
  expect_that(maiorComunidade(grafoInicial), is_less_than(maxsize+1))
  
})

