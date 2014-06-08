library(igraph)

grafoInicial = criarGrafoInicial()


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

g = born(g)
g = born(g)

test_that("Usando a função born 3 vezes",{
  idcomu = max(V(g)$p)
  tamcomu = length(V(g)$p[V(g)$p==idcomu])
  
  expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(g)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(g)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(g), is_more_than(minsize-1))
  expect_that(maiorComunidade(g), is_less_than(maxsize+1))
})

g = born(g)
g = born(g)


test_that("Usando a função born 5 vezes",{
  idcomu = max(V(g)$p)
  tamcomu = length(V(g)$p[V(g)$p==idcomu])
  
  expect_that(mean(degree(g)), equals(avgdegree, tolerance=1))
  expect_that(max(degree(g)),is_less_than(maxdegree+1))
  mixingReal = calculaMixing(g)
  expect_that(mixingReal, equals(mixing, tolerance = toleranciamixing, scale=1))
  expect_that(menorComunidade(g), is_more_than(minsize-1))
  expect_that(maiorComunidade(g), is_less_than(maxsize+1))
})



# aux = T
# while (aux){
#   g = born(g)
#   m = calculaMixing(g)
#   print(m)
#   if ((m>(mixing+toleranciamixing)) || (m < (mixing-toleranciamixing)) ){
#     aux=F
#   }
# }




