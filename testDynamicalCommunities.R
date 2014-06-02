library(igraph)


test_that("",{
  grafoInicial = criarGrafoInicial()
  
  expect_that(is.igraph(grafoInicial), is_true())
  
  
  
})