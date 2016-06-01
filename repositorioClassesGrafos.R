
rm(list = ls())
pasta = "C:\\Users\\sandy\\Dropbox\\2014_Sandy\\DynamicCommunities\\DynamicCommunities\\"
setwd(pasta)
library(igraph)
library(testthat)
library(stringr)

nclasses = 4
nparametros = 6
classesGrafos = matrix(rep(0,nclasses*nparametros),nrow=nclasses)
colnames(classesGrafos) = c("nv","maxs","mins","maxd", "avgd", "mix" )
classesGrafos[1,] = c(300,150,30,60,30,5)
classesGrafos[2,] = c(300,150,30,60,30,20)
classesGrafos[3,] = c(600,300,60,120,60,5)
classesGrafos[4,] = c(600,300,60,120,60,20)


pathParametros <- function(classe){
  path = paste(pasta,"ClassesGrafos\\",sep="")
  aux = ""
  for (i in 1:nparametros){
    aux = paste(aux,colnames(classesGrafos)[i],classesGrafos[classe,i],sep="")
  }
  path = paste(path,aux,"\\",sep="")
  return(path)
}


msgDebug = T
testarErro = T
testeCompleto = F
listaGrafos = F
arquivoErro = paste(pasta,"Erros.dat",sep="")
pastaGrafos = paste(pasta, "ListaGrafos\\",sep="")
ngraph = 1

# if(!testarErro){
#   seed = sample(1:1000,1)
#   set.seed(seed)
# }else{
#   seed = read.table(arquivoErro,sep="\t")
#   seed = seed[nrow(seed),1]
#   set.seed(seed)
# }

nseeds = 100
seeds = sample(1:1e+05,nseeds,replace=F)
tamanhoErro = read.table(arquivoErro,sep="\t")
tamanhoErro = nrow(tamanhoErro)
probfuncao = c(0,0,1,0)
probfuncao = rep(1,4)/4

for(seed in seeds){
  source("DynamicalCommunities.R")
  set.seed(seed)
  cat("Seed:",seed,"\n")
  classe = sample(c(1,2,3,4),1)
  path = pathParametros(classe)
  
  nvertices = as.numeric(classesGrafos[classe,"nv"])
  avgdegree = as.numeric(classesGrafos[classe,"avgd"])
  maxdegree = as.numeric(classesGrafos[classe,"maxd"])
  mixing = as.numeric(classesGrafos[classe,"mix"]/100)
  if (classe %% 2 !=0){
    toleranciamixing = 0.01
  }else{
    toleranciamixing = 0.03
  }
  toleranciagrau = 2
  minsize = as.numeric(classesGrafos[classe,"mins"])
  maxsize = as.numeric(classesGrafos[classe,"maxs"])
  
  cat("Classe:",classe,"\n")
#   cat("Path:",path,"\n")
#   cat("nv:",nvertices,"\n")
#   cat("maxs:",maxsize,"\n")
#   cat("mins:",minsize,"\n")
#   cat("maxd:",maxdegree,"\n")
#   cat("avgd:",avgdegree,"\n")
#   cat("mix:",mixing,"\n")
  resultadosTestes = test_file("testDynamicalCommunities.R",reporter = "summary")
  
  ntestes = length(resultadosTestes)
  for(i in 1:ntestes){
    nresults = length(resultadosTestes[[i]]$results)
    
    for (j in 1:nresults){
      txt = as.character(resultadosTestes[[i]]$results[[j]])
      if (grepl("failure",txt)){
        st = resultadosTestes[[i]]$results[[j]]$message
        st = str_replace_all(st,"\n",". ")
        aux = c(seed,st)
        write(aux,arquivoErro,ncolumns=2,append=T,sep="\t")
      }
    }
  }
}

novoTamanhoErro = read.table(arquivoErro,sep="\t")
novoTamanhoErro = nrow(novoTamanhoErro)

if(novoTamanhoErro > tamanhoErro){
  cat("Erros encontrados:",novoTamanhoErro-tamanhoErro,"\n")
}


