
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
testarErro = F
testeCompleto = F
arquivoErro = paste(pasta,"Erros.dat",sep="")
if(!testarErro){
  seed = sample(1:1000,1)
  set.seed(seed)
}else{
  seed = read.table(arquivoErro,sep="\t")
  seed = seed[nrow(seed),1]
  set.seed(seed)
}

classe = sample(c(1,3,2,4),1)
path = pathParametros(classe)

nvertices = as.numeric(classesGrafos[classe,"nv"])
avgdegree = as.numeric(classesGrafos[classe,"avgd"])
maxdegree = as.numeric(classesGrafos[classe,"maxd"])
mixing = as.numeric(classesGrafos[classe,"mix"]/100)
toleranciamixing = 0.03
minsize = as.numeric(classesGrafos[classe,"mins"])
maxsize = as.numeric(classesGrafos[classe,"maxs"])

cat("Classe:",classe,"\n")
cat("Path:",path,"\n")
cat("nv:",nvertices,"\n")
cat("maxs:",maxsize,"\n")
cat("mins:",minsize,"\n")
cat("maxd:",maxdegree,"\n")
cat("avgd:",avgdegree,"\n")
cat("mix:",mixing,"\n")
source("DynamicalCommunities.R")


