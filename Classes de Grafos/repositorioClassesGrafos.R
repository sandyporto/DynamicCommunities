library(igraph)

pasta = "C:/Users/Sandy/Dropbox/2014_Sandy/DynamicCommunities/DynamicCommunities/"

nverticeslist = c(300,600)
avgdegreelist = c()
maxdegreelist = c()
mixinglist = c(0.05,0.2)
minsizelist = c()
maxsizelist = c()

numeropastas = 0

for (nv in nverticeslist){
  maxsizelist = c(round(nv*0.5))
  for (maxs in maxsizelist){
    minsizelist = c(round(nv*0.05))
    for (mins in minsizelist){
      maxdegreelist = c(round(nv*0.2))
      for (maxd in maxdegreelist){
        avgdegreelist = c(round(nv*0.1))
        for (avgd in avgdegreelist){
          for (mix in mixinglist){
            nomepasta = paste("nv",nv,
                              "maxs",maxs,
                              "mins",mins,
                              "maxd",maxd,
                              "avgd",avgd,
                              "mix",(mix*100),sep="")
            print(nomepasta,quote=F)
            numeropastas = numeropastas + 1
          }
        }
      }
    }
  }
}

print(numeropastas)