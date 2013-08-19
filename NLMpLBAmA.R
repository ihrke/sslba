rm(list=ls())
dat.name <- "sdat.seen.list"
type="nlm"; mind=.1
model.class="pLBA"
model.name="A"
source(paste(model.class,".R",sep=""))
source(paste(model.class,"m",model.name,".R",sep=""))
source("fitter.R")
cores <- 12

start=c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
load("sleep.RData")
datlist <- get(dat.name)
for (i in names(datlist)) 
  attr(datlist[[i]],"p") <- start 

library(snowfall)
sfInit(parallel=TRUE, cpus=cores, type="SOCK")
sfClusterSetupRNG()
sfExportAll()
sfLibrary(truncnorm)
fits <- sfLapply(datlist,dat.fit,fn=objective,type=type,mind=mind)
sfStop()
save(fits,file=paste(model.class,model.name,type,".RData",sep=""))
       
