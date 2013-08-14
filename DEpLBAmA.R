rm(list=ls())
type="de"; mind=0.1
model.class="pLBA"
model.name="A"
source(paste(model.class,".R",sep=""))
source(paste(model.class,"m",model.name,".R",sep=""))
source("fitter.R")
cores <- 12

ostart=c(ter=.2,A=.2,By=.8,Bn=.3,V=2,v=0)
astart=c(ter=.2,A=.2,By=.3,Bn=.8,V=2,v=0)
load("datlist.RData")
for (i in names(datlist)) 
  if (attr(datlist[[i]],"AO")=="and") 
    attr(datlist[[i]],"p") <- trans(astart) else
      attr(datlist[[i]],"p") <- trans(ostart)

# load("datlist.RData")
# load("LBApUCIPAnlm.RData"); nlm<- fits; rm(fits) 
# load("LBApUCIPApso.RData"); pso <- fits; rm(fits) 
# best <- apply(rbind(unlist(lapply(nlm,function(x){x$value})),
#                     unlist(lapply(pso,function(x){x$value}))),2,which.min)
# for (i in names(datlist)) {
#   if (best[[i]]==1) start <- nlm[[i]]$par else start <- pso[[i]]$par
#   attr(datlist[[i]],"p") <- start
# }

library(snowfall)
sfInit(parallel=TRUE, cpus=cores, type="SOCK")
sfClusterSetupRNG()
sfExportAll()
sfLibrary(truncnorm)
fits <- sfLapply(datlist,dat.fit,type=type,mind=mind,psotype=psotype)
sfStop()
save(fits,file=paste(model.class,model.name,type,".RData",sep=""))
       
