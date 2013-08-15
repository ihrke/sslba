# install.packages("pso")        # Dont need as fitter contains requried code
# install.packages("snowfall")   # DEoptim doesnt need
# install.packages("rlecuyer")   # Used by snowfall
# install.packages("truncnorm")  # Only required for positive LBA
# install.packages("compiler")   # Maybe built in and dont need
# install.packages("DEoptim")    # 

######################  READ IN REALL DATA
load("sdat.RData") # produced by Andrew's parser
names(sdat)[3] <- "D"
sdat$SSD[!is.finite(sdat$SSD)] <- Inf
levels(sdat$D) <- c("normal","deprived")
sdat.seen <- sdat[sdat$SEEN==1,-11]
sdat.all <- sdat[,-11]
snams <- levels(sdat$s)
sdat.all.list <- vector(mode="list",length=length(snams))
names(sdat.all.list) <- snams
sdat.seen.list <- sdat.all.list 
for (i in snams) {
  sdat.seen.list[[i]] <- sdat.seen[sdat.seen$s==i,c("D","S","SSD","R","C","RT")] 
  sdat.all.list [[i]] <- sdat.all[sdat.all$s==i,c("D","S","SSD","R","C","RT")]
}
save(sdat.seen,sdat.all,sdat.seen.list,sdat.all.list,file="sleep.RData")

#########################  LNR model
rm(list=ls())
source("fitter.R")
source("pLBA.R")
source("pLBAmA.R")
start=c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
sv=1; pgf=1e-6; ptf=1e-6

start=start
untrans(trans(start))
untrans(lower)
untrans(upper)
make.parlists(trans(start))

adat <- sim.a(1e4,astart)
odat <- sim.o(1e4,ostart)

dat <- adat; start <- astart; objective <- objective.and
dat <- odat; start <- ostart; objective <- objective.or

round(100*tapply(dat$C,list(dat$S),mean))
round(tapply(dat$RT[dat$C],list(dat$S[dat$C]),mean),2)
round(tapply(dat$RT[!dat$C],list(dat$S[!dat$C]),mean),2)
round(tapply(dat$RT[dat$C],list(dat$S[dat$C]),min),2)
round(tapply(dat$RT[dat$C],list(dat$S[dat$C]),max),2)
objective(trans(start),dat)
fit <- fit.one(trans(start),dat,type="nlm",obj=objective)


