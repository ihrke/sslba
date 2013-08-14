# install.packages("pso")        # Dont need as fitter contains requried code
# install.packages("snowfall")   # DEoptim doesnt need
# install.packages("rlecuyer")   # Used by snowfall
# install.packages("truncnorm")  # Only required for positive LBA
# install.packages("compiler")   # Maybe built in and dont need
# install.packages("DEoptim")    # 

#########################  LNR model
rm(list=ls())
source("fitter.R")
source("pLBA.R")
source("pLBAmA.R")
start=c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
sv=1; pgf=1e-6; ptf=1e-6

start=astart
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


