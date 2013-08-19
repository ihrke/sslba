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

#########################  ssLBA model A
rm(list=ls())
source("fitter.R")
source("pLBA.R")
source("pLBAmA.R")

start=c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
# sv=1; pgf=0; ptf=0

# check you get back what you started with!
untrans(trans(start)) 
# check upper and lower look right
untrans(lower)
untrans(upper)
# Check parlist looks right
make.parlists(trans(start))

# This may bug if n is to small to get stop trials
dat <- sim.stop(n=1e4,p=start)

############## Check simulated data looks right
nsig=3
######## look at the distribution 
round(min(dat$RT,na.rm=T),nsig)
round(max(dat$RT,na.rm=T),nsig)
# you will usually get a few very big values
tail(sort(dat$RT))
# ignore them
upper <- 3
plot(density(dat$RT[!is.na(dat$RT) & dat$RT<upper]),main="alldata")
round(mean(dat$RT>upper,na.rm=TRUE),nsig)

######### Non-response trials: dat$SSD!=Inf picks out stop trials
round(tapply(is.na(dat$R),list(stop=dat$SSD!=Inf,D=dat$D,S=dat$S),mean),nsig)
########## Pick out stop data
sdat=dat[dat$SSD!=Inf,]
######### stop proportions at each SSD
round(tapply(is.na(sdat$R),list(D=sdat$D,SSD=sdat$SSD,S=sdat$S),mean),nsig)
######### Percent correct on non-strop trials
round(tapply(dat$C,list(stop=dat$SSD!=Inf,D=dat$D,S=dat$S),mean,na.rm=TRUE),nsig)
########## RT on trials with responses
round(tapply(dat$RT,list(stop=dat$SSD!=Inf,D=dat$D,S=dat$S),mean,na.rm=TRUE),nsig)
########## Correct trials
cdat <- dat[dat$C,]
round(tapply(cdat$RT,list(stop=cdat$SSD!=Inf,D=cdat$D,S=cdat$S),mean,na.rm=TRUE),nsig)
########### Error trials
edat <- dat[!dat$C,]
round(tapply(edat$RT,list(stop=edat$SSD!=Inf,D=edat$D,S=edat$S),mean,na.rm=TRUE),nsig)

# Does objective evaluate to a finite number?
objective(trans(start),dat)
# what does likelihood look like against data
View(cbind(dat,L=objective(trans(start),dat,lvals=TRUE)))
# Does objective evaluate to a finite number at bounds (needed for pso and de)
objective(lower,dat)
objective(upper,dat)


# PARAMTER RECOVERY from true value
fit <- fit.one(trans(start),dat,fn=objective)
# What is the estimate
untrans(fit$par)
# How good is recovery?
untrans(fit$par)-start

###### Check out performance of different fitters

# Make some small data
dat <- sim.stop(n=1e3,p=start)

# NLM should be fastest and often most accurate
nlm.time <- system.time(fit.nlm <- fit.one(trans(start),dat,fn=objective)) # nlm by default

# Simplex is usually slower (when you multi-fit, 
# if you dont multi-fit it is usually less accurate)
simplex.time.one <- system.time(fit.simplex <- fit.one(trans(start),dat,fn=objective,type="simplex"))[3]
# Try a multi-fit with simplex, uses default of a < .1 decrease to stop
# this defualt may need to be played with in real data!
simplex.time.multi <- system.time(fit.simplex.multi <- fit.multi(trans(start),dat,fn=objective,type="simplex"))[3]


# PSO and DE are very slow, set tracing to see if it is stuck!
# Set back to trace=0 for PSO to not do this (there must be a similar option
# in DE but I didnt set it)
# There are LOTS of options to play with for these algorithms, I havent used
# much more thatn the defaults, for PSO "SPSO2011" can really slow things down
# (and is really slow already) but can be better. There are at least 6 DE 
# algorithms, again I have only tried the default

# These use true values as start point
pso.time.start <- system.time(fit.pso.start <- 
    fit.one(p=trans(start),dat=dat,fn=objective,type="pso",control=list(trace=1)))[3]

################  UP TO HERE

# Note you could run de (but not pso) across multiple cores for a single fit, 
# I have not done so 
de.time.start <- system.time(fit.de.start <- 
    fit.one(trans(start),dat,fn=objective,type="de"))[3]

# An advantage of DE and PSO is they dont need start values
pso.time <- system.time(fit.pso <- 
  fit.one(NA,dat,fn=objective,type="pso",control=list(trace=1)))[3]
de.time <- system.time(fit.de <- fit.one(NA,dat,fn=objective,type="de"))[3]


