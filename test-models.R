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
sv=1; pgf=0; ptf=0

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

# Check out performance of different fitters
dat <- sim.stop(n=1e3,p=start)
nlm.time <- system.time(fit.nlm <- fit.one(trans(start),dat,fn=objective)) # nlm by default
simplex.time.one <- system.time(fit.simplex <- fit.one(trans(start),dat,fn=objective,type="simplex"))[3]
# Try a multi-fit with simplex
simplex.time.multi <- system.time(fit.simplex.multi <- fit.multi(trans(start),dat,fn=objective,type="simplex"))[3]


#####################################  UP TO HERE
###############  PROBLEMS WITH UPPER AND LOWER BOUNDS CAUSE ERROR IN DEoptim
# These use true values as start point
pso.time.start <- system.time(fit.simplex <- fit.one(trans(start),dat,fn=objective,type="pso"))[3]
de.time.start <- system.time(fit.simplex <- fit.one(trans(start),dat,fn=objective,type="de"))[3]
# These dont use start values
pso.time.nostart <- system.time(fit.simplex <- fit.one(NA,dat,fn=objective,type="pso"))[3]
de.time.nostart <- system.time(fit.simplex <- fit.one(NA,dat,fn=objective,type="de"))[3]


