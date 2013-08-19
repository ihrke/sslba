# TESTS OF DENSITY FUNCITONS

source("ssLBA.R")

ster=.1;ter=.2;A=.2;Bs=.5;B=.8;Vs=2;V=1;v=0; sv=1; pgf=0; ptf=0

####################### GO
par(mfrow=c(1,2))
x=seq(.21,2,.01); upper=Inf

pl=list(ter=ter,A=rep(A,2),sv=rep(sv,2),b=c(A+B,A+B),v=c(V,v)) # Go left
plot(x,dfun.1(x,pl),type="l",ylab="Density",main="Go left stimulus")
legend("topright",c("correct","error"),lty=c(1,3))
lines(x,dfun.2(x,pl),lty=3)
pLGL=integrate(dfun.1,pl$ter,upper,pl=pl)$value
pRGL=integrate(dfun.2,pl$ter,upper,pl=pl)$value
cat("Probability go left stiim (left and right response and sum\n")
print(c(pLGL,pRGL,pLGL+pRGL))

pl=list(ter=ter,A=rep(A,2),sv=rep(sv,2),b=c(A+B,A+B),v=c(v,V)) # Go right
plot(x,dfun.2(x,pl),type="l",ylab="Density",main="Go right stimulus")
legend("topright",c("correct","error"),lty=c(1,3))
lines(x,dfun.1(x,pl),lty=3)
pLGL=integrate(dfun.1,pl$ter,upper,pl=pl)$value
pRGL=integrate(dfun.2,pl$ter,upper,pl=pl)$value
cat("Probability go right stim (left and right response and sum\n")
print(c(pLGL,pRGL,pLGL+pRGL))

####################### STOP
par(mfrow=c(1,2))
x=seq(.21,2,.01); upper=Inf
SSD=.25

# Stop left stimulus
pl=list(ter=c(ster,ter),A=rep(A,3),sv=rep(sv,3),b=c(A+Bs,A+B,A+B),v=c(Vs,V,v)) 
plot(x,dfun.1s(x,pl,SSD),type="l",ylab="Density",main="Stop left stimulus")
legend("topright",c("correct","error"),lty=c(1,3))
lines(x,dfun.2s(x,pl,SSD),lty=3)
pLSL=integrate(dfun.1s,pl$ter[2],upper,pl=pl,SSD=SSD)$value
pRSL=integrate(dfun.2s,pl$ter[2],upper,pl=pl,SSD=SSD)$value
pstop=dfun.stop(pl,SSD)
cat("Probability stop left stim (left and right response, (stop, and sum\n")
print(c(pLSL,pRSL,pstop,pLSL+pRSL+pstop))

#  Stop right stimulus
pl=list(ter=c(ster,ter),A=rep(A,3),sv=rep(sv,3),b=c(A+Bs,A+B,A+B),v=c(Vs,v,V)) 
plot(x,dfun.2s(x,pl,SSD),type="l",ylab="Density",main="Stop right stimulus")
legend("topright",c("correct","error"),lty=c(1,3))
lines(x,dfun.1s(x,pl,SSD),lty=3)
pLSL=integrate(dfun.1s,pl$ter[2],upper,pl=pl,SSD=SSD)$value
pRSL=integrate(dfun.2s,pl$ter[2],upper,pl=pl,SSD=SSD)$value
pstop=dfun.stop(pl,SSD)
cat("Probability stop left stimulus (left and right response, (stop, and sum\n")
print(c(pLSL,pRSL,pstop,pLSL+pRSL+pstop))

# TESTS COMPARING RANDOM AND DENSITY FUNCTIONS

p <- c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
SSD=.25 # Graphs only handle one SSD

################################ ??? WHEN I SET N=LARGE (I.E.,1E5 SEEMS TO FAIL
################################ LOOKS LIKE A BUG IN THE R density FUNCTION)
################################ CAN SOMEONE CHECK THIS OUT? WORKS FINE FOR n=1E4

# set n big to get a good match
sim <- sim.stop(1e4,p,ssds=SSD)

# simulated probability of correct response and of non-response
pC <- tapply(sim$C,list(D=sim$D,S=sim$S,GO=!is.finite(sim$SSD)),mean,na.rm=TRUE)
pNA <- tapply(is.na(sim$R),list(D=sim$D,S=sim$S,GO=!is.finite(sim$SSD)),mean)
pC <- pC*(1-pNA)
pCC <- 1-pNA-pC
pl <- make.parlists(trans(p))
x=seq(.21,2,.01); upper=Inf

# GO: RED=simulated, BLACK=exact, solid=CORRECT, dashed=error
quartz()
par(mfrow=c(2,2))
plot(x,(1-pl$D1S1$gf)*dfun.1(x,pl$D1S1),type="l",ylab="Density",main="GO:D1S1")
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="deprived" & sim$S=="left" & sim$R=="left","RT"])
dns$y <- dns$y*pC["deprived","left","TRUE"]; lines(dns,col="red")
lines(x,(1-pl$D1S1$gf)*dfun.2(x,pl$D1S1),lty=2)
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="deprived" & sim$S=="left" & sim$R=="right","RT"])
dns$y <- dns$y*(pCC["deprived","left","TRUE"]); lines(dns,col="red")
plot(x,(1-pl$D1S2$gf)*dfun.2(x,pl$D1S2),type="l",ylab="Density",main="GO:D1S2")
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="deprived" & sim$S=="right" & sim$R=="right","RT"])
dns$y <- dns$y*pC["deprived","right","TRUE"]; lines(dns,col="red")
lines(x,(1-pl$D1S2$gf)*dfun.1(x,pl$D1S2),lty=2)
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="deprived" & sim$S=="right" & sim$R=="left","RT"])
dns$y <- dns$y*(pCC["deprived","right","TRUE"]); lines(dns,col="red")
plot(x,d(1-pl$D2S1$gf)*fun.1(x,pl$D2S1),type="l",ylab="Density",main="GO:D2S1")
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="normal" & sim$S=="left" & sim$R=="left","RT"])
dns$y <- dns$y*pC["normal","left","TRUE"]; lines(dns,col="red")
lines(x,(1-pl$D2S1$gf)*dfun.2(x,pl$D2S1),lty=2)
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="normal" & sim$S=="left" & sim$R=="right","RT"])
dns$y <- dns$y*(pCC["normal","left","TRUE"]); lines(dns,col="red")
plot(x,(1-pl$D2S2$gf)*dfun.2(x,pl$D2S2),type="l",ylab="Density",main="GO:D2S2")
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="deprived" & sim$S=="right" & sim$R=="right","RT"])
dns$y <- dns$y*pC["normal","right","TRUE"]; lines(dns,col="red")
lines(x,(1-pl$D2S2$gf)*dfun.1(x,pl$D2S2),lty=2)
dns <- density(sim[!is.finite(sim$SSD) & sim$D=="deprived" & sim$S=="right" & sim$R=="left","RT"])
dns$y <- dns$y*(pCC["normal","right","TRUE"]); lines(dns,col="red")

# STOP: RED=simulated, BLACK=exact, solid=CORRECT, dashed=error
quartz()
par(mfrow=c(2,2))
plot(x,(1-pl$D1S1s$gf)*(pl$D1S1s$tf + (1-pl$D1S1s$tf)*dfun.1s(x,pl$D1S1s,SSD)),type="l",ylab="Density",main="STOP:D1S1")
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="deprived" & tmp$S=="left" & tmp$R=="left","RT"])
dns$y <- dns$y*pC["deprived","left","FALSE"]; lines(dns,col="red")
lines(x,(1-pl$D1S1s$gf)*(pl$D1S1s$tf + (1-pl$D1S1s$tf)*dfun.2s(x,pl$D1S1s,SSD)),lty=2)
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="deprived" & tmp$S=="left" & tmp$R=="right","RT"])
dns$y <- dns$y*(pCC["deprived","left","FALSE"]); lines(dns,col="red")
plot(x,(1-pl$D1S2s$gf)*(pl$D1S2s$tf + (1-pl$D1S2s$tf)*dfun.2s(x,pl$D1S2s,SSD)),type="l",ylab="Density",main="STOP:D1S2")
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="deprived" & tmp$S=="right" & tmp$R=="right","RT"])
dns$y <- dns$y*pC["deprived","right","FALSE"]; lines(dns,col="red")
lines(x,(1-pl$D1S2s$gf)*(pl$D1S2s$tf + (1-pl$D1S2s$tf)*dfun.1s(x,pl$D1S2s,SSD)),lty=2)
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="deprived" & tmp$S=="right" & tmp$R=="left","RT"])
dns$y <- dns$y*(pCC["deprived","right","FALSE"]); lines(dns,col="red")
plot(x,(1-pl$D2S1s$gf)*(pl$D2S1s$tf + (1-pl$D2S1s$tf)*dfun.1s(x,pl$D2S1s,SSD)),type="l",ylab="Density",main="STOP:D2S1")
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="normal" & tmp$S=="left" & tmp$R=="left","RT"])
dns$y <- dns$y*pC["normal","left","FALSE"]; lines(dns,col="red")
lines(x,(1-pl$D2S1s$gf)*(pl$D2S1s$tf + (1-pl$D2S1s$tf)*dfun.2s(x,pl$D2S1s,SSD)),lty=2)
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="normal" & tmp$S=="left" & tmp$R=="right","RT"])
dns$y <- dns$y*(pCC["normal","left","FALSE"]); lines(dns,col="red")
plot(x,(1-pl$D2S2s$gf)*(pl$D2S2s$tf + (1-pl$D2S2s$tf)*dfun.2s(x,pl$D2S2s,SSD)),type="l",ylab="Density",main="STOP:D2S2")
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="normal" & tmp$S=="right" & tmp$R=="right","RT"])
dns$y <- dns$y*pC["normal","right","FALSE"]; lines(dns,col="red")
lines(x,(1-pl$D2S2s$gf)*(pl$D2S2s$tf + (1-pl$D2S2s$tf)*dfun.1s(x,pl$D2S2s,SSD)),lty=2)
tmp <- sim[!is.na(sim$R),]; dns <- density(tmp[is.finite(tmp$SSD) & tmp$D=="deprived" & tmp$S=="right" & tmp$R=="left","RT"])
dns$y <- dns$y*(pCC["normal","right","FALSE"]); lines(dns,col="red")

# simulated - predicted go fails
pNA[,,"TRUE"] - c(pl$D1S1$gf,pl$D2S1$gf,pl$D1S2$gf,pl$D2S2$gf)
# simulated - predicted stops
pNA[,,"FALSE"] - c(pl$D1S1s$gf + (1-pl$D1S1s$gf)*(1-pl$D1S1s$tf)*dfun.stop(pl$D1S1s,SSD),
                   pl$D2S1s$gf + (1-pl$D2S1s$gf)*(1-pl$D2S1s$tf)*dfun.stop(pl$D2S1s,SSD),
                   pl$D1S2s$gf + (1-pl$D1S2s$gf)*(1-pl$D1S2s$tf)*dfun.stop(pl$D1S2s,SSD),
                   pl$D2S2s$gf + (1-pl$D2S2s$gf)*(1-pl$D2S2s$tf)*dfun.stop(pl$D2S2s,SSD))
