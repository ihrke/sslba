# LBA (any sort) density functions for go and stop archiectures

##################################### 2 accumulator go architecture
# respond 1st accumulator
dfun.1 <- function(t,pl) {
  t <- pmax(t-pl$ter[1],1e-5) # doesnt allow a different ter for each accumulator
  out <- p1(t,A=pl$A[1],b=pl$b[1],v=pl$v[1],sv=pl$sv[1])*
    (1-c1(t,A=pl$A[2],b=pl$b[2],v=pl$v[2],sv=pl$sv[2]))
  out[!is.finite(out) | out<0] <- 0
  out
}

# respond 2nd accumulator
dfun.2 <- function(t,pl) {
  t <- pmax(t-pl$ter[1],1e-5) # doesnt allow a different ter for each accumulator
  out <- p1(t,A=pl$A[2],b=pl$b[2],v=pl$v[2],sv=pl$sv[2])*
    (1-c1(t,A=pl$A[1],b=pl$b[1],v=pl$v[1],sv=pl$sv[1]))
  out[!is.finite(out)  | out<0] <- 0
  out
}

##################################### 3 accumulator stop architecture

# Respond 2nd accumulator (1st response) in stop architecture
dfun.1s <- function(t,pl,SSD) {
  tstop <- pmax(t-pl$ter[1]-SSD,1e-5) 
  t <- pmax(t-pl$ter[2],1e-5) # doesnt allow a different ter for each accumulator
  out <- p1(t,A=pl$A[2],b=pl$b[2],v=pl$v[2],sv=pl$sv[2])*
    (1-c1(tstop,A=pl$A[1],b=pl$b[1],v=pl$v[1],sv=pl$sv[1]))*
    (1-c1(t,A=pl$A[3],b=pl$b[3],v=pl$v[3],sv=pl$sv[3]))
  out[!is.finite(out)  | out<0] <- 0
  out
}

# Respond 3rd accumulator (2nd response) in stop architecture
dfun.2s <- function(t,pl,SSD) {
  tstop <- pmax(t - pl$ter[1] - SSD,1e-5) 
  t <- pmax(t - pl$ter[2],1e-5) # doesnt allow a different ter for each accumulator
  out <- p1(t,A=pl$A[3],b=pl$b[3],v=pl$v[3],sv=pl$sv[3])*
    (1-c1(tstop,A=pl$A[1],b=pl$b[1],v=pl$v[1],sv=pl$sv[1]))*
    (1-c1(t,A=pl$A[2],b=pl$b[2],v=pl$v[2],sv=pl$sv[2]))
  out[!is.finite(out)  | out<0] <- 0
  out
}

# Sucessful stop 
dfun.stop <- function(pl,SSD) {

  tmpf=function(t,pl,SSD) p1(z=pmax(t - pl$ter[1] - SSD,1e-5),A=pl$A[1],b=pl$b[1],v=pl$v[1],sv=pl$sv[1])* 
                          (1-c1(pmax(t - pl$ter[2],1e-5),A=pl$A[2],b=pl$b[2],v=pl$v[2],sv=pl$sv[2]))*
                          (1-c1(pmax(t - pl$ter[3],1e-5),A=pl$A[3],b=pl$b[3],v=pl$v[3],sv=pl$sv[3])) 
  pstop=integrate(f=tmpf,lower=pl$ter[1]+SSD,upper=Inf,pl=pl,SSD=SSD)$value
  pmax(pmin(pstop,1),0) # protect from numerical errors
}


# TESTS OF DENSITY FUNCITONS

# ter=.2;A=.2;sv=1;By=.8;Bn=.3;V=2;v=0
# 
# par(mfrow=c(2,2))
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(v,v,V,V)) # NT and
# x=seq(.21,1.5,.01); upper=Inf
# plot(x,dfun.or.no(x,pl),type="l",ylab="Density",main="NT or")
# lines(x,dfun.or.yes(x,pl),lty=3)
# pyNT=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnNT=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(V,v,v,V)) # T1 and
# plot(x,dfun.or.yes(x,pl),type="l",ylab="Density",main="T1 or")
# lines(x,dfun.or.no(x,pl),lty=3)
# pyT1=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnT1=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(v,V,V,v)) # T1 and
# plot(x,dfun.or.yes(x,pl),type="l",ylab="Density",main="T2 or")
# lines(x,dfun.or.no(x,pl),lty=3)
# pyT2=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnT2=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(V,V,v,v)) # DT and
# plot(x,dfun.or.yes(x,pl),type="l",ylab="Density",main="DT or")
# lines(x,dfun.or.no(x,pl),lty=3)
# pyDT=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnDT=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# c(pyNT,pnNT,pyNT+pnNT)
# c(pyT1,pnT1,pyT1+pnT1)
# c(pyT2,pnT2,pyT2+pnT2)
# c(pyDT,pnDT,pyDT+pnDT)
# 
# 
# ter=.2;A=.2;sv=1;By=.3;Bn=.8;V=2;v=0
# 
# par(mfrow=c(2,2))
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(v,v,V,V)) # NT and
# x=seq(.21,1.5,.01); upper=Inf
# plot(x,dfun.and.no(x,pl),type="l",ylab="Density",main="NT and")
# lines(x,dfun.and.yes(x,pl),lty=3)
# pyNT=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnNT=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(V,v,v,V)) # T1 and
# plot(x,dfun.and.no(x,pl),type="l",ylab="Density",main="T1 and")
# lines(x,dfun.and.yes(x,pl),lty=3)
# pyT1=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnT1=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(v,V,V,v)) # T1 and
# plot(x,dfun.and.no(x,pl),type="l",ylab="Density",main="T2 and")
# lines(x,dfun.and.yes(x,pl),lty=3)
# pyT2=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnT2=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# pl=list(ter=ter,A=rep(A,4),sv=rep(sv,4),b=c(A+By,A+By,A+Bn,A+Bn),v=c(V,V,v,v)) # DT and
# plot(x,dfun.and.yes(x,pl),type="l",ylab="Density",main="DT and")
# lines(x,dfun.and.no(x,pl),lty=3)
# pyDT=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnDT=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# c(pyNT,pnNT,pyNT+pnNT)
# c(pyT1,pnT1,pyT1+pnT1)
# c(pyT2,pnT2,pyT2+pnT2)
# c(pyDT,pnDT,pyDT+pnDT)
# 


# TESTS COMPARING RANDOM AND DESNITY FUNCTIONS
# p=c(ter=.2,A=.2,oBy=.8,oBn=.3,aBy=.3,aBn=.8,V=2,v=0); sv=1
# sim <- sim.ao(1e4,p)
# pC <- tapply(sim$C,list(sim$S,sim$AO),mean)
# par(mfrow=c(2,2))
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["oBy"],p["A"]+p["oBy"],p["A"]+p["oBn"],p["A"]+p["oBn"]),
#         v=c(p["v"],p["v"],p["V"],p["V"])) # NT and
# x=seq(.21,1.5,.01); upper=Inf
# pyNT=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnNT=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.or.no(x,pl),type="l",ylab="Density",main="NT or")
# lines(x,dfun.or.yes(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="or" & sim$R=="no" & sim$S=="NT","RT"]); tmp$y <- tmp$y*pC["NT","or"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="or" & sim$R=="yes" & sim$S=="NT","RT"]); tmp$y <- tmp$y*(1-pC["NT","or"]); lines(tmp,lty=3,col="red")
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["oBy"],p["A"]+p["oBy"],p["A"]+p["oBn"],p["A"]+p["oBn"]),
#         v=c(p["V"],p["v"],p["v"],p["V"])) # NT and
# pyT1=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnT1=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.or.yes(x,pl),type="l",ylab="Density",main="T1 or")
# lines(x,dfun.or.no(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="or" & sim$R=="yes" & sim$S=="T1","RT"]); tmp$y <- tmp$y*pC["T1","or"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="or" & sim$R=="no" & sim$S=="T1","RT"]); tmp$y <- tmp$y*(1-pC["T1","or"]); lines(tmp,lty=3,col="red")
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["oBy"],p["A"]+p["oBy"],p["A"]+p["oBn"],p["A"]+p["oBn"]),
#         v=c(p["v"],p["V"],p["V"],p["v"])) # NT and
# pyT2=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnT2=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.or.yes(x,pl),type="l",ylab="Density",main="T2 or")
# lines(x,dfun.or.no(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="or" & sim$R=="yes" & sim$S=="T2","RT"]); tmp$y <- tmp$y*pC["T2","or"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="or" & sim$R=="no" & sim$S=="T2","RT"]); tmp$y <- tmp$y*(1-pC["T2","or"]); lines(tmp,lty=3,col="red")
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["oBy"],p["A"]+p["oBy"],p["A"]+p["oBn"],p["A"]+p["oBn"]),
#         v=c(p["V"],p["V"],p["v"],p["v"])) # NT and
# pyDT=integrate(dfun.or.yes,pl$ter,upper,pl=pl)$value
# pnDT=integrate(dfun.or.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.or.yes(x,pl),type="l",ylab="Density",main="DT or")
# lines(x,dfun.or.no(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="or" & sim$R=="yes" & sim$S=="DT","RT"]); tmp$y <- tmp$y*pC["DT","or"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="or" & sim$R=="no" & sim$S=="DT","RT"]); tmp$y <- tmp$y*(1-pC["DT","or"]); lines(tmp,lty=3,col="red")
# c(pyNT,pnNT,pyNT+pnNT)
# c(pyT1,pnT1,pyT1+pnT1)
# c(pyT2,pnT2,pyT2+pnT2)
# c(pyDT,pnDT,pyDT+pnDT)
# 
# 
# par(mfrow=c(2,2))
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["aBy"],p["A"]+p["aBy"],p["A"]+p["aBn"],p["A"]+p["aBn"]),
#         v=c(p["v"],p["v"],p["V"],p["V"])) # NT and
# x=seq(.21,1.5,.01); upper=Inf
# pyNT=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnNT=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.and.no(x,pl),type="l",ylab="Density",main="NT and")
# lines(x,dfun.and.yes(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="and" & sim$R=="no" & sim$S=="NT","RT"]); tmp$y <- tmp$y*pC["NT","and"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="and" & sim$R=="yes" & sim$S=="NT","RT"]); tmp$y <- tmp$y*(1-pC["NT","and"]); lines(tmp,lty=3,col="red")
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["aBy"],p["A"]+p["aBy"],p["A"]+p["aBn"],p["A"]+p["aBn"]),
#         v=c(p["V"],p["v"],p["v"],p["V"])) # NT and
# pyT1=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnT1=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.and.no(x,pl),type="l",ylab="Density",main="T1 and")
# lines(x,dfun.and.yes(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="and" & sim$R=="no" & sim$S=="T1","RT"]); tmp$y <- tmp$y*pC["T1","and"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="and" & sim$R=="yes" & sim$S=="T1","RT"]); tmp$y <- tmp$y*(1-pC["T1","and"]); lines(tmp,lty=3,col="red")
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["aBy"],p["A"]+p["aBy"],p["A"]+p["aBn"],p["A"]+p["aBn"]),
#         v=c(p["v"],p["V"],p["V"],p["v"])) # NT and
# pyT2=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnT2=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.and.no(x,pl),type="l",ylab="Density",main="T2 and")
# lines(x,dfun.and.yes(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="and" & sim$R=="no" & sim$S=="T2","RT"]); tmp$y <- tmp$y*pC["T2","and"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="and" & sim$R=="yes" & sim$S=="T2","RT"]); tmp$y <- tmp$y*(1-pC["T2","and"]); lines(tmp,lty=3,col="red")
# pl=list(ter=p["ter"],A=rep(p["A"],4),sv=rep(sv,4),b=c(p["A"]+p["aBy"],p["A"]+p["aBy"],p["A"]+p["aBn"],p["A"]+p["aBn"]),
#         v=c(p["V"],p["V"],p["v"],p["v"])) # NT and
# pyDT=integrate(dfun.and.yes,pl$ter,upper,pl=pl)$value
# pnDT=integrate(dfun.and.no,pl$ter,upper,pl=pl)$value
# plot(x,dfun.and.yes(x,pl),type="l",ylab="Density",main="DT and")
# lines(x,dfun.and.no(x,pl),lty=3)
# tmp <- density(sim[sim$AO=="and" & sim$R=="yes" & sim$S=="DT","RT"]); tmp$y <- tmp$y*pC["DT","and"]; lines(tmp,col="red")
# tmp <- density(sim[sim$AO=="and" & sim$R=="no" & sim$S=="DT","RT"]); tmp$y <- tmp$y*(1-pC["DT","and"]); lines(tmp,lty=3,col="red")
# c(pyNT,pnNT,pyNT+pnNT)
# c(pyT1,pnT1,pyT1+pnT1)
# c(pyT2,pnT2,pyT2+pnT2)
# c(pyDT,pnDT,pyDT+pnDT)

