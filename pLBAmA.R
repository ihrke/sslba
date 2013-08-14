# 8 parameter model: c("ster","ter","A","Bs","B","Vs","V","v")
# ster=stop accumulator, ter=go accumlator ter (same for both)
# A=common for all accumulators
# Bs=stop accumlator b-A, B=go accumulator b-A (same for both)
# Vs=stop accumulator rate
# V=matching go accumulator rate (same for both) 
# v=mismatching go accumulator rate (same for both)
#
# Fixed
# sv=1
# gf=1e-6 (probability of go fail)
# tf=1e-6 (probability of stop trigger fail)

#########################  Model A

trans <- function(p) {
  p <- c(log(p[1:5]),p[6:8])
  names(p) <- c("ster","ter","A","Bs","B","Vs","V","v") 
  p  
}

untrans <- function(p) {
  p <- c(exp(p[1:5]),p[6:8])
  names(p) <- c("ster","ter","A","Bs","B","Vs","V","v") 
  p
}

# USED BY PSO and DE
lower=trans(c(ster=1e-5,ter=1e-5,A=1e-5,Bs=1e-5,B=1e-5,Vs=-5,V=-5,v=-5))
upper=trans(c(ster=1,   ter=1,   A=5,   Bs=5,   B=5,   Vs=5, V=5, v=5))

# p is TRANSFORMED version
make.parlists <- function(p,sv=1,gf=1e-6,tf=1e-6) {
  # Accumulators in stop, response 1, response 2 
  p <- untrans(p)     
  list(
      D1S1=list(ter=p["ter"],
        A=rep(p["A"],2),
        sv=rep(sv,2),
        b=p["A"]+rep(p["B"],2),
        v=c(p["V"],p["v"]),
        gf=gf,tf=tf),
      D1S2=list(ter=p["ter"],
        A=rep(p["A"],2),
        sv=rep(sv,2),
        b=p["A"]+rep(p["B"],2),
        v=c(p["v"],p["V"]),
        gf=gf,tf=tf),
      D2S1=list(ter=p["ter"],
        A=rep(p["A"],2),
        sv=rep(sv,2),
        b=p["A"]+rep(p["B"],2),
        v=c(p["V"],p["v"]),
        gf=gf,tf=tf),           
      D2S2=list(ter=p["ter"],
        A=rep(p["A"],2),
        sv=rep(sv,2),
        b=p["A"]+rep(p["B"],2),
        v=c(p["v"],p["V"]),
        gf=gf,tf=tf),
      D1S1s=list(ter=p[c("ster","ter")],
        A=rep(p["A"],3),
        sv=rep(sv,3),
        b=p["A"]+rep(p["B"],3),
        v=p[c("Vs","V","v")],
        gf=gf,tf=tf),
      D1S2s=list(ter=p[c("ster","ter")],
        A=rep(p["A"],3),
        sv=rep(sv,3),
        b=p["A"]+rep(p["B"],3),
        v=p[c("Vs","v","V")],
        gf=gf,tf=tf),
      D2S1s=list(ter=p[c("ster","ter")],
        A=rep(p["A"],3),
        sv=rep(sv,3),
        b=p["A"]+rep(p["B"],3),
        v=p[c("Vs","V","v")],
        gf=gf,tf=tf),           
      D2S2s=list(ter=p[c("ster","ter")],
        A=rep(p["A"],3),
        sv=rep(sv,3),
        b=p["A"]+rep(p["B"],3),
        v=p[c("Vs","v","V")],
        gf=gf,tf=tf)
  )
}


