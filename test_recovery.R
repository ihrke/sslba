# test parameter recovery with the pLBA model A
source("pLBA.R") 
source("objsim.R") 
pars=c(ster=.1,ter=.2,A=.2,Bs=.5,B=.8,Vs=2,V=1,v=0)
start=c(ster=.15,ter=.25,A=.15,Bs=.6,B=.5,Vs=2.5,V=1.3,v=0.4)

sv=1; pgf=1e-6; ptf=1e-6
ntrials=200

dat=sim.stop(ntrials, pars, stopp=1/3.)


# USED BY PSO and DE
lower=trans(c(ster=1e-5,ter=1e-5,A=1e-5,Bs=1e-5,B=1e-5,Vs=-5,V=-5,v=-5))
upper=trans(c(ster=1,   ter=1,   A=5,   Bs=5,   B=5,   Vs=5, V=5, v=5))

fit=fit.one(trans(start), dat, objective, 'simplex',control=list(maxit=500, trace=10))