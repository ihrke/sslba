# OBJECTIVE FUNCITON with GO and TRIGGER failure

# pgf  : probability of a go failure
# ptf  : probability of a trigger failure
# pstop: probability of a sucessful stop given no go or trigger fail 
# Lg(t): likelihood of response at time t on go trial given no go fail
# Ls(t): likelihood of response at time t on stop trial given no go or trigger fail
# 
# GO(NA)  : pgf
# STOP(NA): pgf + (1-pgf)*(1-ptf)*pstop
# GO(t)   : (1-pgf)*Lg(t)
# STOP(t) : (1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]

# p=trans(start)
# p=pl.bad
objective <- function(p,dat,lvals=FALSE) {
  
  L <- numeric(dim(dat)[1])
  pls <- make.parlists(p)
  no.response <- is.na(dat$R)
  stop.trial <- is.finite(dat$SSD)
  for (d in levels(dat$D)) for (s in levels(dat$S)) {

    # Go fail probability 
    pgf <- switch(d,
      "normal"=  switch(s,"left"=pls$D1S1s$gf,"right"=pls$D1S2s$gf),
      "deprived"=switch(s,"left"=pls$D2S1s$gf,"right"=pls$D2S2s$gf))
    
    # Trigger fail probability
    ptf <- switch(d,
      "normal"  =switch(s,"left"=pls$D1S1s$tf,"right"=pls$D1S2s$tf),
      "deprived"=switch(s,"left"=pls$D2S1s$tf,"right"=pls$D2S2s$tf))
    
    # Failed go: GO(NA) = pgf
    is.in <- !stop.trial & no.response & dat$D==d & dat$S==s
    L[is.in] <-rep(pgf,sum(is.in))    
    
    # Sucessful stops
    nstop <- table(dat[stop.trial & no.response & dat$D==d & dat$S==s,"SSD"])
    for (j in as.numeric(names(nstop))) { 
      is.in <- stop.trial & no.response & dat$D==d & dat$S==s & dat$SSD == j
      # STOP(NA): pgf + (1-pgf)*(1-ptf)*pstop
      L[is.in] <- pgf + (1-pgf)*(1-ptf)*rep(switch(d,
        "normal"  =switch(s,"left"=dfun.stop(pls$D1S1s,j),"right"=dfun.stop(pls$D1S2s,j)),
        "deprived"=switch(s,"left"=dfun.stop(pls$D2S1s,j),"right"=dfun.stop(pls$D2S2s,j))),
        nstop[as.character(j)])    
    }
    
    for (r in levels(dat$R)) { # GO(t) and STOP(t)
      
      # Go reponses
      is.in <- !stop.trial & !no.response & dat$D==d & dat$S==s & dat$R==r
      t <- dat[is.in,"RT"]
      if (length(t)>0) {
        # GO(t)=(1-pgf)*Lg(t)
        L[is.in] <- (1-pgf)*switch(d,
          "normal"=switch(s,
            "left" =switch(r,"left"=dfun.1(t,pls$D1S1),"right"=dfun.2(t,pls$D1S1)),
            "right"=switch(r,"left"=dfun.1(t,pls$D1S2),"right"=dfun.2(t,pls$D1S2))),
          "deprived"=switch(s,
            "left" =switch(r,"left"=dfun.1(t,pls$D2S1),"right"=dfun.2(t,pls$D2S1)),
            "right"=switch(r,"left"=dfun.1(t,pls$D2S2),"right"=dfun.2(t,pls$D2S2)))
        )
      }
      
      # STOP trials with a response
      is.in <- stop.trial & !no.response & dat$D==d & dat$S==s & dat$R==r
      t <- dat[is.in,"RT"]
      SSD <- dat[is.in,"SSD"]
      if (length(t)>0) {
        # STOP(t)=(1-pgf)*[ptf*Lg(t) + (1-ptf)*Ls(t)]
        L[is.in] <- (1-pgf)*(ptf*switch(d,
            "normal"=switch(s,   
              "left" =switch(r,"left"=dfun.1(t,pls$D1S1),"right"=dfun.2(t,pls$D1S1)),
              "right"=switch(r,"left"=dfun.1(t,pls$D1S2),"right"=dfun.2(t,pls$D1S2))),
            "deprived"=switch(s, 
              "left" =switch(r,"left"=dfun.1(t,pls$D2S1),"right"=dfun.2(t,pls$D2S1)),
              "right"=switch(r,"left"=dfun.1(t,pls$D2S2),"right"=dfun.2(t,pls$D2S2))))
                     + (1-ptf)*switch(d,
            "normal"=switch(s,
              "left" =switch(r,"left"=dfun.1s(t,pls$D1S1s,SSD),"right"=dfun.2s(t,pls$D1S1s,SSD)),
              "right"=switch(r,"left"=dfun.1s(t,pls$D1S2s,SSD),"right"=dfun.2s(t,pls$D1S2s,SSD))),
            "deprived"=switch(s,
              "left" =switch(r,"left"=dfun.1s(t,pls$D2S1s,SSD),"right"=dfun.2s(t,pls$D2S1s,SSD)),
              "right"=switch(r,"left"=dfun.1s(t,pls$D2S2s,SSD),"right"=dfun.2s(t,pls$D2S2s,SSD)))
        ))
      }
    } # end of response loop      
  } # end of stimulus loop
  if (lvals) L else -2*sum(log(pmax(L,1e-10)))
}


# FULL DESIGN SIMULATION FUNCTION

# nsim=1e2;p=start
# stopp=.25;ssds=c(.1,.2,.3,.4,.5); dl=c("normal","deprived");sl=c("left","right")
sim.stop <- function(nsim,p,stopp=.25,ssds=c(.1,.2,.3,.4,.5),
  dl=c("deprived","normal"),sl=c("left","right")) {
  parlists <- make.parlists(trans(p))
  nssd = length(ssds)
  if (length(stopp)==1) stopp <- rep(stopp/nssd,nssd)
  if (nssd !=length(stopp)) stop("stopp and ssds not same length")
  ssdn <- round(sum(stopp)*nsim/nssd)
  nstop <- nssd*ssdn
  ngo <- nsim-(nstop)  
  
  # Go data
  dat.go <- cbind.data.frame(
    D=factor(rep(dl,each=ngo*2)),
    S=factor(rep(sl,each=ngo,times=2)),
    SSD=rep(Inf,ngo*4),
    rbind(rfun(ngo,parlists$D1S1),
          rfun(ngo,parlists$D1S2),
          rfun(ngo,parlists$D2S1),
          rfun(ngo,parlists$D2S2)
    )
  )
  names(dat.go)[4:5] <- c("R","RT")
  dat.go$R <- factor(dat.go$R,levels=1:2,labels=c("left","right"))
  dat.go$C <- dat.go$R==dat.go$S
  SSD <- rep(ssds,each=ssdn)
  
  # Stop data
  dat.stop <- cbind.data.frame(
    D=factor(rep(dl,each=nstop*2)),
    S=factor(rep(sl,each=nstop,times=2)),
    SSD=rep(ssds,each=ssdn,times=4),
    rbind(rfun(nstop,parlists$D1S1s,SSD),
          rfun(nstop,parlists$D1S2s,SSD),
          rfun(nstop,parlists$D2S1s,SSD),
          rfun(nstop,parlists$D2S2s,SSD)
    )
  )
  names(dat.stop)[4:5] <- c("R","RT")
  dat.stop$R <- factor(dat.stop$R,levels=1:2,labels=c("left","right"))
  dat.stop$C <- dat.stop$R==dat.stop$S 
  rbind(dat.go,dat.stop)[,c("D","S","SSD","R","C","RT")]
}

