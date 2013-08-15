source("ssLBA.R")  # race model dfs
source("objsim.R")
require(truncnorm,quietly=T)

############################  Postive LBA NODE FUNCTIONS

# MAKE SURE PNORM DOESNT PRODUCE FUNNY RESULTS
pnormP  <- function(x,mean=0,sd=1,lower.tail=T){ifelse(abs(x)<7,pnorm(x,mean,sd,lower.tail),ifelse(x<0,0,1))}

dnormP <- function(x,mean=0,sd=1,lower.tail=T){ifelse(abs(x)<7,dnorm(x,mean,sd),0)}

# Density for one unit
p1=function(z,A,b,v,sv) {
  # pdf for a single unit
  if (A<1e-10) # LATER solution 
    return( pmax(0,((b/z^2)*dnormP(b/z,mean=v,sd=sv))/pnormP(v/sv)) ) # posdrifts
  zs=z*sv ; zu=z*v ; bminuszu=b-zu
  bzu=bminuszu/zs ; bzumax=(bminuszu-A)/zs
  pmax(0,((v*(pnormP(bzu)-pnormP(bzumax)) +
             sv*(dnormP(bzumax)-dnormP(bzu)))/A)/pnormP(v/sv)) # posdrifts
}

# Cumulative density for one unit
c1=function(z,A,b,v,sv) {
  # cdf for a single unit
  if (A<1e-10) # LATER solution 
    return( pmin(1,pmax(0,(pnormP(b/z,mean=v,sd=sv,lower.tail=F))/pnormP(v/sv))) ) # posdrifts
  zs=z*sv ; zu=z*v ; bminuszu=b-zu ; xx=bminuszu-A
  bzu=bminuszu/zs ; bzumax=xx/zs
  tmp1=zs*(dnormP(bzumax)-dnormP(bzu))
  tmp2=xx*pnormP(bzumax)-bminuszu*pnormP(bzu)
  pmin(pmax(0,(1+(tmp1+tmp2)/A)/pnormP(v/sv)),1) # posdrifts
}


############### Positive ssLBA RANDOM FUNCTION

# Simulate n trials using parameters in partlist
# SSD=NA simulates a go trial, othewise a stop trial
# If simulating stop trial first accumulator is assumed the stop accumulator
# and given index 0, choice accumulators given index 1, 2... on all trials
# For stop set SSD (one number or one per trial) to simulate SSD
# stoprt: save stoprt (else set to NA)
rfun <- function(n,parlist,SSD=NA,stoprt=FALSE) {
  if (any(is.na(SSD))) { # go trial
    SSD <- 0
    go <- TRUE
  } else go <- FALSE
  # set up SSD
  if (length(SSD)==1) SSD <- rep(SSD,n)
  if (length(SSD)!=n) stop("SSD wrong length")
  
  # Accumulators in order stop, response 1, response 2 ...
  nacc <- length(parlist$v)
  vs <- matrix(rtruncnorm(a=0,n=nacc*n,mean=parlist$v,sd=parlist$sv),nrow=nacc) 
  rts <- (parlist$b-runif(nacc*n,0,parlist$A))/vs # divide by v
  if (!go) {
    rts[1,] <- rts[1,]+parlist$ter[1]+SSD  # stop unit's ter + SSD
    rts[-1,] <- rts[-1,]+parlist$ter[2]  # go unit's ter
  } else rts <- rts+parlist$ter[1]  # go unit's ter
  # pick winner and save in row1=winner index, row2=winner rt matrix  
  winner <- apply(rts,2,function(x){which.min(x)})
  out <- rbind(winner,rts[cbind(winner,1:n)])
  gf <- as.logical(rbinom(dim(out)[2],1,parlist$gf))
  out[,gf] <- NA
  if (!go) {
    tf <- as.logical(rbinom(dim(out)[2],1,parlist$tf))
    
    if (!stoprt) out[2, !tf & winner==1] <- NA
    out[1,] <- out[1,]-1    
  }
  t(out)
}

