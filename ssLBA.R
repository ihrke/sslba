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

# pl=pls$D1S1s;SSD=j
dfun.stop <- function(pl,SSD) {

  tmpf=function(t,pl,SSD) p1(z=pmax(t - pl$ter[1] - SSD,1e-5),A=pl$A[1],b=pl$b[1],v=pl$v[1],sv=pl$sv[1])* 
                          (1-c1(pmax(t - pl$ter[2],1e-5),A=pl$A[2],b=pl$b[2],v=pl$v[2],sv=pl$sv[2]))*
                          (1-c1(pmax(t - pl$ter[2],1e-5),A=pl$A[3],b=pl$b[3],v=pl$v[3],sv=pl$sv[3])) 
  pstop=try(integrate(f=tmpf,lower=pl$ter[1]+SSD,upper=Inf,pl=pl,SSD=SSD)$value)
  # debugging
  if (class(pstop)=="try-error") {
    print(SSD)
    print(pl$ter[1]+SSD)
    print(pl)
    lower.bad <<- pl$ter[1]+SSD; SSD.bad <<- SSD; pl.bad <<- pl 
    stop()
  }
  pmax(pmin(pstop,1),0) # protect from numerical errors
}


