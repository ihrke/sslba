# Parameters to psoptim
#
# MG: par   - vector: length defines number of dimensions of the problem
#                     the dimensionthe parameters of the function
#                     NOTE: In optim, names of variables get passed through
#                     to fn and gr, but not their values.
#
# MG: fn    - function: the function to be optimised 
#                       (by default minimised, but see fnscale)
# MG: gr    - function: a gradient function (for cases where PSO passes to optim)
# MG: lower - vector?: lower bounds on variables
# MG: upper - vector?: upper bounds on variables
# MG: control - control parameters:
#             - fnscale - "SPSO2007" (default) or "SPSO2011" 
#             - type -
#             - maxit -   maximum number of iterations (default 1000)
#             - maxf -
#             - abstol -
#             - reltol -
#             - s -       Swarm size (defaults to floor(10+2*sqrt(length(par))) 
#                         unless type is "SPSO2011" in which case default is 40
#             - k -       exponent for calculating number of informants (default 3)
#             - p -       average % of informants for each particles (1 = fully informed)
#             - w -       exploitation constant (1 or 2 values)
#             - c.p -     local exploration constant
#             - c.g -     global exploration constant
#             - d -
#             - v.max -
#             - rand.order -
#             - max.restart -
#             - maxit.stagnate -
#             - vectorize -
#             - hybrid -
#             - hybrid.control -
#             - REPORT -
#             - trace   - integer: if positive, tracing information is produced
#             - trace.stats -
#
# NOTE: p is the structure into which all arguments are copied
#

# par is a vector or matrix, initial pars for 1 
# source ("psoptimcl.R")
# guesses <- c( 1,  2,  3,  4,
#               5,  6,  7,  8,
#               9, 10, 11, 12)
#  
# P <- t(matrix (guesses, ncol=4, nrow=3, byrow=TRUE))
# out <- psoptimcl (par=P, fn=objfunc, lower=0, upper=18)
#  
# For convenience, I've placed the parameters in columns and different particles in rows 
# (e.g. particle 1 parameters above are 1,2,3,4). This requires transposing the matrix before 
# passing it to psoptimcl (let me know if I should be doing this differently).
#  
# You can also specify a vector of lower and upper bounds if different parameters 
# require different bounds. There's additional parameters you can use (refer to the 
# psoptim package documentation). But for example, to maximise a function using 
# SPSO2011 with 80 particles and a maximum of 1500 iterations:
#  
# out <- psoptimcl (par=P, fn=objfunc, lower=0, upper=8, 
#   control=list (type="SPSO2011", s=80, maxit=1500, fnscale=-1))
# 

U_psoptimcl <- function (par, fn, gr = NULL, ..., lower=-1, upper=1,
                       control = list()) 
{
    # MG: Scales the function by a fnscale parameter (as per optim)
    fn1 <- function(par) fn(par, ...)/p.fnscale
  
    # MG: Creates an n x m matrix containing uniform randoms with upper and lower bounds 
    mrunif <- function(n,m,lower,upper) 
    {
        return(matrix(runif(n*m,0,1),nrow=n,ncol=m)*(upper-lower)+lower)
    }
  
    # MG: ??  Element-wise sqrt(sum(x^2)) Used by rsphere.unif 
    norm <- function(x) sqrt(sum(x*x))
  
    # MG: ...? Spherical uniform random? Used by??? Is this to do with particle topology?
    rsphere.unif <- function(n,r) 
    {
        temp <- runif(n)
        return((runif(1,min=0,max=r)/norm(temp))*temp)
    }
  
    # MG: ?? Create a vector with n elements of a, replace element b with k
    svect <- function(a,b,n,k) 
    {
        temp <- rep(a,n)
        temp[k] <- b
        return(temp)
    }
  
    # MG...? Matrix form of rsphere.unif ?
    mrsphere.unif <- function(n,r) 
    {
        m <- length(r)
        temp <- matrix(runif(n*m),n,m)
        return(temp%*%diag(runif(m,min=0,max=r)/apply(temp,2,norm)))
    }

    # MAIN CODE
    # Initialise variables
    
    # MG: Modified to allow passing in multiple 'sets' of initial parameters
#     if ((length(par) %% num_init) != 0)
#     {
#       stop (paste("Number of parameters [", length(par), "] is not a multiple of the number of parameter sets passed in! [", num_init ,"]", sep=""))
#     }
    # MG: handle matrix and vector arguments
    if (is.matrix(par)) {
      npar <- dim(par)[1]
      num_init <- dim(par)[2]
    }
    else {
      npar <- length(par)
      num_init <- 1      
    }
      
    lower <- as.double(rep(lower, ,npar))  # MG: Is this if upper/lower are not specified for each variable?  
    upper <- as.double(rep(upper, ,npar))   
  
    # Set the default control attributes (and then replace with those passed in via 'control')
    con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
                abstol = -Inf, reltol = 0, REPORT = 10,
                s = NA, k = 3, p = NA, w = 1/(2*log(2)),
                c.p = .5+log(2), c.g = .5+log(2), d = NA,
                v.max = NA, rand.order = TRUE, max.restart=Inf,
                maxit.stagnate = Inf,
                vectorize=FALSE, hybrid = FALSE, hybrid.control = NULL,
                trace.stats = FALSE, type = "SPSO2007")
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

    ## Argument error checks
    if (any(upper==Inf | lower==-Inf))
        stop("fixed bounds must be provided")

    # Use a structure, p, to store parameters
    p.type <- pmatch(con[["type"]],c("SPSO2007","SPSO2011"))-1
    if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
  
    p.trace   <- con[["trace"]]>0L           # provide output on progress?
    p.fnscale <- con[["fnscale"]]            # scale funcion by 1/fnscale
    p.maxit   <- con[["maxit"]]              # maximal number of iterations
    p.maxf    <- con[["maxf"]]               # maximal number of function evaluations
    p.abstol  <- con[["abstol"]]             # absolute tolerance for convergence
    p.reltol  <- con[["reltol"]]             # relative minimal tolerance for restarting
    p.report  <- as.integer(con[["REPORT"]]) # output every REPORT iterations
    p.s       <- ifelse(is.na(con[["s"]]),ifelse(p.type==0,floor(10+2*sqrt(npar)),40),
                  con[["s"]])                # swarm size
    p.p       <- ifelse(is.na(con[["p"]]),1-(1-1/p.s)^con[["k"]],con[["p"]]) # average % of informants
    p.w0      <- con[["w"]]                  # exploitation constant
  
    if (length(p.w0)>1) 
    { 
        p.w1    <- p.w0[2]
        p.w0    <- p.w0[1]
    } 
    else 
    {
        p.w1    <- p.w0
    }
    p.c.p         <- con[["c.p"]]            # local exploration constant
    p.c.g         <- con[["c.g"]]            # global exploration constant
    p.d           <- ifelse(is.na(con[["d"]]),norm(upper-lower),con[["d"]]) # domain diameter
    p.vmax        <- con[["v.max"]]*p.d      # maximal velocity
    p.randorder   <- as.logical(con[["rand.order"]]) # process particles in random order?
    p.maxrestart  <- con[["max.restart"]]    # maximal number of restarts
    p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
    p.vectorize   <- as.logical(con[["vectorize"]]) # vectorize?
    if (is.character(con[["hybrid"]])) 
    {
        p.hybrid <- pmatch(con[["hybrid"]],c("off","on","improved"))-1
        if (is.na(p.hybrid)) stop("hybrid should be one of \"off\", \"on\", \"improved\"")
    }
    else
    {
        p.hybrid <- as.integer(as.logical(con[["hybrid"]])) # use local BFGS search
    }
    p.hcontrol    <- con[["hybrid.control"]]  # control parameters for hybrid optim
    if ("fnscale" %in% names(p.hcontrol))
        p.hcontrol["fnscale"] <- p.hcontrol["fnscale"]*p.fnscale
    else
        p.hcontrol["fnscale"] <- p.fnscale
    p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?
  
    if (p.trace) 
    {
        message("S=",p.s,", K=",con[["k"]],", p=",signif(p.p,4),", w0=",
                signif(p.w0,4),", w1=",
                signif(p.w1,4),", c.p=",signif(p.c.p,4),
                ", c.g=",signif(p.c.g,4))
        message("v.max=",signif(con[["v.max"]],4),
                ", d=",signif(p.d,4),", vectorize=",p.vectorize,
                ", hybrid=",c("off","on","improved")[p.hybrid+1])
        if (p.trace.stats) 
        {
            stats.trace.it <- c()
            stats.trace.error <- c()
            stats.trace.f <- NULL
            stats.trace.x <- NULL
        }
    }

    ## Initialization (MG: ?? if PSO)
    if (p.reltol!=0) p.reltol <- p.reltol*p.d
    if (p.vectorize) 
    {
        lowerM <- matrix(lower,nrow=npar,ncol=p.s)
        upperM <- matrix(upper,nrow=npar,ncol=p.s)
    }
      
    # MG: Old code
    # Copy ALL parameters into matrix X[,1] only if they're in range  
    # These become the initial positions for particle 1
    # Why only particle 1??
    # X <- mrunif(npar,p.s,lower,upper)    
    # if (!any(is.na(par)) && all(par>=lower) && all(par<=upper)) X[,1] <- par    
    
    # MG: New code
    if (!any(is.na(par)) && all(par>=lower) && all(par<=upper)) {
          
      # MG: Make a matrix with the input parameter sets 
      X <- matrix (par, npar, num_init)
      
      if (num_init > p.s) {
        # MG: Could just increase p.s?
        stop (paste("Passed more parameter sets [", num_init, "] than there are particles to initialise [", p.s, "]\n", sep=""))
      }
    
      # MG: Pad to nparticles length with noise
      X <- cbind (X, mrunif(npar, p.s-num_init, lower, upper))
      
#       cat ("Number of parameters:", npar, "\n")
#       cat ("Number of initialised particles:", num_init, "\n")
#       cat ("Total num of parameters:", length(X),"\n")
#       for (z in 1:dim(X)[1]) {
#         cat ("Parameter matrix:", X[z,], "\n")
#       }
    }
    else
    {
      # MG: Initialise with random parameters
      X <- mrunif(npar,p.s,lower,upper)          
    }

    if (p.type==0) 
    {
        V <- (mrunif(npar,p.s,lower,upper)-X)/2
    } 
    else 
    { 
        ## p.type==1
        V <- matrix(runif(npar*p.s,min=as.vector(lower-X),max=as.vector(upper-X)),npar,p.s)
        p.c.p2 <- p.c.p/2 # precompute constants
        p.c.p3 <- p.c.p/3
        p.c.g3 <- p.c.g/3
        p.c.pg3 <- p.c.p3+p.c.g3
    }
  
    if (!is.na(p.vmax)) 
    { 
        # scale to maximal velocity
        temp <- apply(V,2,norm)
        temp <- pmin.int(temp,p.vmax)/temp
        V <- V%*%diag(temp)
    }
  
    # First iteration
    f.x <- apply(X, 2, fn1) # first evaluations (MG: ! Call to fn1)
    stats.feval <- p.s
    P <- X
    f.p <- f.x
    P.improved <- rep(FALSE,p.s)
    i.best <- which.min(f.p)
    error <- f.p[i.best]
    init.links <- TRUE
  
    if (p.trace && p.report==1) 
    {
        message("It 1: fitness=",signif(error,4))
        if (p.trace.stats) 
        {
            stats.trace.it <- c(stats.trace.it,1)
            stats.trace.error <- c(stats.trace.error,error)
            stats.trace.f <- c(stats.trace.f,list(f.x))
            stats.trace.x <- c(stats.trace.x,list(X))
        }
    }
  
    # Iteration loop
    stats.iter <- 1
    stats.restart <- 0
    stats.stagnate <- 0
  
#     cat ("stats.iter:", stats.iter, "\n")
#     cat ("p.maxit:", p.maxit, "\n")
#     cat ("stats.feval:", stats.feval, "\n")
#     cat ("p.maxf:", p.maxf, "\n")
#     cat ("error:", error, "\n")
#     cat ("p.abstol:", p.abstol, "\n")
#     cat ("stats.restart:", stats.restart, "\n")
#     cat ("p.maxrestart:", p.maxrestart, "\n")
#     cat ("stats.stagnate:", stats.stagnate, "\n")
#     cat ("p.maxstagnate) :", p.maxstagnate, "\n")

    # OUTER LOOP
    while (stats.iter<p.maxit && stats.feval<p.maxf && error>p.abstol &&
           stats.restart<p.maxrestart && stats.stagnate<p.maxstagnate) 
    {
        stats.iter <- stats.iter+1
        if (p.p!=1 && init.links) 
        {
            links <- matrix(runif(p.s*p.s,0,1)<=p.p,p.s,p.s)
            diag(links) <- TRUE   # diag - replace matrix diagonals
        }
      
        ## The swarm moves
        if (!p.vectorize) 
        {
            if (p.randorder) 
            {
                index <- sample(p.s)
            } 
            else
            {
                index <- 1:p.s
            }
          
            # MG: INNER LOOP
            for (i in index) 
            {
                # Communication between particles - set j
                if (p.p==1) j <- i.best
                      else  j <- which(links[,i])[which.min(f.p[links[,i]])] # best informant

                temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
                V[,i] <- temp*V[,i] # exploration tendency

                if (p.type==0) # SPSO 2007 
                {
                    V[,i] <- V[,i]+runif(npar,0,p.c.p)*(P[,i]-X[,i]) # exploitation
                    if (i!=j) V[,i] <- V[,i]+runif(npar,0,p.c.g)*(P[,j]-X[,i])
                } 
                else 
                { 
                    # SPSO 2011
                    if (i!=j) temp <- p.c.p3*P[,i]+p.c.g3*P[,j]-p.c.pg3*X[,i] # Gi-Xi
                         else temp <- p.c.p2*P[,i]-p.c.p2*X[,i] # Gi-Xi for local=best
                    V[,i] <- V[,i]+temp+rsphere.unif(npar,norm(temp))
                }
                # Finished using j (end of inter-particle chattery)
                
                if (!is.na(p.vmax)) 
                {
                    temp <- norm(V[,i])
                    if (temp>p.vmax) V[,i] <- (p.vmax/temp)*V[,i]
                }
                X[,i] <- X[,i]+V[,i]
        
                ## Check bounds
                temp <- X[,i]<lower
        
                if (any(temp)) 
                {
                    X[temp,i] <- lower[temp]
                    V[temp,i] <- 0
                }
      
                temp <- X[,i]>upper
        
                if (any(temp)) 
                {
                    X[temp,i] <- upper[temp]
                    V[temp,i] <- 0
                }
      
                ## Evaluate function
                if (p.hybrid==1) 
                {
                    # MG: CALLS OPTIM (if hybrid)
                    temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower, upper=upper,control=p.hcontrol)
                    V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
                    X[,i] <- temp$par
                    f.x[i] <- temp$value
                    stats.feval <- stats.feval+as.integer(temp$counts[1])
                } 
                else
                {
                    # MG: EVALUATE FUNCTION
                    f.x[i] <- fn1(X[,i])  # MG: Calls fn1
                    stats.feval <- stats.feval+1
                }
      
                if (f.x[i]<f.p[i]) 
                { 
                    # improvement
                    P[,i] <- X[,i]
                    f.p[i] <- f.x[i]
                    if (f.p[i]<f.p[i.best])
                    {
                        i.best <- i
                        if (p.hybrid==2) 
                        {
                            # MG: CALLS OPTIM
                            temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower,
                                          upper=upper,control=p.hcontrol)
                            V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
                            X[,i] <- temp$par
                            P[,i] <- temp$par
                            f.x[i] <- temp$value
                            f.p[i] <- temp$value
                            stats.feval <- stats.feval+as.integer(temp$counts[1])
                        }
                    }
                }
                if (stats.feval>=p.maxf) break
            }
            # MG: END INNER LOOP
        } 
        else 
        {
            # MG: If p.vectorize          
            if (p.p==1) j <- rep(i.best,p.s)
                  else  j <- sapply(1:p.s,function(i) # best informant
        
            which(links[,i])[which.min(f.p[links[,i]])]) 
            temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
            V <- temp*V # exploration tendency
      
            if (p.type==0) 
            {
                V <- V+mrunif(npar,p.s,0,p.c.p)*(P-X) # exploitation
                temp <- j!=(1:p.s)
                V[,temp] <- V[,temp]+mrunif(npar,sum(temp),0,p.c.p)*(P[,j[temp]]-X[,temp])
            }
            else 
            { 
                # SPSO 2011
                temp <- j==(1:p.s)
                temp <- P%*%diag(svect(p.c.p3,p.c.p2,p.s,temp))+
                P[,j]%*%diag(svect(p.c.g3,0,p.s,temp))-X%*%diag(svect(p.c.pg3,p.c.p2,p.s,temp)) # G-X
                V <- V+temp+mrsphere.unif(npar,apply(temp,2,norm))
            }
    
            if (!is.na(p.vmax)) 
            {
                temp <- apply(V,2,norm)
                temp <- pmin.int(temp,p.vmax)/temp
                V <- V%*%diag(temp)
            }
            X <- X+V
    
            ## Check bounds
            temp <- X<lowerM
            if (any(temp)) 
            {
                X[temp] <- lowerM[temp] 
                V[temp] <- 0
            }
            temp <- X>upperM
            if (any(temp)) 
            {
                X[temp] <- upperM[temp]
                V[temp] <- 0
            }
    
            ## Evaluate function
            if (p.hybrid==1) # not really vectorizing 
            { 
                # MG: Loop through (what?) calling optim
                for (i in 1:p.s) 
                {
                    temp <- optim(X[,i],fn,gr,...,method="L-BFGS-B",lower=lower, upper=upper,control=p.hcontrol)
                    V[,i] <- V[,i]+temp$par-X[,i] # disregards any v.max imposed
                    X[,i] <- temp$par
                    f.x[i] <- temp$value
                    stats.feval <- stats.feval+as.integer(temp$counts[1])
                }
            } 
            else 
            {
                f.x <- apply(X,2,fn1)
                stats.feval <- stats.feval+p.s
            }
          
            temp <- f.x<f.p
            if (any(temp)) 
            { 
                # improvement
                P[,temp]  <- X[,temp]
                f.p[temp] <- f.x[temp]
                i.best    <- which.min(f.p)
                if (temp[i.best] && p.hybrid==2) 
                { 
                    # overall improvement
                    temp <- optim(X[,i.best],fn,gr,...,method="L-BFGS-B",
                                  lower=lower, upper=upper,control=p.hcontrol)
                    V[,i.best] <- V[,i.best]+temp$par-X[,i.best] # disregards any v.max imposed
                    X[,i.best] <- temp$par
                    P[,i.best] <- temp$par
                    f.x[i.best] <- temp$value
                    f.p[i.best] <- temp$value
                    stats.feval <- stats.feval+as.integer(temp$counts[1])
                }
            }
            if (stats.feval>=p.maxf) break
        }
        
        if (p.reltol!=0) 
        {
            d <- X-P[,i.best]
            d <- sqrt(max(colSums(d*d)))
            if (d<p.reltol) 
            {
                X <- mrunif(npar,p.s,lower,upper)
                V <- (mrunif(npar,p.s,lower,upper)-X)/2
                if (!is.na(p.vmax)) 
                {
                    temp <- apply(V,2,norm)
                    temp <- pmin.int(temp,p.vmax)/temp
                    V <- V%*%diag(temp)
                }
                stats.restart <- stats.restart+1
                if (p.trace) message("It ",stats.iter,": restarting")
            }
        }
    
        init.links <- f.p[i.best]==error # if no overall improvement
        stats.stagnate <- ifelse(init.links,stats.stagnate+1,0)
        error <- f.p[i.best]
    
        if (p.trace && stats.iter%%p.report==0) 
        {
            if (p.reltol!=0) message("It ",stats.iter,": fitness=",signif(error,4), ", swarm diam.=",signif(d,4))
                        else message("It ",stats.iter,": fitness=",signif(error,4))
            
            if (p.trace.stats) 
            {
                stats.trace.it <- c(stats.trace.it,stats.iter)
                stats.trace.error <- c(stats.trace.error,error)
                stats.trace.f <- c(stats.trace.f,list(f.x))
                stats.trace.x <- c(stats.trace.x,list(X))
            }
        }
    }
  
    if (error<=p.abstol) 
    {
        msg <- "Converged"
        msgcode <- 0
    } 
    else if (stats.feval>=p.maxf) 
    {
        msg <- "Maximal number of function evaluations reached"
        msgcode <- 1
    } 
    else if (stats.iter>=p.maxit) 
    {
        msg <- "Maximal number of iterations reached"
        msgcode <- 2
    }
    else if (stats.restart>=p.maxrestart) 
    {
        msg <- "Maximal number of restarts reached"
        msgcode <- 3
    }
    else 
    {
        msg <- "Maximal number of iterations without improvement reached"
        msgcode <- 4
    }
    if (p.trace) message(msg)
  
    o <- list(par=P[,i.best],value=f.p[i.best],
              counts=c("function"=stats.feval,"iteration"=stats.iter, "restarts"=stats.restart),
              convergence=msgcode,message=msg) 
    if (p.trace && p.trace.stats) o <- c(o,list(stats=list(it=stats.trace.it,
                                                  error=stats.trace.error,
                                                  f=stats.trace.f,
                                                  x=stats.trace.x)))
    return(o)
}

library(compiler)
psoptimcl <- cmpfun (U_psoptimcl)

require("DEoptim")

fit.one <- function(p,dat,fn,type="nlm",control=list()) {
  if (type=="de" & !any(is.na(p))) {
    if (is.null(control$NP)) 
      NP <- 10*length(lower) else
      NP <- control$NP
    control$initialpop <- matrix(ncol=length(lower),nrow=NP)
    if (!is.matrix(p)) {
      if (length(p)!=length(lower)) stop("Wrong start length")
      ntodo <- 2
      control$initialpop[1,] <- p
    } else {
      if (dim(p)[1]!=length(lower) | dim(p)[2]>NP) stop("Wrong start dim")
      ntodo <- dim(p)[2]+1
      control$initialpop[1:(ntodo-1),] <- p
    }
    if (ntodo > NP) for (i in ntodo:NP) 
      control$initialpop[i,] <- runif(length(lower),lower,upper)
  }
  out <- switch(type,
                nlm=nlm(p=p,f=fn,dat=dat),
                pso=psoptimcl(par=matrix(p,ncol=1),fn=fn,dat=dat,
                              lower=lower,upper=upper,control=control),
                simplex=optim(par=p,fn=fn,dat=dat,control=control),
                de={
                  tmp <- DEoptim.control()
                  for (i in names(control)) tmp[[i]] <- control[[i]]
                  DEoptim(fn=fn,dat=dat,lower=lower,upper=upper,control=tmp)
                }
  )
  if (type=="de") {
    out$par <- out$optim$bestmem
    out$value <- out$optim$bestval    
  }
    if (type=="nlm") {
    out$par <- out$estimate
    out$value <- out$minimum
  }
  out
}

# If mind<=0 fits only once

# p=start;dat=data;type="simplex";fn=objective;mind=10;control=list(trace=100,maxit=5)
fit.multi <- function(p,dat,fn,type="nlm",mind=.1,max.rep=100,control=list()) {
  fit <- fit.one(p,dat,fn,type,control)
  vals <- fit$value
  rep <- 1
  if (mind>0) repeat {
    tmp <-fit.one(fit$par,dat,fn,type,control)
    vals <- c(vals,tmp$value)
    dvalue <- fit$value-tmp$value
    if (dvalue>0) fit <- tmp
    rep <- rep+1
    if (dvalue<mind) break
    if (rep>=max.rep) break
  }  
  fit$nrep <- rep
  fit$vals <- vals
  fit
}
