# testing seasonality

ns <- 4
na <- 10
ny <- 10
hh <- 0.75
mat <- rep(0,na)
mat[-c(1:5)] <- 1
M <- 0.1
F <- 0.025
R <- 1000
rs <- 1

###################
# 1. eqm (no SRR) #
###################

N <- array(dim=c(ny,ns,na))

# initial year

for(s in 1:ns) {

  if(s < rs) N[1,s,1] <- 0
  if(s == rs) N[1,s,1] <- R
  if(s > rs) N[1,s,1] <- N[1,s-1,1]*exp(-M-F)

}

for(a in 2:na) {

  N[1,1,a] <- N[1,ns,a-1]*exp(-M-F)
  for(s in 2:ns) N[1,s,a] <- N[1,s-1,a]*exp(-M-F)

}

# years

for(y in 2:ny) {

  # rec

  for(s in 1:ns) {

    if(s < rs) N[y,s,1] <- 0
    if(s == rs) N[y,s,1] <- R
    if(s > rs) N[y,s,1] <- N[1,s-1,1]*exp(-M-F)

  } 

  for(a in 2:na) {

    N[y,1,a] <- N[y-1,ns,a-1]*exp(-M-F)
    for(s in 2:ns) N[y,s,a] <- N[y,s-1,a]*exp(-M-F)
  
  }
} 


################
# 2. eqm (SSR) #
################

ss <- ifelse(rs-1>0,rs-1,ns)

# SPR0

neqm <- array(dim=c(ns,na))

for(s in 1:ns) {

  if(s < rs) neqm[s,1] <- 0
  if(s == rs) neqm[s,1] <- 1
  if(s > rs) neqm[s,1] <- neqm[s-1,1]*exp(-M)

}

for(a in 2:na) {

  neqm[1,a] <- neqm[ns,a-1]*exp(-M)
  for(s in 2:ns) neqm[s,a] <- neqm[s-1,a]*exp(-M)

}

spr0 <- sum(neqm[ss,]*mat)

# SPRf

for(s in 1:ns) {

  if(s < rs) neqm[s,1] <- 0
  if(s == rs) neqm[s,1] <- 1
  if(s > rs) neqm[s,1] <- neqm[s-1,1]*exp(-M-F)

}

for(a in 2:na) {

  neqm[1,a] <- neqm[ns,a-1]*exp(-M-F)
  for(s in 2:ns) neqm[s,a] <- neqm[s-1,a]*exp(-M-F)

}

sprf <- sum(neqm[ss,]*mat)

R0 <- 100
B0 <- R0*spr0
alp <- 4*hh/(spr0*(1-hh))
bet <- (5*hh-1)/(B0*(1-hh))
Bbar <- max((alp*sprf-1)/bet,0)
Rbar <- Bbar/sprf

# initial year

SSB <- rep(0,ny)

for(s in 1:ns) {

  if(s < rs) N[1,s,1] <- 0
  if(s == rs) N[1,s,1] <- Rbar
  if(s > rs) N[1,s,1] <- N[1,s-1,1]*exp(-M-F)

}

for(a in 2:na) {

  N[1,1,a] <- N[1,ns,a-1]*exp(-M-F)
  for(s in 2:ns) N[1,s,a] <- N[1,s-1,a]*exp(-M-F)

}

SSB[1] <- sum(N[1,ss,]*mat)

# years

for(y in 2:ny) {

  # rec

  for(s in 1:ns) {

    if(s < rs) N[y,s,1] <- 0
    if(s == rs) N[y,s,1] <- alp*SSB[y-1]/(1+bet*SSB[y-1])
    if(s > rs) N[y,s,1] <- N[1,s-1,1]*exp(-M-F)

  } 

  for(a in 2:na) {

    N[y,1,a] <- N[y-1,ns,a-1]*exp(-M-F)
    for(s in 2:ns) N[y,s,a] <- N[y,s-1,a]*exp(-M-F)
  
  }

  SSB[y] <- sum(N[y,ss,]*mat) 

}

##############################################
# initial condition function solver function #
##############################################

nf <- 2

# input: catch share per fishery by season

cshare <- matrix(nrow=ns,ncol=nf)
cshare[] <- 1/nf

# parameters: F by season by fishery
# objective function: SPR ratio and catch share by fishery per season

Fv <- matrix(0.05,nrow=ns,ncol=nf)

sprfunc <- function(Fx) {

  Fs <- apply(Fx,1,sum)

  # SPR0

  neqm <- array(dim=c(ns,na))
  ceqm <- array(0,dim=c(ns,nf))

  for(s in 1:ns) {

    if(s < rs) neqm[s,1] <- 0
    if(s == rs) neqm[s,1] <- 1
    if(s > rs) neqm[s,1] <- neqm[s-1,1]*exp(-M)

  }

  for(a in 2:na) {

    neqm[1,a] <- neqm[ns,a-1]*exp(-M)
    for(s in 2:ns) neqm[s,a] <- neqm[s-1,a]*exp(-M)

  } 

  spr0 <- sum(neqm[ss,]*mat)

  # SPRf

  for(s in 1:ns) {

    if(s < rs) neqm[s,1] <- 0
    if(s == rs) neqm[s,1] <- 1
    if(s > rs) neqm[s,1] <- neqm[s-1,1]*exp(-M-Fs[s])
    for(f in 1:nf) ceqm[s,f] <- ceqm[s,f]+(Fv[s,f])/(M+Fv[s,f])*(1-exp(-M-Fs[s]))*neqm[s,1]

  }

  for(a in 2:na) {

    neqm[1,a] <- neqm[ns,a-1]*exp(-M-Fs[ns])
    for(s in 2:ns) neqm[s,a] <- neqm[s-1,a]*exp(-M-Fs[s-1])
    for(s in 2:ns) 
      for(f in 1:nf) ceqm[s,f] <- ceqm[s,f]+(Fv[s,f])/(M+Fv[s,f])*(1-exp(-M-Fs[s]))*neqm[s,a] 

  }

  sprf <- sum(neqm[ss,]*mat)

  # relative catch by fleet and season

  peqm <- t(apply(ceqm,1,function(x){x/sum(x)}))
 
  # SPR ratio 

  rho <- sprf/spr0

  return(res)

}

