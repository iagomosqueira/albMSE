###########################################################
# how to use MCMC results and CPP code to run projections #
###########################################################
# R. Hillary & I. Mosqueira (2024) ########################
###########################################################

library(Rcpp)
library(FLCore)
library(ggplotFL)
library(parallel)
library(mvtnorm)
source("utilities.R")

load('data/image/alb_abc_run6b.rda')

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# projection control file

prj.ctrl <- list(TAC=30000, nyproj=10)

# years of catch data to compute seasonal/fleet catch distribution

yref <- 2017:2020
cref <- C[as.character(yref),,]
cbar <- apply(cref,c(2,3),mean)
pcbar <- cbar/sum(cbar)

# Cproj [y, s, f]
Cproj <- array(dim=c(prj.ctrl$nyproj,ns,nf))
for(y in 1:prj.ctrl$nyproj) Cproj[y,,] <- prj.ctrl$TAC*pcbar

for(y in 1:prj.ctrl$nyproj) Cproj[y,,] <- prj.ctrl$TAC*pcbar

# abc.proj

abc.proj <- function(k) { # k is k-th MCMC iteration (parallel stuff)

  # mcvars LIST
  tmp <- mcvars[[k]]
  # mcpars matrix
  # parfv [npar + 3]: R0, ... selvx, , h, M, sigmaR
  parv <- mcpars[k,]

  # selex
  selvx <- exp(parv[(ny+2):npar])
  selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
  sela <- get.sel.age(6,5,selidx,selpars)
  # SRR h
  hh <- parv[npar+1]
  # M
  M <- parv[npar+2]
  # SRR sigmaR
  sigmarx <- parv[npar+3]
  nvec <- as.vector(mcvars[[k]]$N[dim(mcvars[[k]]$N)[1],,,]) # initialise in final year conditioning
  cvec <- as.vector(Cproj)
  R0 <- exp(parv[1])
  B0 <- tmp$SSB[1]/tmp$dep[1]
  spr0 <- B0/R0
  dms <- c(prj.ctrl$nyproj,ns,na,nbins,nf)
  epsrx <- rnorm(dms[1],0,sigmarx) # generate future rec. devs.
 
  # dms 
  onp <- list(dms,srec,R0,hh,psi,epsrx,spr0,M,as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fcpue)

  resx <- pdynlfcpue(dms,srec,R0,hh,psi,epsrx,spr0,M,as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fcpue)
  
  resx$inp <- list(dms,srec,R0,hh,psi,epsrx,spr0,M,as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fcpue)

  resx$dep <- resx$S/B0

  return(resx)

}

dms
srec
R0
hh
psi
epsrx
spr0
M
as.vector(mata)
as.vector(wta)
as.vector(sela)
nvec
cvec
as.vector(pla)
fcpue

# abc.om(om, ctrl)
# dim(om) [a, y, 1, 4, 1, i]
# [q, y, u, s, a, i]

# INPUT:
# OUTPUT: S (ssb), N (stock.n), H (hr), I (index), LF (lf), dep (depletion)

tes <- abc.proj(1)

reszz <- mclapply(1:nits,abc.proj,mc.cores=ncore)

dep.prj <- S.prj <- array(dim=c(nits,prj.ctrl$nyproj,ns))
for(k in 1:nits) {

  dep.prj[k,,] <- matrix(reszz[[k]]$dep,nrow=prj.ctrl$nyproj,ncol=ns)
  S.prj[k,,] <- matrix(reszz[[k]]$dep,nrow=prj.ctrl$nyproj,ncol=ns)

}

depx <- dep.prj[,,srec-1]
ssbx <- S.prj[,,srec-1]

dq <- apply(depx,2,quantile,c(0.025,0.5,0.975))
sq <- apply(ssbx,2,quantile,c(0.025,0.5,0.975))

plot(2020:2029,dq[2,],xlab='year',ylab='Depletion',ylim=c(0,max(dq[3,])),lwd=1.5,type='l',col='blue')
lines(2020:2029,dq[1,],lty=2,lwd=1.5,col='blue')
lines(2020:2029,dq[3,],lty=2,lwd=1.5,col='blue')

# glue it all together

dep.all <- ssb.all <- array(dim=c(nits,ny+prj.ctrl$nyproj-1))
for(k in 1:nits) {

  dep.all[k,] <- c(mcvars[[k]]$dep,depx[k,-1])
  ssb.all[k,] <- c(mcvars[[k]]$S,ssbx[k,-1])

}

dtq <- apply(dep.all,2,quantile,c(0.025,0.5,0.975))
stq <- apply(ssb.all,2,quantile,c(0.025,0.5,0.975))

mlab <- paste("TAC:",prj.ctrl$TAC,sep=" ")
plot(2000:2029,dtq[2,],xlab='year',ylab='Depletion',ylim=c(0,max(dtq[3,])),lwd=1.5,type='l',col='blue',main=mlab)
lines(2000:2029,dtq[1,],lty=2,lwd=1.5,col='blue')
lines(2000:2029,dtq[3,],lty=2,lwd=1.5,col='blue')
abline(h=0,lty=2)

