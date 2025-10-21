# utilities.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/v2/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# mc.output {{{

mc.output <- function(x, C) {

  # DIMENSIONS
  nits <- length(x)
  dmns <- list(age=0:14, year=2000:2020, season=1:4, unit=c('F', 'M'))

  out <- list()

  # - FLQuant [a, y, u, s, 1, i]

  # stock.n - N (y, a, s, u)
  out$stock.n <- Reduce(combine, lapply(x, function(i)
   FLQuant(aperm(i$N, c(2,1,4,3)), dimnames=dmns, units='1000') / 1000
  ))

  # m - M
  out$m <- expand(FLQuant(unlist(lapply(x, '[[', 'M')),
    quant='age',dim=c(1,1,1,1,1,nits)),
    age=0:14, year=2000:2020, season=1:4, unit=c('F', 'M'))

  # index.hat - Ihat
  out$index.hat <- Reduce(combine, lapply(x, function(i)
   FLQuant(c(i$Ihat), dimnames=list(age='all', year=2000:2020, season=1:4))
  ))

  # hr - H [y, s, f]
  out$hr <- FLQuant(unlist(lapply(x, '[[', 'H')),
    dimnames=list(age='all', year=2000:2020, season=1:4, area=1:6,
    iter=seq(nits)), units='hr')
  out$hrs <- areaSums(out$hr)

  # catch.sel - sela (a, s, u, f)
  out$catch.sel <- Reduce(combine, lapply(x, function(i) {
    # a, u, s, f
    res <- FLQuant(c(aperm(i$sela, c(1,3,2,4))), dimnames=list(age=0:14,
      unit=c('F', 'M'), season=1:4, area=1:6))

    res <- expand(res, year=2000:2020, fill=TRUE)

    return(res)
    }
  ))

  # deviances - epsrx (y)
  out$deviances <- Reduce(combine, lapply(x, function(i)
   FLQuant(c(1, exp(i$epsrx)), dimnames=list(age=0, year=2000:2020, season=1:4,
      unit=c('F', 'M')), units='')))

  out$deviances[,,,1:3] <- NA

  # HR @age[a,y,s,f,u]
  out$hra <- expand(out$hr, age=0:14, unit=c('F', 'M'),
    fill=TRUE) %*% out$catch.sel

  # catches (y, s, f)
  caf <- FLQuant(dimnames=list(year=2018:2020, season=1:4, area=1:6))
  caf[] <- C[19:21,,]
  out$cap <- caf %/% areaSums(caf)
 
  # catch.n
  out$catch.n <- out$stock.n %*% out$hra

  # - FLPar

  # SB0
  SB0 <- unlist(lapply(x, '[[', 'B0'))

  # R0, value in thousands
  R0 <- unlist(lapply(x, '[[', 'R0')) / 1000

  # SBMSY
  Bmsy <- unlist(lapply(x, '[[', 'Bmsy'))

  # h [i]
  h <- unlist(lapply(x, '[[', 'h'))
  
  # srpars
  out$srpars <-  FLPar(v=SB0, R0=R0, s=h)
  
  # hmsy [s, i]
  Hmsy <- unlist(lapply(x, '[[', 'hmsy'))

  # Cmsy
  cmsy <- unlist(lapply(x, '[[', 'Cmsy'))

  # refpts
  out$refpts <- FLPar(NA, dimnames=list(params=c('SB0', 'R0', 'HRMSY', 'SBMSY',
    'MSY'), season=1:4, iter=seq(nits)))
  out$refpts$SB0[,4,] <- SB0
  out$refpts$R0[,4,] <- R0
  out$refpts$HRMSY <- Hmsy
  out$refpts$SBMSY[,4,] <- Bmsy
  out$refpts$MSY[,1,] <- cmsy

  # rho
  out$rho <- unlist(lapply(x, '[[', 'rho'))

  # rho
  out$sigmar <- unlist(lapply(x, '[[', 'sigmar'))

  # - FLQuant

  # SSB
  out$ssb <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$SSB, dimnames=list(age='all', year=2000:2020))
  ))

  # Rtot
  out$rec <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$Rtot, dimnames=list(age='0', year=2000:2020))
  ))

  # dep
  out$dep <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$dep, dimnames=list(age='all', year=2000:2020))
  ))

  # index.q
  lnq <- unlist(lapply(x, '[[', 'lnq'))
  out$index.q <- FLQuant(exp(lnq), dimnames=list(age='all', season=1:4,
    iter=seq(length(rho))))

  # stock.n, m, catch.n, catch.sel, ssb, dep, index.q, srpars, refpts,
  # hr, rec, index.hat, hra
  return(out)
}
# }}}

# buildOM (stk,vars) {{{

buildOM <- function(stk, vars) {

  # WINDOW
  stk <- window(stk, start=2000)

  # ASSIGN stock.n and m
  stock.n(stk) <- vars$stock.n
  m(stk) <- vars$m

  # SIMPLIFY

}
# }}}

# FUNCTIONS {{{

logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

# VECTORIZE
aget.sel.age <- function(nf=6,nselg=5,selidx=c(1,2,3,4,5,5),selpars) {

  seltmp <- rep(NA,20)
  sela <- array(dim=c(na,ns,2,nf))
  
  # mula, sdla [a,s,u], 

  almin <- array(pmax(0, mula - sdla * 1.96), dim=dim(mula))
  almax <- array(pmax(0, mula + sdla * 1.96), dim=dim(mula))
  # n, a, s, u
  alref <- aperm(array(c(mapply(function(x,y) seq(x, y, length=20),
    x=almin, y=almax, SIMPLIFY=TRUE)), dim=c(20, 15, 4, 2)),
    c(2,3,4,1))

  adl <- array(dnorm(c(alref), rep(mula, 20), rep(sdla, 20)),
    dim=c(dim(sdla), 20))
  adl <- adl / as.numeric(array(apply(adl, 1:3, sum, simplify=FALSE),
    dim=dim(adl)))

  bb <-   lapply(1:5, function(f) ifelse(lref < selpars[f,1],
    2^((lref - selpars[f,1]) / selpars[f,2] ^ 2),
    2^(-((lref - selpars[f,1]) / selpars[f,3])^2)))

        for(f in 1:nf) {
          fref <- selidx[f]
          for(l in 1:20) seltmp[l] <- ifelse(lref[l] < selpars[fref,1],2^{-((lref[l]-selpars[fref,1])/selpars[fref,2])^2},2^{-((lref[l]-selpars[fref,1])/selpars[fref,3])^2})

          sela[a,s,g,f] <- sum(seltmp*dl)
        }

  return(sela)
}


get.sel.age <- function(nf=6,nselg=5,selidx=c(1,2,3,4,5,5),selpars)
{

  seltmp <- rep(NA,20)
  sela <- array(dim=c(na,ns,2,nf))
  for(g in 1:2) {
    for(s in 1:4) {
      for(a in 1:na) {

        lmin <- max(0,mula[a,s,g]-sdla[a,s,g]*1.96)
        lmax <- mula[a,s,g]+sdla[a,s,g]*1.96
        lref <- seq(lmin,lmax,length=20)
        dl <- dnorm(lref,mula[a,s,g],sdla[a,s,g])
        dl <- dl/sum(dl)

        for(f in 1:nf) {
        
          fref <- selidx[f]
          for(l in 1:20) seltmp[l] <- ifelse(lref[l] < selpars[fref,1],2^{-((lref[l]-selpars[fref,1])/selpars[fref,2])^2},2^{-((lref[l]-selpars[fref,1])/selpars[fref,3])^2})

          sela[a,s,g,f] <- sum(seltmp*dl)
        }
      }
    }
  }

  return(sela)
}

objfn.init <- function(theta,targv,sela) {

  hxinit <- 1 / (1 + exp(-theta))
  resx <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
  as.vector(wta), as.vector(sela), hxinit) 

  cx <- resx$C
  px <- cx/sum(cx)
  tmpv <- c(resx$rho, as.vector(px))
  objv <- logit(tmpv)

  return(fnscale * (sum((objv - targv) ^ 2)))

}

msyfn <- function(H,ph,sela) {

  hx <- H * ph
  resx <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),hx)
    
  return(sum(resx$C))
}
# }}}

# mcmc.abc {{{
mcmc.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar)
  acp <- rep(0,ngibbs)

  # get initial guess discrepancy

  xx <- sim(R0,dep,h,M,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        wtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- parvecold

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc2.abc {{{
mcmc2.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+2)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
    
  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        } 

        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
            
      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      }  

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc2a.abc {{{
mcmc2a.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+2)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    } 

    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  if(length(yfmsy) == 1) {

    sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean)  

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 
       
      if(length(yfmsy) == 1) {

        sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc3.abc {{{
mcmc3.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+3)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmarold,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]  

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmarold,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # conditional posterior for sigmaR

    res.tmp <- sum(0.5*epsrx^2)
    atmp <- alpR+length(epsrx)/2
    btmp <- betR+res.tmp
    sigmarold <- sqrt(1/rgamma(1,atmp,btmp))

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold,sigmarold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc3a.abc {{{
mcmc3a.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+3)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  if(length(yfmsy) == 1) {

    sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmarold,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yfmsy,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      if(length(yfmsy) == 1) {

        sprior <- sprior+dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(hmsyrat,mufmsy,sdfmsy,TRUE)) 

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmarold,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # conditional posterior for sigmaR

    res.tmp <- sum(0.5*epsrx^2)
    atmp <- alpR+length(epsrx)/2
    btmp <- betR+res.tmp
    sigmarold <- sqrt(1/rgamma(1,atmp,btmp))

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold,sigmarold)

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc4.abc {{{
mcmc4.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+2)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    }

    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    if(qtrend) {

      resq <- log(I[,,fcpue]/(xx$I*qt))

    } else {

      resq <- log(I[,,fcpue]/xx$I) 

    } 

    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue.mcmc <- matrix(nrow=nits,ncol=prod(dim(resq)))

  dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
  dcpue <- sum(dcpue.loo)

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yof,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  if(length(yof) == 1) {

    zof <- max(hmsyrat-1,0)
    dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
    sprior <- sprior+dof

  } else {

    zof <- hmsyrat-1
    zof[zof < 0] <- 0
    dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
    sprior <- sprior+sum(dof) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        if(qtrend) {

          resq <- log(I[,,fcpue]/(xx$I*qt))

        } else {

          resq <- log(I[,,fcpue]/xx$I) 

        }

        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
      dcpue <- sum(dcpue.loo)

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yof,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean)  

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      if(length(yof) == 1) {

        zof <- max(hmsyrat-1,0)
        dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
        sprior <- sprior+dof

      } else {

        zof <- hmsyrat-1
        zof[zof < 0] <- 0
        dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
        sprior <- sprior+sum(dof) 

      } 
        

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) {

      theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold)
      dcpue.mcmc[(n-burn)/thin,] <- as.vector(dcpue.loo)

    }

  }

  return(list(pars=theta.mcmc,cpuelogl=dcpue.mcmc,acp=acp))
} 
# }}}

# mcmc5.abc {{{
mcmc5.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar+3)
  acp <- rep(0,ngibbs)
  acphm <- 0

  # get initial guess discrepancy

  xx <- sim(R0,dep,hold,Mold,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue.mcmc <- matrix(nrow=nits,ncol=prod(dim(resq))) 
  dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
  dcpue <- sum(dcpue.loo)

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]
  hmsy <- xx$Hmsy
  hy <- apply(xx$H[yof,,],c(1,2),sum)
  hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }   

  if(length(yof) == 1) {

    zof <- max(hmsyrat-1,0)
    dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
    sprior <- sprior+dof

  } else {

    zof <- hmsyrat-1
    zof[zof < 0] <- 0
    dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
    sprior <- sprior+sum(dof) 

  } 

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmarold,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    # resample (h,M) from pi(h,M)

    zval <- rbinom(1,1,acphmu)
    if(zval == 1) {
    
      xnew <- rmvnorm(1,c(hmu,Mmu),Sigma)
      hold <- xnew[1,1]
      Mold <- xnew[1,2]

    }

    # resample parameters conditional on (h,M)

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,hold,Mold,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))
 
      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue.loo <- dnorm(resq,0,sdcpue,TRUE)
      dcpue <- sum(dcpue.loo) 

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]  
      hmsy <- xx$Hmsy
      hy <- apply(xx$H[yof,,],c(1,2),sum)
      hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      if(length(yof) == 1) {

        zof <- max(hmsyrat-1,0)
        dof <- dnorm(0,0,sdof,TRUE)-dnorm(zof,0,sdof,TRUE)
        sprior <- sprior+dof

      } else {

        zof <- hmsyrat-1
        zof[zof < 0] <- 0
        dof <- dnorm(zof,0,sdof,TRUE)-dnorm(0,0,sdof,TRUE)
        sprior <- sprior+sum(dof) 

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmarold,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }
      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # conditional posterior for sigmaR

    res.tmp <- sum(0.5*epsrx^2)
    atmp <- alpR+length(epsrx)/2
    btmp <- betR+res.tmp
    sigmarold <- sqrt(1/rgamma(1,atmp,btmp))

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) {

      theta.mcmc[(n-burn)/thin,] <- c(parvecold,hold,Mold,sigmarold)
      dcpue.mcmc[(n-burn)/thin,] <- as.vector(dcpue.loo)
    }

  }

  return(list(pars=theta.mcmc,cpuelogl=dcpue.mcmc,acp=acp))
} 
# }}}

# sim {{{
sim <- function(R0=1e6, dep=0.5, h=0.75, M=0.075, selpars, epsr, dms, pctarg,selidx) {

  na <- dms[1]
  ns <- dms[2]
  nf <- dms[3]
  nselg <- dms[4]

  # SPR ratio at exploited eqm given steepness and depletion
  # (based on derivation: dep = (4*h*rho+h-1)/(5*h-1))
  rhotarg <- (dep*(5*h-1)+1-h)/(4*h)

  # create selectivity-at-age

  sela <- get.sel.age(nf,nselg,selidx,selpars) 
    
  # target vector (rho+pc)
  targv <- logit(c(rhotarg,as.vector(pctarg)))

  # wrapper objective function to solve (minimise)

  # Estimate initial Fs
  
  hinit <- array(0.15*pctarg, dim=c(ns, nf)) 

  # scalar to give objective function some bite

  theta <- logit(as.vector(hinit))

  zz <- optim(theta,objfn.init,targv=targv,sela=sela,method=("L-BFGS-B"),control=list(trace=0))

  hinit <- array(ilogit(zz$par), dim=c(ns, nf))
  resinit <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), as.vector(hinit)) 

  # relative H-split for MSY calcs
  ph <- as.vector(hinit[] / sum(hinit))

  msy <- optimise(msyfn,interval=c(0,0.9),ph=ph,sela=sela,maximum=TRUE)
  Hmsy <- msy$maximum
  Cmsy <- msy$objective
  resmsy <- msypdyn(c(ns,na,nf), srec, R0, h, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), Hmsy * ph)
  
  Bmsy <- resmsy$Bmsy
  spr0 <- resinit$spr0
  B0 <- R0*spr0
  alp <- 4*h/(spr0*(1-h))
  bet <- (5*h-1)/(B0*(1-h))

  Bratio <- Bmsy/B0
  Rratio <- resmsy$Rmsy/R0
  hmsyv <- apply(array(Hmsy*ph,dim=c(ns,nf)),1,sum)

  if(!all.equal(Rratio, (4*h*Bratio)/(h*(5*Bratio-1)+1-Bratio)))
    warning("B-H invariant check - should be same as Rratio")

  # set up initial numbers-at-age for input to population dynamics

  Rinit <- R0*(4*h*dep)/(h*(5*dep-1)+1-dep)
  Ninit <- array(resinit$N,dim=c(na,ns,2))
  Ninit[] <- Ninit[]*Rinit
  nvec <- as.vector(Ninit)

  # expected catch @ hinit
  zinit <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),
    as.vector(wta),as.vector(sela),as.vector(hinit))
  Cinit <- array(zinit$C,dim=c(ns,nf))

  # main population stuff

  cvec <- as.vector(C)

  ## generating predicted LF and CPUE 
    
  # fishery for CPUE generation

  resp2 <- pdynlfcpue(c(ny,ns,na,nbins,nf),srec,R0,h,psi,epsr,spr0,M,
    as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fcpue)

  N <- array(resp2$N,dim=c(ny,na,ns,2))
  S <- array(resp2$S,dim=c(ny,ns))
  H <- array(resp2$H,dim=c(ny,ns,nf))
  LFhat <- array(resp2$LF,dim=c(ny,nbins,ns,nf))
  Ihat <- array(resp2$I,dim=c(ny,ns))

  # predicted LF distro for relevant fisheries

  LFhat <- apply(LFhat,c(2,4),sum)
  phat <- LFhat[,flf]
  phat <- apply(phat,2,function(x){x <- x/sum(x)})

  return(list(N=N,S=S,H=H,LF=phat,I=Ihat,Bmsy=Bmsy,Cmsy=Cmsy,Hmsy=hmsyv,B0=B0))

}
# }}}

# get.mcmc.vars {{{

get.mcmc.vars <- function(parsmat) {

  varlist <- list()                    
  nnits <- dim(parsmat)[1]
  for(nn in 1:nnits) {

    R0x <- exp(mcpars[nn,1])
    depx <- ilogit(mcpars[nn,2])
    epsrx <- mcpars[nn,3:(ny+1)]
    selvx <- exp(mcpars[nn,(ny+2):npar])
    selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
    xx <- sim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

    varlist[[nn]] <- list()
    varlist[[nn]][['Rtot']] <- apply(xx$N[,1,srec,],1,sum)
    varlist[[nn]][['SSB']] <- xx$S[,srec-1]
    varlist[[nn]][['dep']] <- xx$S[,srec-1]/xx$B0
    varlist[[nn]][['dbmsy']] <- xx$S[,srec-1]/xx$Bmsy
    varlist[[nn]][['Cmsy']] <- xx$Cmsy
    hmsy <- xx$Hmsy
    hy <- apply(xx$H[,,],c(1,2),sum)
    hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 
    varlist[[nn]][['hmsyrat']] <- hmsyrat
    varlist[[nn]][['Ihat']] <- xx$I
    varlist[[nn]][['LFhat']] <- xx$LF 

    if(nn %% 100 == 0) cat("Iteration",nn,"of",nnits,"\n")
  }

  return(varlist)

}
# }}}

# get.mcmc2.vars {{{

get.mcmc2.vars <- function(mcpars) {

  nnits <- dim(mcpars)[1]
  varlist <- vector("list", nnits)
  for(nn in 1:nnits) {

    # R0, pars[,1]
    R0x <- exp(mcpars[nn,1])
    # dep, pars[,2]
    depx <- ilogit(mcpars[nn,2])
    # epsrsx, pars[3]
    epsrx <- mcpars[nn,3:(ny+1)]
    selvx <- exp(mcpars[nn,(ny+2):npar])
    selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
    hx <- mcpars[nn,npar+1]
    Mx <- mcpars[nn,npar+2]
    xx <- sim(R0x,depx,hx,Mx,selparsx,epsrx,dms,pctarg,selidx)
    # list(R0x,depx,hx,Mx,selparsx,epsrx,dms,pctarg,selidx)

    # N Rtot SSB dep dbmsy Cmsy hmsyrat H Ihat LFhat B0 R0 M h sela
    varlist[[nn]] <- setNames(nm=list("N", "Rtot", "SSB", "dep", "dbmsy",
      "Cmsy", "hmsy", "hmsyrat", "H", "Hy", "Ihat", "LFhat", "B0", "R0",
      "M", "h", "sela"))

    varlist[[nn]][['N']] <- xx$N
    varlist[[nn]][['Rtot']] <- apply(xx$N[,1,srec,],1,sum)
    varlist[[nn]][['SSB']] <- xx$S[,srec-1]
    varlist[[nn]][['dep']] <- xx$S[,srec-1]/xx$B0
    varlist[[nn]][['dbmsy']] <- xx$S[,srec-1]/xx$Bmsy
    varlist[[nn]][['Bmsy']] <- xx$Bmsy
    varlist[[nn]][['Cmsy']] <- xx$Cmsy
    # HMSY per season[s,f]
    hmsy <- xx$Hmsy
    varlist[[nn]][['hmsy']] <- hmsy
    # ADD over f
    hy <- apply(xx$H[,,],c(1,2),sum)
    varlist[[nn]][['hy']] <- hy
    hmsyrat <- apply(apply(hy,1,function(x,hmsy){x <- x/hmsy},hmsy),2,mean) 
    # Annual hrmsy
    varlist[[nn]][['hmsyrat']] <- hmsyrat 
    varlist[[nn]][['H']] <- xx$H
    varlist[[nn]][['Ihat']] <- xx$I
    resq <- log(I[,,fcpue] / xx$I)
    varlist[[nn]][['lnq']] <- apply(resq, 2, mean)
    varlist[[nn]][['LFhat']] <- xx$LF 
    varlist[[nn]][['B0']] <- xx$B0
    varlist[[nn]][['R0']] <- R0x
    varlist[[nn]][['M']] <- Mx
    varlist[[nn]][['h']] <- hx
    varlist[[nn]][['rho']] <- c(rho(FLQuant(epsrx)))
    varlist[[nn]][['sigmar']] <- sd(epsrx)
    varlist[[nn]][['epsrx']] <- epsrx
    varlist[[nn]][['sela']] <- get.sel.age(nf,nselg,selidx,selparsx) 

    if(nn %% 100 == 0) cat("Iteration",nn,"of",nnits,"\n")
  }

  return(varlist)

}
# }}}

# plot.mcmc.vars {{{
plot.mcmc.vars <- function(varlist,ptype) {

  nnits <- length(varlist)

  if(ptype == 'dep') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dep
    vmin <- 0
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmax <- max(vq[3,])
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='SSB depletion',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

  if(ptype == 'bmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dbmsy
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmin <- 0 
    vmax <- max(vq) 
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='Bmsy ratio',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

  if(ptype == 'hmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$hmsyrat
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmin <- 0
    vmax <- max(vq) 
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab=expression(H[y]/H[msy]),col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  } 

  if(ptype == 'rec') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$Rtot
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    vmin <- 0
    vmax <- max(vq) 
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='Recruitment',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

}
# }}}

# {{{ plot.mcmc.cpue
plot.mcmc.cpue <- function(varlist) {

  nnits <- length(varlist)

  vv <- array(dim=c(nnits,ny,ns))
   for(nn in 1:nnits) {
     
     tmpv <- varlist[[nn]]$Ihat
     iobs <- I[,,fcpue]
     if(seasonq) {

       if(qtrend) {

        resq <- log(iobs/(tmpv*qt))

       } else {

         resq <- log(iobs/tmpv) 

       }

       lnq <- apply(resq,2,mean)
       if(qtrend) {

         vv[nn,,] <- t(apply(tmpv*qt,1,function(x,lnq){x <- x*exp(lnq)},lnq))*rlnorm(ny*ns,0,sdcpue)

       } else {

         vv[nn,,] <- t(apply(tmpv,1,function(x,lnq){x <- x*exp(lnq)},lnq))*rlnorm(ny*ns,0,sdcpue) 

       }

     } else {

       if(qtrend) {

        resq <- log(iobs/(tmpv*qt))

       } else {

         resq <- log(iobs/tmpv) 

       } 

       lnq <- mean(resq)
       if(qtrend) {

         vv[nn,,] <- tmpv*qt*exp(lnq)*rlnorm(ny*ns,0,sdcpue)  

       } else {

         vv[nn,,] <- tmpv*exp(lnq)*rlnorm(ny*ns,0,sdcpue) 

       }

     }
   }   
   
   vq <- apply(vv,c(2,3),quantile,c(0.025,0.5,0.975))

   # ggplot the sumbitch

   vdf <- expand.grid(year=yrs,season=1:ns,obs=NA,hat=NA,lq=NA,uq=NA)
   vdf$obs <- as.vector(iobs)
   vdf$hat <- as.vector(vq[2,,])
   vdf$lq <- as.vector(vq[1,,]) 
   vdf$uq <- as.vector(vq[3,,]) 
   ggplot(vdf)+geom_line(aes(x=year,y=hat),colour='blue')+geom_line(aes(x=year,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=year,y=uq),colour='blue',linetype='dashed')+geom_point(aes(x=year,y=obs),colour='magenta')+facet_wrap(~season)+ylab("CPUE")+theme_bw()

}
# }}}

# {{{ plot.mcmc.lf
plot.mcmc.lf <- function(varlist) {

  nnits <- length(varlist)
  vv <- array(dim=c(nnits,nbins,nselg)) 
    for(nn in 1:nnits) vv[nn,,] <- varlist[[nn]]$LF
    vq <- apply(vv,c(2,3),quantile,c(0.025,0.5,0.975))
    vdf <- expand.grid(length=mulbins,fishery=1:nselg,obs=NA,hat=NA,lq=NA,uq=NA) 
    vdf$obs <- as.vector(pobs)
    vdf$hat <- as.vector(vq[2,,])
    vdf$lq <- as.vector(vq[1,,]) 
    vdf$uq <- as.vector(vq[3,,]) 
    ggplot(vdf)+geom_line(aes(x=length,y=hat),colour='blue')+geom_line(aes(x=length,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=length,y=uq),colour='blue',linetype='dashed')+geom_point(aes(x=length,y=obs),colour='magenta')+facet_wrap(~fishery)+ylab("Length frequency")+theme_bw()

}
# }}}

# {{{ plot.mcmc.sel 
plot.mcmc.sel <- function(mcpars) {

  nnits <- dim(mcpars)[1]
  selparsx <- exp(mcpars[,paridx[[3]]])
  msel <- array(dim=c(nnits,nbins,nselg))

  for(nn in 1:nnits) {
    for(f in 1:nselg) {
      
      sx50 <- selparsx[nn,1+(f-1)]
      sxL <- selparsx[nn,nselg+f]
      sxR <- selparsx[nn,2*nselg+f]

      for(l in 1:nbins) {
  
        lref <- mulbins[l]
        msel[nn,l,f] <- ifelse(lref<sx50,2^{-(lref-sx50)^2/(sxL^2)},2^{-(lref-sx50)^2/(sxR^2)})
        
      }
    }
  }

  vq <- apply(msel,c(2,3),quantile,c(0.025,0.5,0.975))
  vdf <- expand.grid(length=mulbins,fishery=1:nselg,med=NA,lq=NA,uq=NA) 
  vdf$med <- as.vector(vq[2,,])
  vdf$lq <- as.vector(vq[1,,])
  vdf$uq <- as.vector(vq[3,,]) 
  ggplot(vdf)+geom_line(aes(x=length,y=med),colour='blue')+geom_line(aes(x=length,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=length,y=uq),colour='blue',linetype='dashed')+facet_wrap(~fishery)+ylab("Size selectivity")+theme_bw()

}

# }}}

# get.mcmc.vars {{{

get.mcmc.vars <- function(varlist,vtype) {

  nnits <- length(varlist)

  if(vtype == 'dep') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dep
    vmin <- 0
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  if(vtype == 'bmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dbmsy
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  if(vtype == 'hmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$hmsyrat
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  if(vtype == 'rec') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$Rtot
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975)) 

  }

  return(round(vq,2))

} # }}}

# harvest {{{

setMethod('harvest', signature(object="FLmse", catch="missing"),
  function(object) {
    if(length(biol(om(object))) == 1)
      return(attr(om(object), "harvest")[[1]])
    else
      return(attr(om(object), "harvest"))
  })

setMethod('harvest', signature(object="FLombf", catch="missing"),
  function(object) {
    if(length(biol(object)) == 1)
      return(attr(object, "harvest")[[1]])
    else
  return(attr(object, "harvest"))
  })

setReplaceMethod('harvest', signature(object="FLombf", value="FLQuant"),
  function(object, biol=1, value) {
    attr(object, "harvest")[[biol]] <- value
    return(object)
})

setReplaceMethod('harvest', signature(object="FLombf", value="FLQuants"),
  function(object, value) {
    attr(object, "harvest") <- value
    return(object)
})

setMethod('fbar', signature(object="FLombf"),
  function(object, value) {
    areaSums(unitSums(seasonSums(quantMeans(harvest(om)))))
})

# }}}

# fwdabc.om {{{

fwdabc.om <- function(om, ctrl, pcbar, pla, verbose=FALSE, ...) {

  # DIMS
  dom <- dims(om)
  
  yrs <- ac(unique(ctrl$year))
  dyr <- ac(min(ctrl$year) - 1)
  yrs <- c(dyr, yrs)
  nyrs <- length(yrs)
  
  na <- dims(biol(om))$age
  nf <- length(fisheries(om))
  its <- seq(dims(om)$iter)

  # PROPAGATE ctrl
  if(dim(iters(ctrl))[3] == 1 & dims(om)$iter > 1) {
    ctrl <- propagate(ctrl, dims(om)$iter)
  }

  # CLEAN
  n(biol(om))[, yrs[-1]] <- as.numeric(NA)

  for(i in its) {

  # TODO: MODIFY inp every iter

  # ASSEMBLE inputs
  inp <-   list(
    # dms: dimensions: year, season, age, lengths, fishery
    dm_=c(nyrs, 4, 15, 27, 6),
    # srec: rec season
    srec_=4,
    # - R0: SRR R0 in numbers
    R0_=c(params(sr(biol(om)))$R0)[i] * 1000,
    # - hh: SRR H
    hh_=c(params(sr(biol(om)))$s)[i],
    # - psi: sex ratio at birth
    psi_=0.5,
    # - epsrx: future rec devs [year]
    # epsr_=rnorm(nyrs,0,0.355),
    epsr_=log(deviances(om)[, yrs,'F',4,,i]),
    # - spr0: SRR B0/R0
    spr0_=c(params(sr(biol(om)))$v[,i] / (params(sr(biol(om)))$R0[,i] * 1000)),
    # - M: M
    M_=iter(m(biol(om)), i)[[1]],
    # - mata: mat at age [age, season, unit (sex)] USE dyr as no change
    mata_=c(unname(aperm(iter(mat(biol(om)), i)[, dyr,,, 1,, drop=TRUE], c(1, 3, 2)))),
    # - wta: Wt at age (gr) [age, season, unit (sex)] USE dyr as no change
    wta_=c(unname(aperm(iter(wt(biol(om)), i)[, dyr,,, 1,, drop=TRUE], c(1, 3, 2)))) / 1000,
    # - sela: Selex at age [age, season, unit (sex), fishery]
    sela_=c(unname(aperm(abind(lapply(fisheries(om), function(x)
      catch.sel(x[[1]])[, dyr,,,,i]))[drop=TRUE], c(1,3,2,4)))),
    # - nvec: Last year Ns [age, season, unit (sex)], as vector.
    Ninit_=c(aperm(n(biol(om))[, dyr,,,,i, drop=TRUE], c(1,3,2))) * 1000,
    # - cvec: Catch in projection [year, season, fishery], as vector.
    Cb_=c(array(rep(pcbar, each=nyrs), dim=c(nyrs, 4, 6)) * 
      array(c(c(catch(om)[,dyr,,,,i]), iters(ctrl)[, 2, i]), dim=c(nyrs, 4, 6))),
    # - pla: ALK [lengths, age, season, unit (sex)]
    pla_=c(pla),
    # - fcpue: Index of fleet to generate CPUE
    fref_=1)

  # CALL pdynlfcpue, rei: S (ssb_fy), N, H, I, LF
  rei <- do.call(pdynlfcpue, inp)

  # EXTRACT n - N [y,a,s,u]
  n(biols(om)[[1]])[, yrs[-1],,,, i] <- 
    aperm(array(rei$N, dim=c(nyrs, na, 4, 2))[-1,,,, drop=FALSE], c(2,1,4,3)) / 1000

  # GET hr - H [y,s,f]
  hri <- divide(FLQuant(rei$H, dimnames=list(year=yrs, season=seq(4), 
    area=names(fisheries(om)))), 5)

  sei <- lapply(fisheries(om), function(x) iter(catch.sel(x[[1]]), i))

  hri <- Map(function(h, s) expand(h, unit=c('F', 'M')) %*%
      s[, yrs], h=hri, s=sei)

  # COMPUTE catch.n_f
  can <- lapply(hri, function(x) x %*% n(biols(om)[[1]])[, yrs,,,, i])

  # TEST: total catch
  # Reduce('+', lapply(can,
  #   function(x) quantSums(unitSums(seasonSums(x * wt(biol(om))[, yrs,,,, i])))))

  # ASSIGN as landings.n per fleet
  for(f in seq(6))
    fisheries(om)[[f]][[1]]@landings.n[, yrs[-1],,,, i] <- c(can[[f]][,-1])

  # TEST: total catch
  # Reduce('+', lapply(fisheries(om), function(x)
  #   unitSums(seasonSums(catch(x)[[1]][, yrs,,,,i]))))

  if(verbose)
    print(i)
  }
  return(list(om=om))
}

# }}}

# hr {{{

hrf <- function(om) {

  # C_f / sum(S_f * N * W)
  lapply(fisheries(om), function(x) {
    # ADD catch across sexes
    res <- unitSums(catch(x[[1]])) /
      # ADD N * WT * S_f across ages and units
      unitSums(quantSums(n(biol(om)) * wt(biol(om)) * catch.sel(x[[1]])))
    # SET to max = 0.9
    res <- ifelse(res > 0.9, 0.9, res)

    units(res) <- 'hr'
    
    return(res)
  })
}

hr <- function(om) {
  # 
  Reduce('+', hrf(om))
}

# }}} 

# cpuescore.ind {{{
# Calculated over mean and sd on refyrs

cpuescore.ind <- function(stk, idx, index = 1, refyrs = NULL, args, tracking) {
  # ARGS
  ay <- args$ay
  dlag <- args$data_lag
  dy <- ay - dlag
  # TODO: CHECK for frq>1

  # GET metric until dy
  met <- seasonMeans(window(biomass(idx[[index]])[1, ], end = dy))

  if (!is.null(refyrs)) {
    ref <- met[, ac(refyrs)]
  } else {
    ref <- met
  }

  ind <- FLQuants(zscore = (met[, ac(dy)] %-% yearMeans(ref)) %/%
    sqrt(yearVars(ref)))

  track(tracking, "cpue.ind", ac(ay)) <- ind$zscore

  return(list(stk = stk, ind = ind, tracking = tracking, cpue = met))
}
# }}}

# metrics {{{

# relhr: HR / HR_MSY
relhr <- function(x) {
  seasonMeans(hr(x) / refpts(x)$HRMSY)
}

# ssb: annual female spawning stock biomass in season 4
setMethod('ssb', signature(object="FLombf"),
  function(object) {
    return(quantSums(n(biol(object))[,,'F',4] * wt(biol(object))[,,'F',4] *
      mat(biol(object))[,,'F',4]))
  })

# tsb: annual total stock biomass in season 4
setMethod('tsb', signature(object="FLombf"),
  function(object) {
    return(unitSums(quantSums(n(biol(object))[,,,1] * wt(biol(object))[,,,4])))
  })

# catch: annual total catch in weight
setMethod('catch', signature(object="FLombf"),
  function(object) {
    unitSums(seasonSums(Reduce('+', catch(fisheries(object)))))
})

setMethod('msy', signature(x="FLombf"),
  function(x) {
    refpts(x)$MSY[,1,]
})

setMethod('fbar', signature(object="FLombf"),
  function(object) {
    relhr(object)
})

mets <- list(SB=ssb, C=catch, HR=relhr)

setMethod('depletion', signature(x='FLombf'),
  function(x)
  ssb(x) / ssb0(x)
  )

# }}}

# refpts {{{

ssb0 <- function(x) {
  refpts(x)$SB0[,4,]
}

setMethod('sbmsy', signature(x="FLombf"),
  function(x) {
    refpts(x)$SBMSY[,4,]
})

hrmsy <- function(x) {
  refpts(x)$HRMSY
}

setMethod('msy', signature(x="FLombf"),
  function(x) {
    refpts(x)$MSY[,1,]
})

# }}}

# statistics {{{

statistics <- list(
  # SB
  SB = list(~yearMeans(SB), name = "SB",
    desc = "Mean spawner biomass"),
  # SB0
  SB0 = list(~yearMeans(SB/SB0), name = "SB/SB[0]",
    desc = "Mean spawner biomass relative to unfished"),
  # minSB0
  minSB0 = list(~apply(SB/SB0, c(1, 3:6), min), name = "min(SB/SB[0])",
    desc = "Minimum spawner biomass relative to unfished"),
  # SBMSY
  SBMSY = list(~yearMeans(SB/SBMSY), name = "SB/SB[MSY]",
    desc = "Mean spawnwer biomass relative to SBMSY"),
  # HRMSY
  HRMSY = list(~yearMeans(HR), name = "HR",
    desc = "Mean annual relative harvest rate"),
  # green
  green = list(~yearSums(FLQuant((SB / SBMSY) > 1 & HR < 1)) / dim(SB)[2],
    name = "P(Green)", desc = "Probability of being in Kobe green quadrant"),
  # orange
  orange = list(~iterSums(FLQuant((SB / SBMSY) >= 1 & HR >= 1)) / dim(SB)[6],
    name = "P(Orange)", desc = "Probability of being in Kobe orange quadrant"),
  # yellow
  yellow = list(~iterSums(FLQuant((SB / SBMSY) < 1 & HR < 1)) / dim(SB)[6],
    name = "P(Yellow)", desc = "Probability of being in Kobe yellow quadrant"),
  # red
  red = list(~yearSums(FLQuant((SB / SBMSY) < 1 & HR > 1)) / dim(SB)[2],
    name = "P(Red)", desc = "Probability of being in Kobe red quadrant"),
  # PSBMSY
  PSBMSY = list(~yearMeans((SB / SBMSY) >= 1), name = "P(SB>=SB[MSY])",
    desc = "Probability of SB greater or equal to SBMSY"),
  # PSBlim
  PSBlim = list(~yearMeans((SB / (SB0 * 0.10)) > 1), name = "P(SB>SB[limit])", 
    desc = "Probability that spawner biomass is above 10% SB0"),
  # C
  C = list(~yearMeans(C), name = "mean(C)", desc = "Mean catch over years"),
  # C/MSY
  CMSY = list(~yearMeans(C/MSY), name = "C/MSY", desc = "Mean proportion of MSY"),
  # AAV
  AAVC = list(~yearMeans(abs(C[, -1] - C[, -dim(C)[2]]) / C[, -dim(C)[2]]),
    name = "AAV(C)", desc = "Average annual variability in catch"),
  # IACC
  IACC = list(~100 * yearSums(abs(C[, -1] - C[, -dim(C)[2]])) / 
    yearSums(C[, -dim(C)[2]]),
  name="IAC(C)", desc="Percentage inter-annual change in catch"),
  # PC0
  PC0 = list(~yearSums(C < 0.01 * MSY) / dim(C)[2], name = "P(shutdown)", 
    desc = "Probability of fishery shutdown")
  )
# }}}

# bufferdelta.hcr {{{

timeMeans <- function(x)
  seasonMeans(yearMeans(x))

zscore <- function(x, mean=yearMeans(x), sd=sqrt(yearVars(x)))
  (x %-% mean) %/% sd

bufferdelta.hcr <- function(stk, ind, target=0, metric='zscore',
  width=1, lim=target - 2 * width, sloperatio=0.20, initac=NULL,
  dupp=NULL, dlow=NULL, all=TRUE, ..., args, tracking) {

  # CONVERT parameters
  bufflow <- target - width
  buffupp <- target + width

  # EXTRACT args
  ay <- args$ay
  iy <- args$iy
  data_lag <- args$data_lag
  man_lag <- args$management_lag
  frq <- args$frq

  # SET data year
  dy <- ay - data_lag
  # SET control years
  cys <- seq(ay + man_lag, ay + man_lag + frq - 1)

  # COMPUTE metric
  met <- mse::selectMetric(metric, stk, ind)
  met <- window(met, start=dy, end=dy)
  
  # COMPUTE HCR multiplier if ...
  # BELOW lim
  hcrm <- ifelse(met <= lim, ((lim/met) ^ 2) / 2,
    # BETWEEN lim and bufflow
    ifelse(met < bufflow,
      (0.5 * (1 + (met - lim) / (bufflow - lim))),
    # BETWEEN bufflow and buffupp
    ifelse(met < buffupp, 1, 
    # ABOVE buffupp
      1 + sloperatio * 1 / (2 * (bufflow - lim)) * (met - buffupp))))

  # GET previous TAC from last hcr ...
  if(is.null(initac)) {
    pre <- tracking[[1]]['hcr', ac(ay)]
    # ... OR catch
    if(all(is.na(pre)))
      pre <- unitSums(areaSums(seasonSums(catch(stk)[, ac(ay - args$data_lag)])))
  } else {
    pre <- FLQuant(initac, iter=args$it)
  }

  # SET TAC as tac = B * (1 - exp(-fm * hcrm * (F / FMSY))
  out <- pre * hcrm

  # TRACK initial target
  track(tracking, "tac.hcr", cys) <- out

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
      out[out > pre * dupp] <- pre[out > pre * dupp] * dupp
    } else {
      out[out > pre * dupp & met < bufflow] <- pre[out > pre * dupp & met <
        bufflow] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
      out[out < pre * dlow] <- pre[out < pre * dlow] * dlow
    } else {
      out[out < pre * dlow & met < bufflow] <- pre[out < pre * dlow & met <
        bufflow] * dlow
    }
  }

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(cys, function(x) list(quant="catch", value=c(out), year=x)))
  )
	
  list(ctrl=ctrl, tracking=tracking)
}
# }}}

# buffer.hcr {{{

#' A buffered and non-linear hockeystick HCR
#'
#' A Harvest Control Rule (HCR) that extends the traditional hockeystick shape
#' by allowing for increasing output when stock status rises above a buffer set
#' around the target.
#' @param lim Point at which the HCR response moves from linear to quadratic in terms of reducing HCR multiplier, numeric.
#' @param bufflow Lower point of "buffer zone" where HCR response is fixed at 1.
#' @param buffupp Upper point of "buffer zone" where HCR response is fixed at 1.
#' @param sloperatio fractional difference
#'
#' so the response of the HCR in this case would be as follows:
#' if(muI > Ilim & muI <= buffl) HCRmult <- (0.5*(1+(muI-Ilim)/(buffl-Ilim)))
#' if(muI > buffl & muI < buffh) HCRmult <- 1
#' if(muI >= buffh) HCRmult <- 1+sloperatio*gr*(muI-buffh) 
#' if(muI <= Ilim) HCRmult <- (muI/Ilim)^2/2
#'
#' @author Original design by R. Hillary (CSIRO). Implemented by I. Mosqueira (WMR).
#' @references
#' Hillary, R. 2020. *Management Strategy Evaluation of the Broadbill Swordfish ETBF harvest strategies*. Working document.

buffer.hcr <- function(stk, ind, target, metric='depletion', lim=0.10,
  bufflow=0.30, buffupp=0.50, sloperatio=0.20, initac=NULL, dupp=NULL, dlow=NULL,
  all=TRUE, ..., args, tracking) {

  # EXTRACT args
  ay <- args$ay
  iy <- args$iy
  data_lag <- args$data_lag
  man_lag <- args$management_lag
  frq <- args$frq

  # SET data year
  dy <- ay - data_lag
  # SET control years
  cys <- seq(ay + man_lag, ay + man_lag + frq - 1)

  # COMPUTE metric
  met <- mse::selectMetric(metric, stk, ind)
  met <- window(met, start=dy, end=dy)

  # TRACK metric & status (HCR segment)
  track(tracking, "metric.hcr", dy) <- met
  track(tracking, "status.hcr", dy) <- cut(c(met),
    breaks=c(-1e-08, lim, bufflow, buffupp, Inf), c(1, 2, 3, 4))

  # COMPUTE HCR multiplier if ...
  # BELOW lim
  hcrm <- ifelse(met <= lim, ((met / lim) ^ 2) / 2,
    # BETWEEN lim and bufflow
    ifelse(met < bufflow,
      (0.5 * (1 + (met - lim) / (bufflow - lim))),
    # BETWEEN bufflow and buffupp
    ifelse(met < buffupp, 1, 
    # ABOVE buffupp
      1 + sloperatio * 1 / (2 * (bufflow - lim)) * (met - buffupp))))

  # GET previous TAC from last hcr ...
  if(is.null(initac)) {
    pre <- tracking[[1]]['hcr', ac(ay),,1]
    # ... OR catch
    if(all(is.na(pre)))
      pre <- unitSums(areaSums(seasonSums(catch(stk)[, ac(ay - args$data_lag)])))
  } else {
    pre <- FLQuant(initac, iter=args$it)
  }

  # SET TAC as tac = B * (1 - exp(-fm * hcrm * (F / FMSY))
  out <- target * hcrm

  # TRACK initial target
  track(tracking, "tac.hcr", cys) <- out

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
    out[out > pre * dupp] <- pre[out > pre * dupp] * dupp
    } else {
    out[out > pre * dupp & met < bufflow] <- pre[out > pre * dupp & met <
      bufflow] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
    out[out < pre * dlow] <- pre[out < pre * dlow] * dlow
    } else {
    out[out < pre * dlow & met < bufflow] <- pre[out < pre * dlow & met <
      bufflow] * dlow
    }
  }

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(cys, function(x) list(quant="catch", value=c(out), year=x)))
  )
	
  list(ctrl=ctrl, tracking=tracking)
}
# }}}

# plots {{{

plotOMRuns <- function(x, ...) {

  xms <- FLQuants(lapply(mets, function(f) window(f(x), end=2023)))

  args <- list(...)

  yms <- lapply(args, function(i)
    FLQuants(lapply(mets, function(f) window(f(om(y)), start=2023))))

  ms <- c(list(OM=xms), yms)

  plotListFLQuants(ms)
}

plotMetrics <- function(...) {

  args <- list(...)

  ms <- lapply(args, function(i) {
    FLQuants(lapply(mets, function(f) f(i)))
  })

  ms <- lapply(ms, function(i) {
    names(i) <- c("SSB (t)", "C (t)", "H/H[MSY]")
    return(i)
  })

  ref <- data.frame(qname='H/H[MSY]', data=1)

  plotListFLQuants(ms) +
    geom_hline(data=ref, aes(yintercept=data), linetype=2, alpha=0.60) +
    scale_y_facet(qname == 'H/H[MSY]', limits = c(0, 5))
}

# }}}
