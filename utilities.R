# utilities.R - DESC
# albMSE/utilities.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(Rcpp)

# SOURCE Rcpp
sourceCpp("utilities/project.cpp")

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

# }}}

# fwdabc.om {{{

fwdabc.om <- function(om, ctrl, split, ...) {

  args <- list(...)

  # DIMS
  dmb <- dim(biol(om))
  na <- dmb[1]
  nu <- dmb[3]
  ns <- dmb[4]
  ni <- dmb[6]
  nf <- length(fisheries(om))
  
  dy <- ac(min(ctrl$year) - 1)
  cyrs <- ac(ctrl$year)
  nyrs <- length(cyrs) + 1

  # dm_ [fy,s,a,f,i], [10,4,15,6,500]
  dm_ <- c(nyrs, ns, na, nf, ni)
  # srec_ rec season
  srec_ <- 4
  # R0_ McMC R0s
  R0_ <- c(params(sr(biol(om)))$R0)
  # hh_ McMC hs
  hh_ <- c(params(sr(biol(om)))$s)
  # psi_3
  psi_ <- 0.5
  # sigmar_ McMC sigmaR TODO: GET from McMC
  sigmar_ <- runif(ni, 0.22, 0.32)
  # spr0_ McMC SPR0 TODO: GET from McMC
  spr0_ <- rep(0.006769, ni)
  # M_ McMc M
  M_ <- c(m(biol(om))[1,1,1,1,1,])
  # mata_ Matage [a,s,u,i]
  mata_ <- aperm(mat(biol(om))[, cyrs[1]], c(1, 4, 3, 6, 2, 5))[,,,1,1,1,drop=TRUE]
  # wta_ Wtatage [a,s,u]
  wta_ <- aperm(wt(biol(om))[, cyrs[1]], c(1, 4, 3, 6, 2, 5))[,,,1,1,1,drop=TRUE]
  # sela_ McMC selatage [a,s,u,f,i]
  sela_ <- aperm(Reduce(abind, lapply(fisheries(om), function(x)
    catch.sel(x[[1]])[,cyrs[1]])), c(1,4,3,5,6,2))[,,,,,1,drop=TRUE]
  # Ninit_ Init N [a,s,u,i]
  Ninit_ <- aperm(n(biol(om))[, dy], c(1,4,3,6,2,5))[,,,,1,1,drop=TRUE]
  # Cb_ Target C [fy,s,f,i]
  carr <- aperm(array(c(c(apply(catch(om)[[1]][, dy], 6, sum)), ctrl$value),
    dim=c(nyrs, ni, ns, nf)), c(1, 3, 4, 2))
  # carr <- aperm(array(ctrl$value, dim=c(nyrs, ni, ns, nf)), c(1, 3, 4, 2))
  parr <- aperm(array(split, dim=c(ns, nf, nyrs, ni)), c(3, 1, 2, 4))
  Cb_ <- c(carr * parr)
  # q_ McMC Index Q [s,i] TODO: GET from oem
  q_ <- array(3e-06, dim=c(ns, ni))
  # fref_ CPUE fishery number
  fref_ <- 3
  
  # CALL Rcpp
  run <- projv2(dm_,  srec_,  R0_,  hh_, psi_, sigmar_, spr0_, M_, mata_,
    wta_, sela_, Ninit_, Cb_, q_, fref_) 

  # ASSIGN output (n, harvest, ssb and index)
  nhat <- array(run$N, dim=c(nyrs, na, ns, nu, ni))
  hhat <- array(run$H, dim=c(nyrs, ns, nf, ni))
  shat <- array(run$S, dim=c(nyrs, ns, ni))
  ihat <- array(run$I, dim=c(nyrs, ns, ni))
  
  # TEST: COMPUTE catches on year 5
  # lapply over its & fleets
  print(unlist(lapply(seq(dm_[5]), \(i) sum(unlist(lapply(seq(6), \(x)
    # N * S * W * HR
    sum(nhat[5,,,,i] * sela_[,,,x,i] * wta_ *
      array(rep(hhat[5,, x, i], each=15), dim=c(15, 4, 2)))))))))

  # ASSIGN nhat[y,a,s,u,i] to n(biol)[a,y,u,s,i]
  n(biol(om))[, cyrs] <- aperm(nhat[-1,,,,,drop=FALSE], c(2,1,4,3,5))

  # CREATE hr by F  TODO: LOAD historical
  hrf <- FLQuant(c(hhat[-1,,,,drop=FALSE]), dimnames=list(age='all', year=cyrs,
    season=seq(ns), area=seq(nf), iter=seq(ni)), units='hr')

  # ASSIGN to harvest(om)
  #harvest(om)[, cyrs] <- Reduce(abind, lapply(seq(nf), function(x)
  attr(om, 'harvest')[[1]][, cyrs] <- Reduce(abind, lapply(seq(nf), function(x)
    catch.sel(fisheries(om)[[x]][[1]])[, cyrs,] %*% expand(hrf[,,,,x], 
    unit=c('F', 'M'))))

  # ASSIGN landings.n
  for(i in seq(nf))
    landings.n(fisheries(om)[[i]][[1]])[, cyrs] <-
      attr(om, 'harvest')[[1]][, cyrs,,,i] * n(biol(om))[,cyrs] 
#       * exp(-m(biol(om))[,cyrs] / 2)

 return(list(om=om))
}

# }}}

# hrbar
hrbar <- function(x) 
  quantMeans(areaSums(unitSums(seasonSums(attr(x, 'harvest')[[1]]))))
