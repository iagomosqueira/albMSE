# utilities.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albMSE/utilities.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(Rcpp)

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# fwdabc.om {{{

fwdabc.om <- function(om, ctrl, pcbar, pla, ...) {

  # DIMS
  dom <- dims(om)
  
  yrs <- ac(unique(ctrl$year))
  yr <- yrs[1]
  dyr <- ac(min(ctrl$year) - 1)
  nyrs <- length(yrs)
  
  na <- dims(biol(om))$age
  nf <- length(fisheries(om))
  its <- seq(dims(om)$iter)

  # PROPAGATE ctrl
  if(dim(iters(ctrl))[3] == 1 & dims(om)$iter > 1) {
    ctrl <- propagate(ctrl, dims(om)$iter)
  }

  n(biol(om))[, yrs] <- as.numeric(NA)

  for(it in its) {

  # ASSEMBLE inputs
  inp <-   list(
    # dms: dimensions: year, season, age, lengths, fishery
    dm_=c(nyrs + 1, 4, 15, 27, 6),
    # srec: rec season
    srec_=4,
    # - R0: SRR R0 in numbers
    R0_=c(params(sr(biol(om)))$R0)[it] * 1000,
    # - hh: SRR H
    hh_=c(params(sr(biol(om)))$s)[it],
    # - psi: sex ratio at birth
    psi_=0.5,
    # - epsrx: future rec devs [year]
    epsr_=log(c(deviances(om)[, yrs,,,,it])),
    # - spr0: SRR B0/R0
    spr0_=c(params(sr(biol(om)))$v[,it] / (params(sr(biol(om)))$R0[,it] * 1000)),
    # - M: M
    M_=m(biol(om))[,,,,,it][[1]],
    # - mata: mat at age [age, season, unit (sex)]
    mata_=c(unname(aperm(mat(biol(om))[, yr,,, 1, it, drop=TRUE], c(1, 3, 2)))),
    # - wta: Wt at age [age, season, unit (sex)]
    wta_=c(unname(aperm(wt(biol(om))[, yr,,, 1, it, drop=TRUE],
      c(1, 3, 2)))) / 1000,
    # - sela: Selex at age [age, season, unit (sex), fishery]
    sela_=c(unname(aperm(abind(lapply(fisheries(om), function(x)
      catch.sel(x[[1]])[, yr,,,,it]))[drop=TRUE], c(1,3,2,4)))),
    # - nvec: Last year Ns [age, season, unit (sex)], as vector.
    Ninit_=c(aperm(n(biol(om))[, dyr,,,,it, drop=TRUE], c(1,3,2))) * 1000,
    # - cvec: Catch in projection [year, season, fishery], as vector.
    Cb_=c(
      array(rep(pcbar, each=nyrs+1), dim=c(nyrs+1, 4, 6)) * 
      array(c(sum(catch(om)[[1]][,dyr]), iters(ctrl)[, 2, it]),,
        dim=c(nyrs+1, 4, 6))
      ),
    # - pla: ALK [lengths, age, season, unit (sex)]
    pla_=c(pla),
    # - fcpue: Index of fleet to generate CPUE
    fref_=1)
  
  # CALL pdynlfcpue
  rei <- do.call(pdynlfcpue, inp)

  # EXTRACT n - N [y,a,s,u]
  n(biols(om)[[1]])[, yrs,,,, it] <- 
    aperm(array(rei$N, dim=c(nyrs+1, na, 4, 2)), c(2,1,4,3))[,-1,,] / 1000

  # EXTRACT harvest - H [y,s,f]
  attr(om, 'harvest')[[1]][, yrs,,,, it]  <- (expand(FLQuant(rei$H,
    dimnames=list(year=c(dyr, yrs), season=1:4, area=1:6)),
    unit=c('F', 'M')) %*% abind(lapply(fisheries(om),
    function(x) catch.sel(x[[1]])[, c(dyr, yrs),,,,it])))[,-1]

  # SUM across fisheries
  attr(om, 'hr')[[1]][, yrs,,,, it]  <- 
    FLQuant(rei$H, dimnames=list(year=c(dyr, yrs), season=1:4, area=1:6))[,-1]

  # COMPUTE catch.
  can <- attr(om, 'harvest')[[1]][, yrs,,,, it] %*%
    (n(biols(om)[[1]])[, yrs,,,, it])

  # ASSIGN as landings.n per fleet
  for(f in seq(6))
    fisheries(om)[[f]][[1]]@landings.n[, yrs,,,, it] <- can[,,,,f]
    
  # PRINT total catch per year
  #print(rowSums(array(rep(pcbar, each=nyrs), dim=c(nyrs, 4, 6)) * 
  #  array(iters(ctrl)[, 2, it], dim=c(nyrs, 4, 6))))
  #print(c(Reduce('+', lapply(fisheries(om), function(x)
  #  unitSums(seasonSums(landings(x[[1]])[,yrs,,,,it]))))))
  #print('# --')
  }
  return(list(om=om))
}

# }}}
