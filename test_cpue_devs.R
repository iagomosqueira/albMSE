# test_cpue_devs.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/test_cpue_devs.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD input(s) from v2
load('../abc_tuna/v2/data.RData')
ori <- FLQuant(c(inp$I[,,1]), dimnames=list(age='all', year=2000:2020, season=1:4))

# LOAD SS3
load('data/base.rda')
nwi <- Reduce(join, lapply(base$ids[1:4], index))[, ac(2000:2020)]

# TEST: COMPARE v2 input and SS3
range(nwi / ori)

# LOAD MC output
load('data/om5b/mcout_abc5b.rda')

# GET estimated index
iha <- out$index.hat

# RE-COMPUTE index
ies <- quantSums(unitSums(out$stock.n * out$catch.sel[,,,,1] *
  wt(base$bio)[, ac(2000:2020)]))

# TEST: COMPARE
range(iha / ies)

# PLOT

plot(FLQuants(EST=iha, OBS=nwi))
plot(FLQuants(EST=iter(iha, 1), OBS=nwi))

# NOTE: HOW to generate index in original scale?


# -- COMPUTE deviances as in mcmc

res <- log(nwi[,ac(2000:2020)] / out$index.hat)

lnr <- mean(res)

des <- res - lnr

plot(exp(des))

# COMPUTE meanlog, sdlog and rho
mea <- yearMeans(des)
sda <- sqrt(yearVars(des))
rha <- rho(des)

plot(FLPar(MU=mea, SD=sda, RHO=rha))

# GENERATE deviances

rlnormar1(n=500, meanlog=mea, sdlog=sda, rho=rha, years=2010:2023)

plot(rlnormar1(n=500, meanlog=mea, sdlog=sda, rho=rha, years=2010:2023))
