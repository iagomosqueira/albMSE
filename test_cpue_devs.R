# test_cpue_devs.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/test_cpue_devs.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# -- TEST: COMPUTE residuals Ihat (RH)

source("config.R")

# LOAD input(s) from v2
load('../abc_tuna/v2/data.RData')
ori <- FLQuant(c(inp$I[,,1]), dimnames=list(age='all', year=2000:2020, season=1:4))

# LOAD SS3
load('data/base.rda')
nwi <- Reduce(join, lapply(base$ids[1:4], index))

# LOAD MC output
load('data/om5b/mcout_abc5b.rda')
omi <- out$index.hat %*% out$index.q

# PLOT CPUEs
plot(nwi, omi) +
  theme(legend.position='bottom') +
  scale_fill_manual('CPUE', labels=c('SS3 NW', 'MC output'),
    values=c('black', 'red')) +
  scale_color_manual('CPUE', labels=c('SS3 NW', 'MC output'),
    values=c('black', 'red'))

# RECONSTRUCT index from N * S * wt
Iest <- quantSums(out$stock.n %*% out$catch.sel[,,,,1] * stock.wt(base$stk)[, ac(2000:2020)])

# COMPARE to Ihat from mcoutput
plot(Iest, out$index.hat)
