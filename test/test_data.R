# test_data.R - TEST hr, catch and catch match
# abc_tuna/v3/test.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

source("utilities.R")

load('data/om5b.rda')

# --- TEST hr and catch.n {{{

# COMPUTE hr_fys
hom <- hrf(om)

dim(hom[[1]]) # [1]   1  21   1   4   1 500

# COMPUTE catch.n by age, year, unit, season, iter
can <- lapply(hom, function(x)
  expand(x, age=0:14, unit=c('F', 'M')) * n(biol(om)))

# MULTIPLY catch.n by catch.wt and selex by fishery
ca <- Reduce('+',Map(function(ca, fi) quantSums(seasonSums(unitSums(ca *
  landings.wt(fi[[1]]) * catch.sel(fi[[1]])))), ca=can, fi=fisheries(om)))

# and calculated from catch.n * catch.wt in om
# INFO: ADD over ages(sexes(seasons(elements of catch by fishery list
quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(om))))))

# DIRECT calculation on F1

# TODO: fishery(om, 1)
print(seasonSums(unitSums(catch(fisheries(om)[[1]][[1]]))))

# H_f = C_f / sum(S_f * N * W)
hr1 <- unitSums(catch(fisheries(om)[[1]][[1]])) /
  quantSums(unitSums(catch.sel(fisheries(om)[[1]][[1]]) * n(biol(om)) * wt(biol(om))))

# COMPUTE catch.n by age, year, unit, season, iter
ca1 <- expand(hr1, age=0:14, unit=c('F', 'M')) * n(biol(om))

# COMPARE computed from HR, ...
quantSums(unitSums(seasonSums(ca1 * wt(biol(om)))))

# and calculated from catch.n * catch.wt in om
unitSums(seasonSums(catch(fisheries(om)[[1]][[1]])))

# }}}

# --- TEST fwdabc.om {{{

library(Rcpp)
library(parallel)
library(mvtnorm)

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# LOAD pcbar (catch proportions by fleet & season)
pcbar <- as.matrix(fread('data/pcbar.dat'))

# LOAD ALK [len, age, season, sex]
load('data/pla.rda')

# LOAD om
load('data/om5b.rda')

its <- seq(50)
iy <- 2017

# - TEST fwd iy:2020 with historical catch

ctrl <- fwdControl(year=iy:2020, quant="catch",
  value=c(unitSums(seasonSums(Reduce('+', catch(om))))[, ac(iy:2020),,,,1]))

system.time(
tes1 <- fwdabc.om(iter(om, its), ctrl, pcbar=pcbar, pla=pla)
)

# COMPARE catch
quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(tes1$om))))))
quantSums(unitSums(seasonSums(Reduce('+', iter(catch(fisheries(om)), 1)))))

# GET hr
plot(Reduce('+', unitMeans(hrf(tes1$om))),

# - TEST larger catch

ctrl <- fwdControl(year=iy:2020, quant="catch",
  value=1.5 * c(unitSums(seasonSums(Reduce('+', catch(om))))[, ac(iy:2020),,,,1]))

tes2 <- fwdabc.om(iter(om, its), ctrl, pcbar=pcbar, pla=pla)

# COMPARE catch
quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(tes2$om))))))
quantSums(unitSums(seasonSums(Reduce('+', iter(catch(fisheries(om)), 1)))))


  Reduce('+', unitMeans(hrf(om))))



# }}}
