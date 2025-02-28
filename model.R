# model.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albMSE/model.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

source('config.R')

# LOAD
load('boot/initial/data/pcbar_2017-2020.rda')
load('data/om6b.rda')

nit <- 100
its <- seq(nit)

om <- iter(om, its)
oem <- iter(oem, its)

df <- 6

# SET harvest attr TODO: SET as slot in FLombf
attr(om, 'harvest') <- FLQuants(ALB=expand(n(biol(om)), area=1:df) %=%
  as.numeric(NA))
units(attr(om, 'harvest')[[1]]) <- 'hr'

om@refpts[[1]]$SBlim <- om@refpts[[1]]$B0 * 0.10

method(projection(om)) <- fwdabc.om
args(projection(om)) <- list(split=pcbar)

method(oem) <- shortcut.oem

# TEST fwdabc.om

ctrl <- propagate(fwdControl(year=2021:2045, quant="catch",
  value=runif(25, 47000, 53000)), nit)

ctrl <- propagate(fwdControl(year=2021:2045, quant="catch",
  value=35000), nit)

tes <- fwdabc.om(om, ctrl, split=pcbar)$om

# - CHECK targets

unitSums(seasonSums(catch(tes)[[1]]))[, ac(2021:2045)]
unitSums(sesonSums(landings(tes)[[1]]))[, ac(2021:2045)]

ctrl

seasonSums(unitSums(Reduce('+', landings(fisheries(tes)))))

Reduce('+', lapply(seq(df), \(x)
  seasonSums(unitSums(Reduce('+', landings(fisheries(tes)[x]))))))

# BUG:
all.equal(
landings(fisheries(tes)),
landings(tes)
  )

it <- 1
it <- 5

Reduce('+', lapply(seq(df), \(x)
  seasonSums(unitSums(quantSums(
    iter(landings.n(fisheries(tes)[[x]][[1]]), it) *
    iter(landings.wt(fisheries(tes)[[x]][[1]]), it))))[, ac(2021:2045)]
))

#

taf.png("om.png")
plot(biol(om))
dev.off()

taf.png("om-annual.png")
plot(window(FLQuants(
  R=unitSums(seasonSums(rec(biol(om)))),
  SB=unitSums(seasonSums(ssb(biol(om)))),
  C=unitSums(seasonSums(catch(om)[[1]]))), end=2020)) +
  ylim(0, NA)
dev.off()

taf.png("om-refpts.png", height=700)
plot(refpts(om)[c(1,2,3), ])
dev.off()




taf.png("biol_tes-C35k.png")
plot(biol(tes)) +
  geom_vline(xintercept=ISOdate(2022,1,1), alpha=0.5, linetype=2)
dev.off()

taf.png("summ_tes-C35k.png")
plot(FLQuants(
  R=unitSums(seasonSums(rec(biol(tes)))),
  SB=unitSums(seasonSums(ssb(biol(tes)))),
  C=unitSums(seasonSums(catch(tes)[[1]])))) +
  ylim(0, NA) +
  geom_vline(xintercept=2021)
dev.off()

# - CHECK proportions

# RUN mp


tdepletion <- function(x, B0) unitSums(ssb(x))[,,,4] / B0

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa, args=list(metric="tdepletion",
    B0=refpts(om)$B0)),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr, args=list(ctrg=45000))
))

# BUG:
bug <- mp(om, oem, ctrl=ctrl, args=list(iy=2020, fy=2032, frq=1))

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa, args=list(metric="catch")),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr, args=list(ctrg=45000))
))


bug <- mp(om, oem, ctrl=ctrl, args=list(iy=2020, fy=2026, frq=3))

plot(om, C30k=bug)

attr(om(bug), 'harvest')[[1]][, ac(2021:2030)]

#
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa, args=list(metric="tdepletion",
    B0=refpts(om)$B0)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=40000, bufflow=0.30, buffupp=0.50, lim=0.10, min=0,
      sloperatio=0.20, metric="tdepletion"))
))

system.time(
bug <- mp(om, oem, ctrl=ctrl, args=list(iy=2020, fy=2035, frq=3))
)

# BUG:

tune05 <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=2021, frq=3), 
  statistic=statistics['PSBlim'], years=2030:2040, metrics=mets,
  tune=list(target=c(2000, 100000)), prob=0.1, tol=0.01, maxit=12)

performance(tune05[[2]], statistics=statistics['PSBlim'], metrics=mets,
  years=2030:2040)

plot(om, tune05)
