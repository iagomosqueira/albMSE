# model_buffer.R - DESC
# abc_tuna/v3/model_buffer.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: TEST mp() runs
# TODO: SKELETON model_.R
# TODO: args(ctrl, 'hcr'), method(ctrl, 'hcr')

source("config.R")

# LOAD pcbar (catch proportions by fleet & season)
pcbar <- as.matrix(fread('data/pcbar.dat'))

# LOAD ALK [len, age, season, sex]
load('data/pla.rda')

# LOAD om
# load('data/om5b.rda')
load('data/om5b_updated.rda')

# ASSIGN args to om@projection
args(projection(om)) <- list(pla=pla, pcbar=pcbar)
method(projection(om)) <- fwdabc.om

om <- iter(om, seq(50))
oem <- iter(oem, seq(50))

# SETUP

iy <- 2023
fy <- 2045

# PLAN
plan(multicore, workers=5)

# --- TEST 0 mp(shortcut.sa + fixedC.hcr) {{{

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=FLQuant(runif(fy-iy+1, min=15000, max=25000),
      dimnames=list(year=seq(iy + 1, fy)))))
))

tes0 <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=1))

catch(tes0)

performance(tes0, statistics=statistics, metrics=mets,
  om="abc5b", type="test", run="fixedC")

performance(om, statistics=statistics['green'], metrics=mets)[year == 2023, mean(data)]

# TODO: args(ctrl, 'hcr'), method(ctrl, 'hcr')

# - DO stock dynamics make sense?
ssb(tes0)
ssb(tes0) / sb0(tes0)

unitSums(rec(tes0)[[1]][,,,4])

# }}}

# --- TUNE shortcut.sa + fixedC {{{

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=25000))
))

system.time(
  tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2026, frq=1))
)

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(ctrg=c(15000, 60000)), prob=0.6, tol=0.05, maxit=12)
)

performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="tune", run="fixedC")

performance(tune)[statistic == 'green' & year %in% seq(2034, 2038), mean(data)]

# C = 43125
args(control(tune)$hcr)$ctrg

# }}}

# TODO:
# [X] RUN hindcasst with new NC
# [X] TUNE fixedC.hcr to 60% kobe
# [X] catch up or down? DOWN to 43000 t
# [ ] COMPARE idx ~ kobe(sa)
# [ ] select refyrs
# [ ] TUNE refyrs, fixed width
# [ ] EXPLORE width vs. upp/low limits

# --- TUNE cpuescore.ind + bufferdelta.hcr(zscore) {{{

# - EXPLORE idx ~ SA

library(ss3om)

sso <- readOutputss3('boot/data/base')

plot(FLQuants(NW=seasonMeans(index(observations(oem)$ALB$idx$NW))[, ac(2000:2020)],
  SB0=extractSSB(sso)[, ac(2000:2020)] / sso$SBzero)) +
  ylim(0, NA)

#
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore.ind,
    args=list(index=1, refyrs=c(2000:2005, 2015:2020))),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=0, width=1.5, sloperatio=0.15, metric="zscore",
      initac=42000))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2029, frq=3))#, .DEBUG=TRUE)

# TODO: plot_buffer.hcr(results=TRUE)

plot(om, tes)

pe <- performance(tes, statistics=statistics, metrics=mets,
  om="abc5b", type="test", run="bufferdelta")

pe[statistic == 'green' & year %in% 2026:2032, mean(data)]

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=1), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(width=c(0.01, 4)), prob=0.6, tol=0.05, maxit=12)
)

performance(tune, statistics=statistics['green'], metrics=mets,
  om="abc5b", type="test", run="bufferdelta")[year %in% 2026:2032, mean(data)]

args(control(tune, 'hcr'))

plot(om, tune)

#
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=25000, bufflow=0.30, buffupp=0.60, sloperatio=0.15,
      metric=function(x) ssb(x)[,,,4] / refpts(om)$SB0[,4], dupp=1.15, dlow=0.85, 
      initac=42000))
))

tes <- mp(om, oem=oem, control=ctrl, args=list(iy=iy, fy=2029, frq=1))

plot(om, tes)

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=1), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(target=c(15000, 50000)), prob=0.6, tol=0.05, maxit=12, .DEBUG=TRUE)
)



# TODO: PERFORMANCE om
performance(om, statistics=statistics['green'], metrics=mets)[year %in% 2000:2020, .(data=mean(data)), by=year]

plot(relhr(om)[,ac(2000:2020)]) + 
  geom_hline(yintercept=1) + ylim(0, NA)

plot(ssb(om)[,ac(2000:2020)] / sb0(om)) + 
  geom_hline(yintercept=0.40) + ylim(0, NA)

plot(ssb(om)[,ac(2000:2020)] / sbmsy(om)) + 
  geom_hline(yintercept=0.40) + ylim(0, NA)

# COMPARE tsb ~ ssb
plot(tsb(om)[,ac(2000:2020)], ssb(om)[,ac(2000:2020)])

# PLOT %TSB/age

# }}}

# TEST mp(cpuescore.ind + buffer.hcr(zscore)) {{{

plot(zscore(index(observations(oem)$ALB$idx[[1]])[, ac(2000:2020)]))

# REF years? From FULL index?
# GET full dataset for JABBA

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore.ind, args=list(index=1)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=25000, bufflow=-1, buffhigh=1,
      sloperatio=0.15, metric="zscore", initac=40000))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=2020, fy=fy, frq=3))

plot(om, tes)

catch(tes)

performance(tes1, statistics=statistics, metrics=mets,
  om="abc5b", type="test", run="bufferdelta")

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=2020, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=seq(2026, fy),
  tune=list(target=c(15000, 50000)), prob=0.6, tol=0.05, maxit=12)
)


plot(om, tune05)

# }}}

plot(
FLQuants(SSB=ssb(om(tes)), R=rec(om(tes))[[1]][,,,4], HR=relhr(om(tes)))
)
