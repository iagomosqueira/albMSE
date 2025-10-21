# model_buffer.R - DESC
# abc_tuna/v3/model_buffer.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: TEST mp() runs
# TODO: SKELETON model_.R
# TODO: args(ctrl, 'hcr'), method(ctrl, 'hcr')
# TODO: plot_buffer.hcr(results=TRUE)

# TODO:
# [X] RUN hindcasst with new NC
# [X] TUNE fixedC.hcr to 60% kobe
# [X] catch up or down? DOWN to 43000 t
# [ ] COMPARE idx ~ kobe(sa)
# [ ] select refyrs
# [ ] TUNE refyrs, fixed width
# [ ] EXPLORE width vs. upp/low limits


source("config.R")

# LOAD om
load('data/om5b_updated.rda')

# EXTRACT
om <- iter(om5b$om, seq(100))
oem <- iter(om5b$oem, seq(100))

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP

iy <- 2023
fy <- 2045
ty <- seq(iy + 11, iy + 15)

# PLAN
plan(multicore, workers=5)

# PROJECTIONS {{{

# FWD(C=C0) 

ctrl <- fwdControl(year=2024:2045, biol=1, quant='catch', value=0)

system.time(
fom_c0 <- fwdabc.om(om, ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# FWD(C=MSY2025) 

ctrl <- fwdControl(year=2024:2045, biol=1, quant='catch', value=ss25$rps$MSY)

system.time(
fom_cmsy <- fwdabc.om(om, ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# SAVE
fom <- list(C2023=fom_c2023, C0=fom_c0, MSY=fom_cmsy)

save(fom, file="model/om5b_fwd.rda", compress="xz")

# }}}

# --- TEST 0 mp(shortcut.sa + fixedC.hcr, frq=1) {{{

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=FLQuant(runif(fy-iy+1, min=15000, max=25000),
      dimnames=list(year=seq(iy + 1, fy)))))
))

tes0 <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=1))

# COMPARE catch

catch(tes0)
args(ctrl$hcr)

# COMPUTE performance
performance(tes0, statistics=statistics, metrics=mets,
  om="abc5b", type="test", run="fixedC")

# - DO stock dynamics make sense?

ssb(tes0)
ssb(tes0) / ssb0(tes0)

unitSums(rec(tes0)[[1]][,,,4])

# }}}

# --- TUNE shortcut.sa + fixedC {{{

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=25000))
))

# RUN
system.time(
  tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2026, frq=3))
)

# TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(ctrg=c(15000, 60000)), prob=0.6, tol=0.05, maxit=12)
)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="tune", run="fixedC")

# CHECK Kobe green
performance(tune)[statistic == 'green' & year %in% ty, mean(data)]

# C = 45938
args(control(tune)$hcr)$ctrg

# }}}

# --- TUNE cpuescore.ind + buffer.hcr(zscore) {{{

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore.ind,
    args=list(index=1, refyrs=c(2000:2005, 2015:2020))),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=25000, lim=-2, bufflow=-1, buffupp=1, sloperatio=0.15,
      dlow=0.85, dupp=1.15, metric="zscore", initac=42000))
))

# RUN
tes <- mp(iter(om, seq(5)), iter(oem, seq(5)), ctrl=ctrl,
  args=list(iy=iy, fy=fy, frq=3))#, .DEBUG=TRUE)

# PLOT
plotMetrics(OM=window(iter(om, seq(5)), end=2023),
  TES=window(om(tes), start=2023))

# TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(target=c(35000, 55000)), prob=0.6, tol=0.01, maxit=12)
)

# PLOT
plotMetrics(OM=window(om, end=2023),
  K60=window(om(tune), start=2023))

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="tune", run="kobe60")

# CHECK P(Kobe=green) in ty

performance(tune)[statistic == 'green' & year %in% ty, mean(data)]

performance(tune)[statistic == 'SBMSY' & year %in% ty, mean(data)]
performance(tune)[statistic == 'HRMSY' & year %in% ty, mean(data)]

plot(mets$HR(om(tune)))

mets$HR(om(tune))[, ac(ty)]

# GET value of tuned argument
args(control(tune, 'hcr'))$target

plot(window(catch(om), end=2023), window(catch(tune), start=2023)) +
  geom_hline(yintercept=args(control(tune, 'hcr'))$target)

# }}}

# ----- STOP HERE -----

# --- TUNE cpuescore.ind + bufferdelta.hcr(zscore) {{{

# - EXPLORE idx ~ SA

load('data/base.rda')

bae$

library(ss3om)

sso <- readOutputss3('boot/data/base')

plot(FLQuants(NW=seasonMeans(index(observations(oem)$ALB$idx$NW))[, ac(2000:2020)],
  SB0=ssb(base$stk)[, ac(2000:2020),'F',1] / base$rps$SB0)) +
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
tes <- mp(iter(om, seq(5)), iter(oem, seq(5)), ctrl=ctrl,
  args=list(iy=iy, fy=2028, frq=1))#, .DEBUG=TRUE)

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
