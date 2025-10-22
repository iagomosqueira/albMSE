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
fy <- 2042
ty <- seq(iy + 11, iy + 15)

# PLAN
plan(multicore, workers=5)

# RESULTS
res <- list()

# PROJECTIONS {{{

load('data/alb_2025_nw.rda')

performance(window(om, end=2023), statistics=statistics['green'],
  metrics=mets)[, .(Pgreen=mean(data)), by=year]

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

# --- TEST: 0 mp(shortcut.sa + fixedC.hcr, frq=1) {{{

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

# DO stock dynamics make sense?

ssb(tes0)
ssb(tes0) / ssb0(tes0)

unitSums(rec(tes0)[[1]][,,,4])

# }}}

# --- TUNE shortcut.sa + fixedC.hcr {{{

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

# EXPLORE ctrg
exp <- mps(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2026, frq=3),
  hcr=list(ctrg=seq(20000, 50000, length=5)))

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(ctrg=c(15000, 60000)), prob=0.6, tol=0.01, maxit=12)
)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="fixedC", run="tune_kobe60")

# WRITE to table
writePerformance(performance(tune))

# STORE in results
res[["om5b_fixedC_tune_kobe60"]] <- tune

# CHECK Kobe green
performance(tune)[statistic == 'green' & year %in% ty, mean(data)]

# C = 41719
args(control(tune)$hcr)$ctrg

# PLOT
plotMetrics(OM=window(om, end=2023), CtrgK60=window(om(tune), start=2023)) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# }}}

# --- TUNE cpuescore.ind + buffer.hcr(C~zscore) {{{

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
plotMetrics(OM=window(om, end=2023), K60=window(om(tune), start=2023))

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="zscore_NW_buffer_C", run="tune_kobe60")

# WRITE to table
writePerformance(performance(tune))

# STORE in results
res[["abc5b_zscore_NW_buffer_C_tune_kobe60"]] <- tune

# GET value of tuned argument
args(control(tune, 'hcr'))$target

# }}}

# TODO: ADD om perfomance to table

# SAVE
save(res, file="model_cpue_buffer.rda", compress="xz")

# --- DOES NOT TUNE -- TUNE cpuescore.ind + bufferdelta.hcr(zscore) {{{

# FIX width and tune for buffupp
# EXPLORE buffupp and bufflow

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore.ind,
    args=list(index=1, refyrs=c(2000:2005, 2015:2020))),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=0, width=1, buffupp=1.5, sloperatio=0.15, dlow=0.85, dupp=1.15,
      metric="zscore", initac=42000))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3))#, .DEBUG=TRUE)

exp <- mps(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3),
  hcr=list(buffupp=seq(-0.5, 3, length=5)))

performance(tune[[1]], statistics=statistics['green'], metrics=mets)[year %in% ty, mean(data)]

# PLOT
plotMetrics(OM=window(om, end=2023),
  TES=window(om(tes), start=2023))

# TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(1, 4)), prob=0.6, tol=0.01, maxit=12)
)

# TEST:
plotMetrics(OM=window(om, end=2023),
  A=window(om(tune[[1]]), start=2023),
  B=window(om(tune[[2]]), start=2023)
)

performance(tune[[1]], statistics=statistics['green'], metrics=mets)[year %in% ty, mean(data)]

performance(tune[[2]], statistics=statistics['green'], metrics=mets)[year %in% ty, mean(data)]

# PLOT
plotMetrics(OM=window(om, end=2023),
  K60=window(om(tune[[1]]), start=2023))

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="tune", run="kobe60")

# CHECK P(Kobe=green) in ty

performance(tune)[statistic == 'green' & year %in% ty, mean(data)]

# }}}
