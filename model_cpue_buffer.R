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
load('data/om5b/om5b_updated.rda')

# EXTRACT
om <- iter(om, seq(100))
oem <- iter(oem, seq(100))

# STORE in db
writePerformance(performance(om, statistics[c("SB", "SB0", "SBMSY", "HRMSY", "C")],
  metrics=mets, years=2000:2024), overwrite=TRUE)

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP

iy <- 2024
fy <- 2042
ty <- seq(iy + 11, iy + 15)

# PLAN
plan(multicore, workers=5)

# RESULTS
res <- list()

# PROJECTIONS {{{

load('data/alb_2025_nw.rda')

performance(window(om, end=2024), statistics=statistics['green'],
  metrics=mets)[, .(Pgreen=mean(data)), by=year]

# FWD(C=C0) 

ctrl <- fwdControl(year=2025:2045, biol=1, quant='catch', value=0)

system.time(
fom_c0 <- fwdabc.om(om, ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# FWD(C=C2024) 

ctrl <- fwdControl(year=2025:2045, biol=1, quant='catch', value=catch(om)[,'2024'][[1]])

system.time(
fom_c2024 <- fwdabc.om(om, ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# FWD(C=MSY2025) 

ctrl <- fwdControl(year=2024:2045, biol=1, quant='catch', value=ss25$rps$MSY)

system.time(
fom_cmsy <- fwdabc.om(om, ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# SAVE
fom <- list(C2024=fom_c2024, C0=fom_c0, MSY=fom_cmsy)

save(fom, file="model/om5b_fwd.rda", compress="xz")

# }}}

# --- TUNE shortcut.sa + fixedC.hcr {{{

# SET control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=catch(om)[,'2024'][[1]]))
))

# TEST <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2026, frq=3), .DEBUG=FALSE)

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(ctrg=c(15000, 60000)), prob=0.6, tol=0.01, maxit=12)
)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="fixedC", run="tune_kobe60")

# WRITE to table
writePerformance(performance(tune))

# STORE in results
save(tune, file="model/runs/om5b_fixedC_tune_kobe60.rda", compress="xz")

# CHECK Kobe green
performance(tune)[statistic == 'green' & year %in% ty, mean(data)]
performance(tune)[statistic == 'green', mean(data), by=year]

# GET tuned C = 40 664 t
args(control(tune)$hcr)$ctrg

# PLOT
plotMetrics(OM=window(om, end=2023), CtrgK60=window(om(tune), start=2023)) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# }}}

# --- EXPLORE bufferdelta.hcr

# --- EXPLORE bufferhcr

# --- TUNE cpuescore.ind + buffer.hcr(C~zscore) {{{

zscore <- function(x, mean=yearMeans(x), sd=sqrt(yearVars(x)))
  exp((x %-% mean) %/% sd)

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpuescore.ind,
    args=list(index=1, refyrs=c(2000:2005, 2015:2020))),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=42070, lim=0.10, bufflow=0.8, buffupp=1.2, sloperatio=0.15,
      dlow=0.85, dupp=1.15, metric="zscore", initac=catch(om)[,'2024'][[1]]))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3))#, .DEBUG=TRUE)

# PLOT

plotMetrics(OM=window(om(tes), end=2023), CtrgK60=window(om(tes), start=2023)) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(target=c(25000, 55000)), prob=0.6, tol=0.01, maxit=12)
)

# PLOT
plotMetrics(OM=window(om, end=2023), K60=window(om(tune), start=2023))

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="zscore_NW_buffer_C", run="tune_kobe60")

performance(tune)[statistic == 'green' & year %in% ty, mean(data)]
performance(tune)[statistic == 'green', mean(data), by=year]

# WRITE to table
writePerformance(performance(tune))

# STORE in results
res[["abc5b_zscore_NW_buffer_C_tune_kobe60"]] <- tune

# GET value of tuned argument
args(control(tune, 'hcr'))$target

# }}}

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
