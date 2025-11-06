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

# [ ] COMPARE idx ~ kobe(sa)
# [ ] select refyrs
# [ ] TUNE refyrs, fixed width
# [ ] EXPLORE width vs. upp/low limits

source("config.R")

# LOAD 100 iter objects
qs_readm("model/om5b.qs2")

# STORE in db
writePerformance(performance(om, statistics[c("SB", "SB0", "SBMSY", "R", "HRMSY", "C")],
  metrics=mets, years=2000:2024), overwrite=TRUE)

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP
iy <- 2024
fy <- 2042
ty <- seq(iy + 11, iy + 15)

# PLAN
plan(multicore, workers=5)

# --- TUNE shortcut.sa + fixedC.hcr {{{

# SET control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=catch(om)[,'2024'][[1]]))
))

# test <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2026, frq=3), .DEBUG=TRUE)

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(ctrg=c(35000, 50000)), prob=0.6, tol=0.01, maxit=12)
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

# GET tuned C = 40 625 t
args(control(tune)$hcr)$ctrg

# PLOT
plot(om, "Constant catch Kobe 60%"=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# }}}

# --- DOES NOT TUNE -- TUNE cpuescore.ind + bufferdelta.hcr(zscore) {{{

load('data/base.rda')

# INDEX
ind <- index(observations(oem)$ALB$idx[[1]])[, ac(2000:2020)]

# COMPUTE SBMSY 2000:2020
sbmsy <- ssb(base$stk)[, ac(2000:2020),1,1] / base$rps$SBMSY

# FIND years where 1 < SBMSY < 2
ref <- seasonMeans(ind[, which(sbmsy > 1 & sbmsy < 2.5)])

yearMeans(ref)

exp(zscore(seasonMeans(ind), mean=yearMeans(ref), sd=sqrt(yearVars(ref))))

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind,
    args=list(index=1, mean=yearMeans(ref), sd=sqrt(yearVars(ref)), nyears=4)),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=0.80, bufflow=0.5, buffupp=1.5, sloperatio=0.15,
      metric="wmean", initac=42000))))

# - TEST run
test <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3), .DEBUG=FALSE)

# PLOT
plotMetrics(OM=window(om, end=2024), TEST=window(om(test), start=2024))

# KOBE performance
performance(test, statistics=statistics['green'], metrics=mets)[year %in% ty, mean(data)]

# - TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(width=c(0.25, 0.45)), prob=0.6, tol=0.05, maxit=16)
)

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(target=c(0.20, 1.5)), prob=0.6, tol=0.05, maxit=16)
)


# PLOT

plotTimeSeries(readPerformance())

plotMetrics(OM=window(om, end=2024), TEST=window(om(tune), start=2024))
plotMetrics(OM=window(om, end=2024), T=window(om(tune[[1]]), start=2024),
  T2=window(om(tune[[2]]), start=2024))

# KOBE performance
performance(tune[[1]], statistics=statistics['green'], metrics=mets)[year %in% ty, mean(data)]

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  om="abc5b", type="tune", run="kobe60")

# }}}

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

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(target=c(25000, 55000)), prob=0.6, tol=0.01, maxit=12)
)

# PLOT
plot(om, CCK60=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

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
