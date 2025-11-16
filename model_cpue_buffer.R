# model_buffer.R - DESC
# albMSE/model_cpue_buffer.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD 100 iter objects
qs_readm("data/om5b.qs2")

# BUG: EMPTY future index to check values
index(observations(oem)$ALB$idx[[1]])[, ac(2025:2045)] <- NA

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP
iy <- 2024
fy <- 2050
ty <- seq(iy + 11, iy + 15)

# PLAN
plan(multicore, workers=10)

# --- 1. TUNE cpues.ind + buffer.hcr(mult~wmean) {{{

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind, args=list(index=1, nyears=4)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(lim=0.25, bufflow=0.35, buffupp=0.60, sloperatio=0.10,
      dlow=0.85, dupp=1.15, metric="wmean", initac=36458))
))

#  tes <- mp(om, oem, control=ctrl, args=list(iy=2024, fy=2032, frq=3), .DEBUG=FALSE)

# - TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(0.50, 0.75)), prob=0.6, tol=0.01)
)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="buffer-wmean", run="kobe60")

writePerformance(performance(tune))

# SAVE
save(tune, file="model/runs/om5b_buffer-wmean_tune_kobe60.rda", compress="xz")

# STORE summary
appendSummary(tune, mp=unique('om5b_buffer-wmean_tune_kobe60'))


# -- CHECK:

# KOBE performance
performance(tune)[statistic=='green' & year %in% ty, .(PKg=mean(data)), by=mp]
performance(tune)[statistic=='green', .(PKg=mean(data)), by=year]

# PLOT time series
p1 <- plot(om, K60=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# PLOT HCR & future observations TODO: ADD decisions
p2 <- plot_buffer.hcr(control(tune)$hcr) +
  geom_point(data=data.table(met=c(index(observations(oem(tune))$ALB$idx[[1]])), out=0),
    alpha=0.01)

# PLOT observed index
p3 <- plot(seasonMeans(index(observations(oem(tune))$ALB$idx[[1]]))) +
  geom_hline(yintercept=c(0.25, 0.35, 0.72), color=c('red', 'black', 'black'),
    linetype=c(1,2,2)) +
  ylab("CPUE LL1 NW") + ylim(0, 1.5)

p1 + (p2 / p3)

# --:

# }}}

# --- 2. TUNE cpues.ind + buffer.hcr(mult~mean) {{{

# SET control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind, args=list(index=1, nyears=4)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(lim=0.25, bufflow=0.35, buffupp=0.60, sloperatio=0.10,
      dlow=0.85, dupp=1.15, metric="mean", initac=36458))
))

# tes <- mp(om, oem, control=ctrl, args=list(iy=2024, fy=2032, frq=3), .DEBUG=FALSE)

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(0.50, 0.75)), prob=0.6, tol=0.01)
)

# PLOT
plot(om, mean=tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="buffer-wmean", run="kobe60")

# WRITE to performance table
writePerformance(performance(tune))

# SAVE
save(tune, file="model/runs/om5b_buffer-mean_tune_kobe60.rda", compress="xz")

# STORE summary
appendSummary(tune, mp=unique('om5b_buffer-mean_tune_kobe60'))

# -- CHECK:

# KOBE performance
performance(tune)[statistic=='green' & year %in% ty, mean(data), by=mp]

# PLOT time series
plot(om, K60=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# PLOT HCR & future observations TODO: ADD decisions
plot_buffer.hcr(control(tune)$hcr) +
  geom_point(data=data.table(met=c(index(observations(oem(tune))$ALB$idx[[1]])), out=0),
    alpha=0.01)

# PLOT observed index
plot(seasonMeans(index(observations(oem(tune))$ALB$idx[[1]]))) +
  geom_hline(yintercept=c(0.25, 0.35, 0.72), color=c('red', 'black', 'black'),
    linetype=c(1,2,2)) +
  ylab("zscore(CPUE LL1 NW)") + ylim(0, 1.5)

# }}}

# --- ROBUSTNESS runs {{{

# LOAD run
load('model/runs/om5b_buffer-wmean_tune_kobe60.rda')

# or LOAD summary

qs_readm("data/om5a.qs2")

tes1 <- mp(om, oem, control=control(tune), args=list(iy=2024, fy=2050, frq=3))

qs_readm("data/om6b.qs2")

tes2 <- mp(om, oem, control=control(tune), args=list(iy=2024, fy=2050, frq=3))

plot(om, TUNE=tune, T1=tes[[1]], T2=tes[[2]])

tes <- FLmses(list('kobe60-om5a-robust'=tes1, 'kobe60-om6b-robust'=tes2),
  statistics=statistics, metrics=mets, type="buffer-wmean")

save(tes, file="model/runs/om5b_buffer-wmean_kobe60-robust.rda", compress="xz")

writePerformance(performance(tes))

# PLOT from performance

library(mseviz)

perf <- readPerformance()

# BUG: REMOVE very large HRMSY for plots
perf[statistic == "HRMSY", data:=ifelse(data > 3, 3, data)]

# RE-LABEL performance elements
perf <- labelPerformance(perf, labels=list(
  'om5b'='OM base',
  'om6b'='OM 1% Q',
  'om5a'='OM SW',
  'om5a_buffer-wmean_kobe60-om5a-robust'="Robust SW",
  'om5b_buffer-wmean_kobe60'="Tuned Kobe 60%",
  'om6b_buffer-wmean_kobe60-om6b-robust'="Robust 1% LL Q"))

plotTimeSeries(perf)

#
tun <- periodsPerformance(perf, list(tune=2034:2038))

plotBPs(tun, c("C", "SBMSY", "HRMSY", "IACC", "green"), show.mean=c('green'))

plotTOs(tun, x="C", y=c("SBMSY", "HRMSY", "IACC", "green"))

kobeMPs(tun,x = "SBMSY", y = "HRMSY", probs = c(0.25, 0.50, 0.75)) +
  ylab(expression(HR/HR[MSY]))

# }}}
