# model_buffer.R - DESC
# albMSE/model_cpue_buffer.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD 100 iter objects list(om, oem)
om5b <- readRDS("data/om5b.rds")

# SPREAD list into workspace
spread(om5b)

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

# --- 1. TUNE cpue.ind + buffer.hcr(mult~wmean) {{{

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind, args=list(index=1, nyears=4)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(lim=0.25, bufflow=0.35, buffupp=0.60, sloperatio=0.10,
      dlow=0.85, dupp=1.15, metric="wmean", initac=36458))
))

system.time(
tes <- mp(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3))
)
  
# - EXPLORE grid

combs <- combinations(lim=0.10,
  bufflow=seq(0.20, 0.35, by=0.05),
  buffupp=seq(0.50, 0.65, by=0.05),
  sloperatio=c(0.10, 0.20))

explore <- mps(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  hcr=combs, parallel=FALSE)

performance(explore) <- performance(explore, statistics=statistics, metrics=mets,
  type="buffer-wmean")

perf <- performance(explore)

save(perf, file='model/explore_performance.rda')

# - TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(0.50, 0.75)), prob=0.6, tol=0.01)
)

# PLOT
plot(om, kobe60=tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="buffer-wmean", run="kobe60")

writePerformance(performance(tune))

# SAVE
saveRDS(tune, file="model/runs/om5b_buffer-wmean_tune_kobe60.rds")

# - TUNE for P(Kobe=green) = 70%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(0.50, 0.75)), prob=0.7, tol=0.01)
)

# PLOT
plot(om, kobe70=tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="buffer-wmean", run="kobe70")

# WRITE to db
writePerformance(performance(tune))

# SAVE
saveRDS(tune, file="model/runs/om5b_buffer-wmean_tune_kobe70.rds")

# }}}

# --- 2. TUNE cpue.ind + buffer.hcr(mult~mean) {{{

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
saveRDS(tune, file="model/runs/om5b_buffer-mean_tune_kobe60.rds")

# }}}

# --- 3. ROBUSTNESS runs cpue.ind + buffer.hcr(mult~wmean) {{{

# -- kobe60

# LOAD run
run <- readRDS("model/runs/om5b_buffer-wmean_tune_kobe60.rds")

# - OM5a

# LOAD om
om5a <- readRDS("data/om5a.rds")

# RUN mp
rob5a <- mp(om5a$om, om5a$oem, control=control(run), args=args(run))

# PLOT
plot(om, TUNE=run, ROB5a=rob5a)

# - OM6b

# LOAD om
om6b <- readRDS("data/om6b.rds")

# RUN mp
rob6b <- mp(om6b$om, om6b$oem, control=control(run), args=args(run))

# PLOT
plot(om, TUNE=run, ROB5a=rob5a, ROB6b=rob6b)

# - SAVE

rob <- FLmses(list('kobe60-om5a-robust'=rob5a, 'kobe60-om6b-robust'=rob6b),
  statistics=statistics, metrics=mets, type="buffer-wmean")

saveRDS(rob, file="model/runs/om5b_buffer-wmean_robust_kobe60.rds")

writePerformance(rob)

# }}}
