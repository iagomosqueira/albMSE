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

# --- 1. TUNE shortcut.sa + fixedC.hcr {{{

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

# PLOT
plot(om, cc=tune)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="fixedC", run="tune_kobe60")

# WRITE to performance table
writePerformance(performance(tune))

# SAVE
saveRDS(tune, file="model/runs/om5b_buffer-mean_tune_kobe60.rds")

# GET tuned C = 40 625 t
args(control(tune)$hcr)$ctrg

# }}}
