# explore_bufferdelta.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/test/explore_bufferdelta.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD om
load('data/om5b/om5b_updated.rda')

# EXTRACT
om <- iter(om, seq(50))
oem <- iter(oem, seq(50))

# bufferdelta.hcr {{{

timeMeans <- function(x)
  seasonMeans(yearMeans(x))

zscore <- function(x, mean=yearMeans(x), sd=sqrt(yearVars(x)))
  (x %-% mean) %/% sd

bufferdelta.hcr <- function(stk, ind, target=1, metric='zscore',
  width=1, bufflow=target - width, buffupp=target + width,
  lim=target - 2 * width, sloperatio=0.20, initac=NULL,
  dupp=NULL, dlow=NULL, all=TRUE, ..., args, tracking) {

  # EXTRACT args
  ay <- args$ay
  iy <- args$iy
  data_lag <- args$data_lag
  man_lag <- args$management_lag
  frq <- args$frq

  # SET data year
  dy <- ay - data_lag
  # SET control years
  cys <- seq(ay + man_lag, ay + man_lag + frq - 1)

  # COMPUTE metric
  met <- mse::selectMetric(metric, stk, ind)
  met <- window(met, start=dy, end=dy)
  
  # COMPUTE HCR multiplier if ...
  # BELOW lim
  hcrm <- ifelse(met <= lim, ((lim/met) ^ 2) / 2,
    # BETWEEN lim and bufflow
    ifelse(met < bufflow,
      (0.5 * (1 + (met - lim) / (bufflow - lim))),
    # BETWEEN bufflow and buffupp
    ifelse(met < buffupp, 1, 
    # ABOVE buffupp
      1 + sloperatio * 1 / (2 * (bufflow - lim)) * (met - buffupp))))

  # GET previous TAC from last hcr ...
  if(is.null(initac)) {
    pre <- tracking[[1]]['hcr', ac(ay)]
    # ... OR catch
    if(all(is.na(pre)))
      pre <- unitSums(areaSums(seasonSums(catch(stk)[, ac(ay - args$data_lag)])))
  } else {
    pre <- FLQuant(initac, iter=args$it)
  }

  # SET TAC as tac = B * (1 - exp(-fm * hcrm * (F / FMSY))
  out <- pre * hcrm

  # TRACK initial target
  track(tracking, "tac.hcr", cys) <- out

  # APPLY limits, always or if met < trigger
  if(!is.null(dupp)) {
    if(all) {
      out[out > pre * dupp] <- pre[out > pre * dupp] * dupp
    } else {
      out[out > pre * dupp & met < bufflow] <- pre[out > pre * dupp & met <
        bufflow] * dupp
    }
  }

  if(!is.null(dlow)) {
    if(all) {
      out[out < pre * dlow] <- pre[out < pre * dlow] * dlow
    } else {
      out[out < pre * dlow & met < bufflow] <- pre[out < pre * dlow & met <
        bufflow] * dlow
    }
  }

  # CONTROL
  ctrl <- fwdControl(
    # TARGET for frq years
    c(lapply(cys, function(x) list(quant="catch", value=c(out), year=x)))
  )
	
  list(ctrl=ctrl, tracking=tracking)
}
# }}}

# SETUP control

ref <- index(observations(oem)$ALB$idx[[1]])[, ac(c(2000:2005, 2015:2020))]

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind,
    args=list(index=1, mean=yearMeans(ref), sd=sqrt(yearVars(ref)))),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(bufflow=0.5, buffupp=1.5, sloperatio=0.15, dlow=0.85, dupp=1.15,
      metric="zscore", initac=42000))
))

# EXPLORE
exp <- mps(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3),
  hcr=combinations(list(bufflow=seq(0, 1, length=3),
    buffupp=seq(1.5, 4, length=3),
    sloperatio=seq(0.10, 0.40, length=3))))

save(exp, file='test/explore_bufferdelta.rda', compress='xz')

plotMetrics(OM=window(om, end=2024), EXP=window(om(exp[[1]]), start=2024))




