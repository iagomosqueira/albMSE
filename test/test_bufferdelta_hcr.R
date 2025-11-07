# test_bufferdelta_hcr.R - DESC
# albMSE/test_bufferdelta_hcr.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD om
load('data/om5b_updated.rda')

om <- iter(om5b$om, seq(100))
oem <- iter(om5b$oem, seq(100))

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP

iy <- 2023
fy <- 2042
ty <- seq(iy + 11, iy + 15)

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind, args=list(index=1, ayears=4)),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=1, width=1, buffupp=3, sloperatio=0.15, dlow=0.85, dupp=1.15,
      metric="zscore", initac=42000))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3))#, .DEBUG=TRUE)

# EXPLORE

# target	lim	 bufflow	buffupp	sloperatio
# 1	      0.25 0.50	    1.20	   0.15
# 1	      0.20 0.40     1.30     0.15
# 1	      0.15 0.30     1.40     0.15
# 1	      0.10 0.20     1.50     0.15

opts <- list(
    target=rep(1, 4),
    lim=seq(0.25, 0.10, length=4),
    bufflow=seq(0.50, 0.20, length=4),
    buffupp=seq(1.20, 1.50, length=4))

exp00 <- mps(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2036, frq=3), hcr=opts)

performance(exp00) <- performance(exp00, statistics=statistics, metrics=mets,
  type='explore_buffer')

performance(exp00)[statistic == 'green' & year %in% ty, mean(data), by=mp]

plot(om, exp00)

# PLOT
plotMetrics(TES=window(om(tes), start=2023))
plotMetrics(OM=window(om, end=2023), TES=window(om(tes), start=2023))

# TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(1, 4)), prob=0.6, tol=0.01, maxit=12)
)
