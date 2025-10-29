# test_bufferdelta_hcr.R - DESC
# albMSE/test_bufferdelta_hcr.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD om
load('data/om5b_updated.rda')

om <- iter(om5b$om, seq(50))
oem <- iter(om5b$oem, seq(50))

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP

iy <- 2023
fy <- 2042
ty <- seq(iy + 11, iy + 15)

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind,
    args=list(index=1, ayears=4)),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=1, width=1, buffupp=3, sloperatio=0.15, dlow=0.85, dupp=1.15,
      metric="wmean", initac=42000))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3))#, .DEBUG=TRUE)

# EXPLORE width
exp <- mps(om, oem, ctrl=ctrl, args=list(iy=iy, fy=fy, frq=3),
  hcr=list(width=seq(0.2, 3, length=15)))

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
