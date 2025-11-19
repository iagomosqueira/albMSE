# output.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/output.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")


# RE-LABEL performance elements

perf <- readPerformance()

# BUG: REMOVE very large HRMSY for plots
perf[statistic == 'HRMSY', data := ifelse(data > 4, 4, data)]

perf <- labelPerformance(perf, labels=list(
  'om5b'='OM base',
  'om6b'='OM 1% Q',
  'om5a'='OM SW',
  'om5b_fixedC_tune_kobe60'="Constant catch",
  'om5b_buffer-wmean_kobe60'="Tuned Kobe 60%",
  'om5b_buffer-wmean_kobe70'="Tuned Kobe 70%",
  'om5a_buffer-wmean_kobe60-om5a-robust'="Robust SW Kobe 60",
  'om6b_buffer-wmean_kobe60-om6b-robust'="2Robust 1% LL Q Kobe 60"))

writePerformance(perf, file="output/performance.dat.gz")

# TABLES



# --- RE-STORE performance in db {{{

# OMs

lapply(list.files('data', pattern="*.rds"), function(x) {

  obj <- readRDS(file.path("data", x))

  writePerformance(performance(obj$om, metrics=mets, years=2000:2024,
    statistics[c("SB", "SB0", "SBMSY", "R", "HRMSY", "C")]))
  }
)


# MPs

lapply(setNames(nm=list.files("model/runs")), function(x) {
  obj <- readRDS(file.path("model", "runs", x))
  writePerformance(performance(obj))
})

# }}}
