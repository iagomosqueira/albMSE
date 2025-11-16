# output.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/output.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# STORE OM perf in db

lapply(list.files('data', pattern="*.qs2"), function(x) {

  obj <- qs_read(file.path("data", x))

  writePerformance(performance(obj$om, metrics=mets, years=2000:2024,
    statistics[c("SB", "SB0", "SBMSY", "R", "HRMSY", "C")]))
  }
)

readPerformance()
