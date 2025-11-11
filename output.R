# output.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/output.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

# STORE in db
writePerformance(performance(om, statistics[c("SB", "SB0", "SBMSY", "R", "HRMSY", "C")],
  metrics=mets, years=2000:2024), overwrite=TRUE)

# PLOTS

