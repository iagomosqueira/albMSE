# model_jabba_buffer.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/model_jabba_buffer.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

source("config.R")

# LOAD 100 iter objects
qs_readm("data/om5b.qs2")

library(ss3om)

ssout <- readOutputss3("boot/data/base")

mvnl <- ss3diags::SSdeltaMVLN(ssout)

stkr <- FLRef::ss2FLStockR(mvnl)

  SBMSY=refpts(om)['SBMSY', 4, drop=TRUE],

omrps <- FLPar(MSY=refpts(om)['MSY', 1, drop=TRUE],
  SBMSY=refpts(om)['SBMSY', 4, drop=TRUE],
  SB0=refpts(om)['SB0', 4, drop=TRUE])

save(ssout, stkr, mvnl, omrps, file="test/iotc_alb_refpts.rda", compress="xz")
