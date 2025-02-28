# config.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/albMSE/config.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

source('utilities.R')

# metrics {{{
mets <- list(
  C=function(x) lapply(catch(x), function(y) seasonSums(unitSums(y))),
  R=function(x) lapply(rec(x), function(y) unitSums(y[,,,4])),
  SB=function(x) lapply(biols(x), function(y)
    unitSums(quantSums(n(y) * wt(y) * mat(y))[,,,4])))
# }}}


data(statistics)
