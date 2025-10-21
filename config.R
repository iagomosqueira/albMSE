# config.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/abc_tuna/v3/config.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(Rcpp)
library(FLCore)
library(ggplotFL)
library(parallel)
library(mvtnorm)
library(mse)
source("utilities.R")

library(TAF)

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")
