# config.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/abc_tuna/v3/config.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(TAF)
library(FLCore)
library(ggplotFL)
library(mse)
library(Rcpp)
library(mvtnorm)
library(parallel)
library(ggh4x)
library(qs2)

source("utilities.R")

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# SETUP global progressr handlers
handlers(global=TRUE)
