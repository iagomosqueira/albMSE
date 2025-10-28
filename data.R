# data.R - DESC
# /home/mosqu003/Active/ALB_MSE-IOTC/albMSE/data.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library("data.table")

# --- IOTC catch data {{{

iotc_catch <- fread("https://iotc.org/sites/default/files/documents/2025/10/IOTC-DATASETS-2025-10-22-NC-SCI_1950-2024.zip")

alb_catch <- iotc_catch[SPECIES_CODE == "ALB" & YEAR %in% 2010:2023]

# NC
nominal_catch <- alb_catch[, .(catch=sum(CATCH)), by=.(year=YEAR)][order(year)]

# SAVE
save(nominal_catch, alb_catch, file='data/iotc_alb_catch.rda', compress='xz')

# }}}

library(ss3om)

# --- BASE 2020 SS3 SA {{{

path <- file.path("boot", "data", "base")

res <- readOutputss3(path)

# FLBiol + FLFisheries
bfs <- buildFLBFss330(res)
sbio <- bfs$biol
sfis <- bfs$fisheries

# FLStock
sstk <- buildFLSss330(res)
range(sstk, c("minfbar", "maxfbar")) <- c(1, 12)

# FLIndices
idss <- buildFLIBss330(res)

# refpts
rpss <- buildFLRPss330(res)

# SRR
srss <- buildFLSRss3(res)

base <- list(bio=sbio, fis=sfis, stk=sstk, ids=idss, rps=rpss, srr=srss)

# SAVE
save(base, file='data/base.rda', compress='xz')

# }}}

# --- 2023 SS3 NW SA {{{

path <- file.path("boot", "data", "ALB_2025_NW_Basecase_24July2025")

res <- readOutputss3(path)

# FLBiol + FLFisheries
bfs <- buildFLBFss330(res)
sbio <- bfs$biol
sfis <- bfs$fisheries

# FLStock
sstk <- buildFLSss330(res)
range(sstk, c("minfbar", "maxfbar")) <- c(1, 12)

# FLIndices
idss <- buildFLIBss330(res)

# refpts
rpss <- buildFLRPss330(res)

# SRR
srss <- buildFLSRss3(res)

ss25 <- list(bio=sbio, fis=sfis, stk=sstk, ids=idss, rps=rpss, srr=srss)

# SAVE
save(ss25, file='data/alb_2025_nw.rda', compress='xz')

# }}}
