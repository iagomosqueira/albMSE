# data.R - LOAD SS3 model and add ABC output
# abc_tuna/om/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


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
