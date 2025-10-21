# report.R - DESC
# abc_tuna/v3/report.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

source("config.R")

mkdir("report")

# --- data.R {{{

# - OM update to 2023

load('data/om5b.rda')


# - OM update to 2023

oom <- window(om, end=2020)

load('data/om5b_updated.rda')

uom <- window(om, end=2023)

# PLOT comparison of conditioned and updated OMs

taf.png("data_oms_compare.png")
plotMetrics(COND=oom, UPDATE=uom)
dev.off()

# PLOT compartison updated OM & SS3 SA

load('data/base.rda')
load('data/alb_2025_nw.rda')

taf.png("data_om_compare_sa.png")
plot(ssb(uom),
  window(ssb(ss25$stk)[,,1,1], start=2000),
  window(ssb(base$stk)[,,1,1], start=2000)) +
  ylim(0, NA) +
  scale_color_manual(name='Regression Model',
    breaks=c('v1', 'v2', 'v3'),
    values=c('v1'='pink', 'v2'=flpalette[2], 'v3'=flpalette[4])) +
  ylab("SSB (t)")
dev.off()

# PLOT OM projections

load("data/om5b/fwd_om5b.rda")

library(ggh4x)

taf.png("data_om_fwd.png")
plotMetrics(OM=window(om, end=2023), C0=window(fom_c0, start=2023),
  MSY=window(fom_cmsy, start=2023))
dev.off()

# }}}
