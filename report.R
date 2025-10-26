# report.R - DESC
# abc_tuna/v3/report.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

source("config.R")

library(mseviz)

mkdir("report")

# --- data.R {{{

# - COMPARE SS3 2022-2025

load('data/base.rda')
load('data/alb_2025_nw.rda')

sas <- FLStocks('2022'=noseason(nounit(base$stk)),
  '2025'=noseason(nounit(ss25$stk)))

# PLOT FLStocks
taf.png("data_sas_compare.png")
plot(sas)
dev.off()

# COMPUTE diffs in rps
(ss25$rps$SB0 / base$rps$SB0) * 100
(ss25$rps$R0 / base$rps$R0) * 100

# COMPARE indices

id1_2020 <- Reduce(join, lapply(base$ids[1:4], index))
id1_2025 <- Reduce(join, lapply(ss25$ids[1:4], index))

taf.png("data_id1_compare.png")
plot(id1_2020, id1_2025) +
  theme(legend.position = "bottom") +
  scale_fill_manual(name='',
    labels=c('2022', '2025'),
    values=c(v1=flpalette[1], v2=flpalette[2])) +
  scale_color_manual(name='NW (LLCPUE1)',
    labels=c('2022', '2025'),
    values=c(v1=flpalette[1], v2=flpalette[2]))
dev.off()
 
# - COMPARE OMs2

load('data/om5b.rda')
load('data/om5b_updated.rda')

# PLOT OM

taf.png("data_om_basecase.png")
plotMetrics('Base case'=window(om, end=2020), maxhr=2)
dev.off()

# PLOT 2022 and 2025 OMs
taf.png("data_oms_compare.png")
plotMetrics('2022'=window(om, end=2020), '2025'=window(om5b$om, end=2023), maxhr=2)
dev.off()

# PLOT comparison updated OM & SS3 SAs

d1 <- as.data.frame(metrics(window(base$stk, start=2000), metrics=list(
  SB=function(x) unitSums(seasonSums(ssb(x)[,,1,1])),
  R=function(x) unitMeans(seasonSums(rec(x)[,,1,1])))))=3

d2 <- as.data.frame(metrics(window(ss25$stk, start=2000), metrics=list(
  SB=function(x) unitSums(seasonSums(ssb(x)[,,1,1])),
  R=function(x) unitMeans(seasonSums(rec(x)[,,1,1])))))

taf.png("data_om_compare_sa.png")
plot(FLQuants(lapply(mets[c(1,4)], function(x) window(x(om5b$om), end=2023)))) +
  ylim(0, NA) +
  geom_flquantiles(data=d1, colour="blue") +
  geom_flquantiles(data=d2, colour="darkgreen")
dev.off()

#

taf.png("data_om_compare_sb0.png")
plot(refpts(om)['SB0',]) +
  geom_vline(xintercept=base$rps$SB0, colour="darkgreen") +
  annotate(geom='point', x=base$rps$SB0, y=0, color="darkgreen", size=3) +
  geom_vline(xintercept=ss25$rps$SB0, colour="blue") +
  annotate(geom='point', x=ss25$rps$SB0, y=0, color="blue", size=3) +
  ggtitle("Virgin SSB")
dev.off()


# PLOT OM projections

load("data/om5b/fwd_om5b.rda")

library(ggh4x)

taf.png("data_om_fwd.png")
plotMetrics(OM=window(om, end=2023), C0=window(fom_c0, start=2023),
  MSY=window(fom_cmsy, start=2023))
dev.off()

# }}}

# --- model.R {{{

# LOAD om
load('data/om5b_updated.rda')
om <- iter(om5b$om, seq(100))

# LOAD results
load("model/model_cpue_buffer.rda")

# LOAD performance
perf <- readPerformance()

# PLOT constant catch tuned MP
taf.png("model_ccatch.png")
plotMetrics(OM=window(om, end=2023),
  'Constant catch'=window(om(res[[1]]), start=2023)) +
  geom_vline(xintercept=ISOdate(c(2034, 2038), 1, 1), linetype=3, alpha=0.8)
dev.off()

# ASSEMBLE dat: mean performance 2034-2038
dat <- perf[year %in% seq(2034, 2038), .(data=mean(data)),
  by=.(statistic, name, desc, mp)]

plotBPs(dat)

# PLOT runs
taf.png("model_buffer_catch.png")
plotMetrics(OM=window(om, end=2023),
  'buffer(C~CPUE)'=window(om(res[[2]]), start=2023)) +
  geom_vline(xintercept=ISOdate(c(2034, 2038), 1, 1), linetype=3, alpha=0.8)
dev.off()

# PLOT runs
taf.png("model_runs_compare.png")
plotMetrics(OM=window(om, end=2023),
  'Constant catch'=window(om(res[[1]]), start=2023),
  'buffer(C~CPUE)'=window(om(res[[2]]), start=2023)) +
  geom_vline(xintercept=ISOdate(c(2034, 2038), 1, 1), linetype=3, alpha=0.8)
dev.off()

# }}}

# RENDER

render('report_wpm_2025.Rmd', output_dir='report',
  output_file='IOTC-2025-WPM16-11_ALB_MSE.pdf')

render('presentation_wpm_2025.Rmd', output_dir='report', output_file='presentation-IOTC-2025-WPM16-11_ALB_MSE.pdf')
