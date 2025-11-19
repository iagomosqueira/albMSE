# report.R - DESC
# abc_tuna/v3/report.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

source("config.R")

library(mseviz)

mkdir("report")

# TIMING
iy <- 2024
ty <- seq(iy + 11, iy + 15)

# VLINES for tuning years (ty)
tperiod <- geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1),
  linetype=3, linewidth=0.75, alpha=0.8)

# SELECT iters for plots worms
worms <- sample(100, 5)

# --- data.R 

# -- McMC output {{{

load("data/om5b/mcvars_abc5b.rda")

# }}}

# -- CONDITIONED OMs

# -- COMPARE SS3 2022-2025 {{{

load('data/base.rda')
load('data/alb_2025_nw.rda')

sas <- FLStocks('2022'=noseason(nounit(base$stk)),
  '2025'=noseason(nounit(ss25$stk)))

# PLOT CATCH proportions by quarter

catch_props <- FLQuants("2022"=unitSums(catch(base$stk))[, ac(2010:2020)] %/%
  seasonSums(unitSums(catch(base$stk))[, ac(2010:2020)]),
  "2025"=unitSums(catch(ss25$stk))[, ac(2010:2023)] %/%
  seasonSums(unitSums(catch(ss25$stk))[, ac(2010:2023)])
)

taf.png("data_sas_catchprops.png")
ggplot(catch_props, aes(x=year, y=data, fill=season)) +
  geom_bar(stat='identity', position='fill') +
  facet_wrap(~qname, ncol=1) +
  xlab("") + ylab("Proportion of annual catch") +
  scale_y_continuous(labels = scales::percent)
dev.off()

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

# }}}
 
# -- EXTENDED OMs {{{

# - OM5b

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

# }}}

# -- OM projections {{{

load("data/om5b/fwd_om5b.rda")

library(ggh4x)

taf.png("data_om_fwd.png")
plotMetrics(OM=window(om, end=2023), C0=window(fom_c0, start=2023),
  MSY=window(fom_cmsy, start=2023))
dev.off()

# }}}

# --- model.R

# LOAD performance
perf <- readPerformance(file="output/performance.dat.gz")

# -- MP 0: om5b + shortcut.sa + fixedC.hcr 60% Kobe green {{{

# GET dataset
dat <- perf[om == 'om5b' & mp %in% c("", "om5b_fixedC_tune_kobe60")]

# PLOT MP run
taf.png("om5b_fixedC_tune_kobe60-timeseries.png")
plotTimeSeries(dat, statistics = c("SB", "R", "C", "HRMSY"), worms=TRUE) +
  tperiod
dev.off()

# }}}

# -- MP 1: om5a + cpue.ind(NW) + buffer.hcr(mult~wmean) 60% Kobe green {{{

# GET results
dat <- perf[om == 'om5b' & mp %in% c("", "om5b_buffer-wmean_kobe60")]

load('model/runs/om5b_buffer-wmean_tune_kobe60.rda')

# PLOT MP run
plotTimeSeries(dat, iters=worms) + tperiod

# PLOT HCR & future observations TODO: ADD decisions
plot_buffer.hcr(sat$control$hcr, xlim=1.5) +
  geom_point(data=data.table(met=c(index(observations(oem(tune))$ALB$idx[[1]])), out=0),
    alpha=0.01)

# PLOT index
plot(observations(oem(tune))$ALB$idx[[1]])

# PLOT Kobe
kobeMPs(dat[mp != ""], x="SBMSY", y="HRMSY") +
  ylim(0, 2)

# PLOT Kobe time series

kobeTS(dat[statistic %in% c('green', 'red', 'orange', 'yellow') & mp != ""])

dat[statistic %in% c('green', 'red', 'orange', 'yellow') & mp != ""]

kobeTS(dat[statistic %in% c('green', 'red', 'orange', 'yellow') & mp != "",
  .(data=mean(data)), by=.(statistic, year, mp)]) +
  geom_hline(yintercept=0.6, linetype=2, color='white')

# }}}

# -- MP 00 {{{

# PLOT MP run

# PLOT HCR

# PLOT index

# PLOPT Kobe

# PLOT Kobe time series

# }}}

# --- RENDER

render('report.Rmd', output_dir='report', output_file='report.pdf')

render('presentation_wpm_2025.Rmd', output_dir='report',
  output_file='presentation-IOTC-2025-WPM16-11_ALB_MSE.pdf')
