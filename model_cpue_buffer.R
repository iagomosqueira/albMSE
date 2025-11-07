# model_buffer.R - DESC
# abc_tuna/v3/model_buffer.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# LOAD 100 iter objects
qs_readm("data/om5b.qs2")

# RESET method JIC
method(projection(om)) <- fwdabc.om

# SETUP
iy <- 2024
fy <- 2042
ty <- seq(iy + 11, iy + 15)

# PLAN
plan(multicore, workers=10)

# --- TUNE shortcut.sa + fixedC.hcr {{{

# SET control
ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=shortcut.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr,
    args=list(ctrg=catch(om)[,'2024'][[1]]))
))

# test <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2026, frq=3), .DEBUG=TRUE)

# TUNE for P(Kobe=green) = 60%
system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3), 
  statistic=statistics["green"], metrics=mets, years=iy + c(11, 15),
  tune=list(ctrg=c(35000, 50000)), prob=0.6, tol=0.01, maxit=12)
)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="fixedC", run="tune_kobe60")

# WRITE to table
writePerformance(performance(tune))

# STORE in results
save(tune, file="model/runs/om5b_fixedC_tune_kobe60.rda", compress="xz")

# CHECK Kobe green
performance(tune)[statistic == 'green' & year %in% ty, mean(data)]
performance(tune)[statistic == 'green', mean(data), by=year]

# GET tuned C = 40 625 t
args(control(tune)$hcr)$ctrg

# PLOT
plot(om, "Constant catch Kobe 60%"=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# }}}

# --- DOES NOT TUNE -- TUNE cpues.ind + bufferdelta.hcr(mult~zscore) {{{

# EXPLORE NW index to get zscore reference years

load('data/base.rda')

nwi <- Reduce(join, lapply(base$ids[1:4], index))

# NW CPUE years with catch ~ CC tuned value
ref <- seasonMeans(nwi[, unitSums(seasonSums(catch(base$stk)[, ac(1975:2020)])) < 45000 & unitSums(seasonSums(catch(base$stk)[, ac(1975:2020)])) > 35000])

yearMeans(ref)
sqrt(yearVars(ref))

exp(zscore(seasonMeans(nwi)[, ac(2000:2020)], mean=yearMeans(ref), 
  sd=sqrt(yearVars(ref))))

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind,
    args=list(index=1, mean=yearMeans(ref), sd=sqrt(yearVars(ref)), nyears=4)),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=1, buffupp=1.8, bufflow=0.8, sloperatio=0.15, lim=0.5,
      metric="zscore", initac=42000))))

# TEST run
test <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2032, frq=3), .DEBUG=FALSE)

# PLOT
plot(om, TEST=test) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# KOBE performance
performance(test, statistics=statistics['green'], metrics=mets)[year %in% ty, mean(data)]

# - TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(buffupp=c(1.5, 3)), prob=0.6, tol=0.02, maxit=16)
)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="tune") #, run="kobe60")

# KOBE performance
performance(tune)[statistic=='green' & year %in% ty, mean(data), by=mp]

# PLOT
plot(om, K60=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# }}}

# --- TUNE cpue.ind + buffer.hcr(C~zscore) {{{

# SET control

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind,
    args=list(index=1, mean=yearMeans(ref), sd=sqrt(yearVars(ref)), nyears=4)),
  # HCR
  hcr = mseCtrl(method=buffer.hcr,
    args=list(target=42070, lim=0.5, bufflow=0.8, buffupp=1.2, sloperatio=0.15,
      dlow=0.85, dupp=1.15, metric="zscore", initac=catch(om)[,'2024'][[1]]))
))

# RUN
tes <- mp(om, oem, ctrl=ctrl, args=list(iy=iy, fy=2032, frq=3))#, .DEBUG=TRUE)

# TUNE for P(Kobe=green) = 60%

system.time(
tune <- tunebisect(om, oem=oem, control=ctrl, args=list(iy=iy, fy=fy, frq=3),
  statistic=statistics["green"], metrics=mets, years=ty,
  tune=list(target=c(25000, 55000)), prob=0.6, tol=0.01, maxit=12)
)

# PLOT
plot(om, CCK60=tune) +
  geom_vline(xintercept=ISOdate(c(ty[1], ty[length(ty)]), 1, 1), linetype=3, alpha=0.8)

# COMPUTE performance
performance(tune) <- performance(tune, statistics=statistics, metrics=mets,
  type="zscore_NW_buffer_C", run="tune_kobe60")

performance(tune)[statistic == 'green' & year %in% ty, mean(data)]
performance(tune)[statistic == 'green', mean(data), by=year]

# WRITE to table
writePerformance(performance(tune))

# STORE in results
res[["abc5b_zscore_NW_buffer_C_tune_kobe60"]] <- tune

# GET value of tuned argument
args(control(tune, 'hcr'))$target

# }}}
