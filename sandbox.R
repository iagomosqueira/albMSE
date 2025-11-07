# sandbox.R - DESC
# albMSE/sandbox.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


source("config.R")

# [ ] select refyrs
# [ ] TUNE refyrs, fixed width
# [ ] EXPLORE width vs. upp/low limits

# --- [ ] plot_buffer.hcr(results=TRUE)

target=1
width=0.5
lim=max(target * 0.10, target - 2 * width)
bufflow=max(lim, target - width)
buffupp=target + width
sloperatio=0.15

args <- list(target=target, width=width, lim=lim,
  bufflow=bufflow, buffupp=buffupp, sloperatio=sloperatio)

plot_buffer.hcr(args)

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind, args=list(index=1, ayears=4)),
  # HCR
  hcr = mseCtrl(method=bufferdelta.hcr,
    args=list(target=1, width=1, buffupp=3, sloperatio=0.15, dlow=0.85, dupp=1.15,
      metric="zscore", initac=42000))
))

plot_buffer.hcr(ctrl)

plot_buffer.hcr(ctrl$hcr)

plot_buffer.hcr(args)


bdargs <- function(target=1, width=0.5) {

  list(target=target, width=width, lim=max(target * 0.10, target - 2 * width),
  bufflow=max(lim * 1.50, target - width), buffupp=target + width, sloperatio=0.15)
}

plot_buffer.hcr(bdargs(target=1))
plot_buffer.hcr(bdargs(target=2))
plot_buffer.hcr(bdargs(target=3))
plot_buffer.hcr(bdargs(target=4))

plot_buffer.hcr(bdargs(width=0.20))
plot_buffer.hcr(bdargs(width=0.50))
plot_buffer.hcr(bdargs(width=0.75))


# --- [X] COMPARE idx ~ kobe(sa)

# LOAD
load('data/base.rda')

out <- ss3om::readOutputss3('boot/initial/data/base/')

# EXTRACT Kobe
kobe_dat <- data.table(out$Kobe)

# GREEN years
kobe_dat[, kobe := B.Bmsy >= 1 & F.Fmsy < 1]

# NW CPUE years wih catch ~ CC tuned value
ref <- index(base$ids[[1]])[, unitSums(seasonSums(catch(base$stk)[, ac(1975:2020)])) < 45000 & unitSums(seasonSums(catch(base$stk)[, ac(1975:2020)])) > 35000]

yearMeans(ref)
sqrt(yearVars(ref))


