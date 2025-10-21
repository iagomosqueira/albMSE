# abc5.R - DESC
# /home/mosquia/Active/Doing/ABC_tuna+iotc/abc_tuna/v2/abc.R
# resample (h,M) from joint distribution
# resample sigmaR using informative conjugate prior

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: SET mpCtrl(est=shortcut.sa) as default
# TODO: SET mseCtrl(function) method


library(Rcpp)
library(FLCore)
library(ggplotFL)
library(parallel)
library(mvtnorm)
library(mse)
source("utilities.R")

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# NC by fleet [y, s, f]
# 

load("boot/data/alb_abcdata.rda")
load("boot/data/hmuprior.rda")
load("data/base.rda")

# --- SETUP {{{

# fnscale 
fnscale <- 1

# years and values for stock status priors

yof <- 1:ny # apply over-fishing prior penalty
sdof <- 0.5 # idea P(hmsyrat > 2) <= 0.05 essentially
ybmsy <- c(ny-1,ny)
mubmsy <- c(2.25,2)
sdbmsy <- c(0.35,0.35)
ydep <- 1
mudep <- 0.5
sddep <- 0.1

# rescale weight to tonnes

wta <- wta[]*1e-3
 
pobs <- apply(LFfits,2,function(x){x <- x/sum(x)}) # what we will fit to

# --- simulator

# - arguments
#   - parameters: R0, dep, h, epsr, selpars
#   - biology
#   - fishery

R0 <- 14e6
dep <- 0.5
h <- 0.8

# number of fisheries

nf <- dim(C)[3]

# number of distinct selectivity groups

nselg <- 5

# selectivity for each fishery

selidx <- c(1,2,3,4,5,5)

# fisheries with LF data

flf <- c(1,2,3,4,6)
nflf <- length(flf)
pobs <- apply(LFfits,2,function(x){x <- x/sum(x)}) # what we will fit to
pobs <- pobs[,flf]

# catch distro targets

pctarg <- C[1,,] / sum(C[1,,])

# set up selectivity parameters (all double normal)

smax <- c(120,125,85,85,115)
sL <- c(20,20,7,7,10)
sR <- c(35,30,25,15,30)
selpars <- cbind(smax,sL,sR)

# recruitment variations (ny-1)

epsr <- rep(0,ny-1)

# informative priors (roughly between 0.2-0.5 mean of 0.3)

musigma2R <- 0.3^2
cvx <- 0.4
mux <- 1/musigma2R
vx <- (mux*cvx)^2
alpR <- mux^2/vx
betR <- mux/vx
psigmaR <- sqrt(1/rgamma(1000,alpR,betR))
round(quantile(psigmaR,c(0.025,0.5,0.975)),2)

# recruitment season

srec <- 4

# sex ratio at birth (fiddy:fiddy)

psi <- 0.5

# dimensions

dms <- c(na,ns,nf,nselg)

## MCMC algorithm 

# set up MCMC controls for unconditional sampling of (h,M)

acphmu <- 0.25 # force acceptance rate at "optimal" MCMC value

# Gibbs sampling parameter groupings
# 1. B0 and dep
# 2. recruitment deviates
# 3. selectivity 

npar <- 2+ny-1+3*nselg
ngibbs <- 3
paridx <- list()
paridx[[1]] <- 1:2
paridx[[2]] <- 3:(ny+1)
paridx[[3]] <- (ny+2):npar
lidx <- unlist(lapply(paridx,length))

# SD in CPUE index

fcpue <- 1
scpue <- 1:4
sd.cpue <- rep(NA,length(scpue))

sdcpue <- mean(sd.cpue)

fcpue <- 1
scpue <- 1:4
sd.cpue <- rep(NA,length(scpue))
for(s in scpue) {
  idf <- data.frame(t=yrs,y=log(I[,s,fcpue]))
  ires <- loess(y~t,idf)
  sd.cpue[s] <- sd(residuals(ires))
}

sdcpue <- mean(sd.cpue)

# KLmax

KLmax <- 0.8 # consistent with minimum Neff = 20 multinomial

# seasonal q for CPUE (T or F)

seasonq <- TRUE

# catchability trend in q

qtrend <- FALSE

# burn-in and thinning factor
 
burn <- 10
thin <- 1

# }}}

# --- SAMPLER {{{

###################
# run the sampler #
###################

# set up initial h and M

hold <- hmu
Mold <- Mmu
sigmarold <- sqrt(musigma2R)

# set up initial guess parameter vector

parvecold <- c(log(R0),logit(dep),epsr,log(as.vector(selpars)))

# RW variance by Gibbs grouping

rwsd <- rep(0,npar)
rwsd[paridx[[1]]] <- c(0.1,0.05)
rwsd[paridx[[2]]] <- 0.08
rwsd[paridx[[3]]] <- 0.025

# TEST
nits1 <- 10 # total number of retained samples
system.time(zzz <- mcmc5.abc(nits1))
zzz$acp/nits1

# parallelised efficient version

parvecold <- zzz$pars[nits1,1:npar]
hold <- zzz$pars[nits1,npar+1]
Mold <- zzz$pars[nits1,npar+2]
sigmarold <- zzz$pars[nits1,npar+3]
nits <- 500
ncore <- 5
thin <- 100
mcnits <- floor(nits/ncore)

# SAVE image
save.image(file="data/om5b/image_abc5b.rda", compress="xz")

# RUN: 1.2 h
system.time(
  mczzz <- mclapply(rep(mcnits,ncore),mcmc5.abc,mc.cores=ncore)
)

# SAVE runs
save(mczzz, file="data/om5b/mcmc_abc5b.rda", compress="xz")

# EXTRACT
# load('data/image/abc5b.rda')
mcpars <- do.call(rbind, lapply(mczzz, '[[', 'pars'))
# BUG: RE-source Rcpp
mcvars <- get.mcmc2.vars(mcpars)

# SAVE output
save(mcpars, mcvars, C, file="data/om5b/mcvars_abc5b.rda", compress="xz")

# }}}

# --- CREATE om & oem {{{

load('data/om5b/mcmc_abc5b.rda')
load('data/om5b/mcvars_abc5b.rda')

# LOAD pcbar (catch proportions by fleet & season)
pcbar <- as.matrix(fread('data/pcbar.dat'))

# LOAD ALK [len, age, season, sex]
load('data/pla.rda')

# EXTRACT output for all iters
system.time(
  out <- mc.output(mcvars, C)
)

save(out, file="data/om5b/mcout_abc5b.rda", compress="xz")

its <- dims(out$m)$iter

# - FLBiol (stock.n, m)

bio <- propagate(window(sbio, start=2000), its)

n(bio) <- out$stock.n
m(bio) <- out$m

sr(bio) <- predictModel(model=bevholtss3()$model, params=out$srpars)

# spwn, mid Q4
spwn(bio)[,,,4] <- 0

# FIX mat BUG: CHECK readFLSss3 +ss3om 
mat(bio)[,,'F',1:4] <- mat(bio)[,,'F',1]

# DEVIANCES
deviances(bio) <- out$deviances

# exp(residuals(rec(bio)[,,,4],
#  expand(predict(sr(bio), ssb=out$ssb) / 2, unit=c('F', 'M')))), 1)

# - FLFisheries (catch.n, catch.sel)

# FLCatch(es)
cas <- Map(function(x, y) FLCatch(landings.n=x, landings.wt=wt(bio),
  catch.sel=y, discards.n=x %=% 0, discards.wt=wt(bio)),
  x=divide(out$catch.n, 5), y=divide(out$catch.sel %*% (out$catch.n %=% 1), 5))

# FLFisheries
fis <- FLFisheries(lapply(cas, function(x)
  FLFishery(effort=unitSums(catch(cas[[1]])) %=% 0, ALB=x)))

names(fis) <- c(paste0("LL", 1:4), "PS", "Other")

om <- FLombf(biols=FLBiols(ALB=bio), fisheries=fis,
  refpts=FLPars(ALB=out$refpts))

method(projection(om)) <- fwdabc.om
args(projection(om)) <- list(pla=pla, pcbar=pcbar)

# - BUILD oem

# idx: FLIndexBiomass by season, with sel.pattern by sex

# BUG: mc.output to output FLQuants by fiosheru, not 'area'
sp <- expand(divide(out$catch.sel, 5)[[1]], year=2000:2020)
dimnames(sp)$area <- 'unique'

NW <- FLIndexBiomass(index=out$index.hat %*% out$index.q,
  index.q=expand(out$index.q, year=2000:2020),
  sel.pattern=sp,
  catch.wt=wt(biol(om)),
  range=c(startf=0.5, endf=0.5))

# stk: no units
oem <- FLoem(observations=list(ALB=list(idx=FLIndices(NW=NW),
  stk=simplify(stock(om)[[1]], 'unit'))), method=sampling.oem)

survey(observations(oem)$ALB$stk, observations(oem)$ALB$idx)

# TODO: verify(oem, om)

# SAVE
save(om, oem, file='data/om5b/om5b-raw.rda', compress='xz')

# }}}

# --- EXTEND {{{

load('data/om5b/om5b-raw.rda')
load('data/om5b/mcout_abc5b.rda')

om <- fwdWindow(om, end=2045)

deviances(om)[, ac(2021:2045),,4] <- rlnormar1(n=dims(om)$it, meanlog=0, 
  sdlog=out$sigmar, rho=out$rho, years=2021:2045)

oem <- fwdWindow(oem, end=2045)

# SAVE
save(om, oem, file='data/om5b.rda', compress='xz')

# }}}

# UPDATE NC {{{

load('data/om5b.rda')
load('data/alb_2025_nw.rda')

# GET new NC

ctrl <- as(FLQuants(catch=unitSums(seasonSums(catch(ss25$stk)))[, ac(2010:2023)]), 
  'fwdControl')

# - UPDATE using fwdabc.om

system.time(
nom <- fwdabc.om(om, ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# CHECK catch match
ctrl$value
catch(nom)[, ac(2010:2023)]

# CHECK indices

onw <- observations(oem)$ALB$idx[[1]][, ac(2000:2023)]
nnw <- survey(stock(nom)[[1]][, ac(2000:2023)],
  observations(oem)$ALB$idx[[1]][, ac(2000:2023)]

plot(index(onw), index(nnw))

# UPDATE

om <- nom

observations(oem)$ALB$idx$NW <- nnw

observations(oem)$ALB$stk[, ac(2010:2023)] <- nounit(stock(nom)[[1]])[, ac(2010:2023)]

# SAVE

om5b <- list(om=om, oem=oem)

system.time(
save(om5b, file="data/om5b_updated.rda", compress="xz")
)
system.time(
saveRDS(om5b, file="data/om5b_updated.rds", compress="xz")
)

system.time(
qs_save(om5b, file="data/om5b_updated.qs2")
)

# }}}

# PROJECTIONS {{{

# FWD(C=C0) 

ctrl <- fwdControl(year=2024:2045, biol=1, quant='catch', value=0)

system.time(
fom_c0 <- fwdabc.om(iter(om, seq(100)), ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# FWD(C=MSY2025) 

ctrl <- fwdControl(year=2024:2045, biol=1, quant='catch', value=ss25$rps$MSY)

system.time(
fom_cmsy <- fwdabc.om(iter(om, seq(100)), ctrl, pcbar=args(projection(om))$pcbar,
  pla=args(projection(om))$pla, verbose=TRUE)$om
)

# SAVE
save(fom_c2023, fom_c0, fom_cmsy, file="data/om5b/fwd_om5b.rda", compress="xz")

# }}}
