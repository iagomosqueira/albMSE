---
title:
author: "Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>"
tags:
---


# 

- boot
  - initial/data: SS3 inputs
  - data.R: RIUN ss3 base
   - data: SS3 model run
- data.R: PREPARE inputs for abc
- model_.R: RUN ABC
- output.R: LOAD base SS3 run & ABC results
 

#

- mcpars
- mcvars
  - N
  - H
  - sela
  - Ihat
  - Rtot
  - SSB
  - dep
  - dbmsy
  - Cmsy
  - hmsyrat
  - LFhat
  - B0
  - R0
  - M
  - h
  - lnq
  - rho
- out
  - stock.n
  - m
  - catch.n
  - catch.sel
  - ssb
  - dep
  - index.q
  - srpars
  - srpars
  - refpts
  - hr
  - rec
  - index.hat
  - hra

# pdynlfcpue

- dms: dimensions [a,s]
- srec: recruitment season 
- R0, hh, spr0: sr(om)
- psi
- epsrx
- M: m(om)
- mata: mat(om) [a, s, u]
- wta: wt(om) [a, s, u]
- sela: catch.sel(om) [a, s, u, f]
- nvec: n(om)
- cvec: catch(om)
- pla [l, a, s, u]
- fcpue

- N
- H [y,s,f]


# fwdabc.om

- CALL at om & oem
- HR / F in object
  - STORE or recompute
- effort in fisheries
  - average HR

# oem
  - CALL pdynlfcpue
  - ADD stock OE
  -

* SEND RH SWO example oem, om, mp + EOP today

* LOAD refpts

# MPS
  - cpue.ind + buffer.hcr?
