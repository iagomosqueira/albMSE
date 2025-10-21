---
title:
author: "Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>"
tags:
---


# How to install

```
install.packages("TAF")
pkgs <-  TFA::deps()

install.packages(c(pkgs), repos = c(
    CRAN = "https://cloud.r-project.org/",
    SPFRMO = "https://sprfmo.r-universe.dev/"))
```

# Structure

- boot
  - initial/data: SS3 inputs
  - data.R: RIUN ss3 base
   - data: SS3 model run
- data.R: PREPARE inputs for abc
- model_.R: RUN ABC
- output.R: LOAD base SS3 run & ABC results
