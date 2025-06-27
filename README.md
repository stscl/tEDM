
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tEDM <img src="man/figures/logo.png" align="right" height="139" alt="https://stscl.github.io/tEDM/">

<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/tEDM)](https://CRAN.R-project.org/package=tEDM)
[![CRAN
Release](https://www.r-pkg.org/badges/last-release/tEDM)](https://CRAN.R-project.org/package=tEDM)
<!-- [![CRAN Checks](https://badges.cranchecks.info/worst/tEDM.svg)](https://cran.r-project.org/web/checks/check_results_tEDM.html) -->
[![Downloads_all](https://badgen.net/cran/dt/tEDM?color=orange)](https://CRAN.R-project.org/package=tEDM)
[![Downloads_month](https://cranlogs.r-pkg.org/badges/tEDM)](https://CRAN.R-project.org/package=tEDM)
[![License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/stscl/tEDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stscl/tEDM/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-20b2aa.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-universe](https://stscl.r-universe.dev/badges/tEDM?color=cyan)](https://stscl.r-universe.dev/tEDM)

<!-- badges: end -->

***T**emporal **E**mpirical **D**ynamic **M**odeling*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=tEDM) with:

``` r
install.packages("tEDM", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/tEDM) with:

``` r
install.packages("tEDM",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/tEDM)
  with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/tEDM",
                         #build_vignettes = TRUE,
                         dep = TRUE)
```
