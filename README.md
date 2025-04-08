
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evalpheno <img src="fig/evalpheno.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15174551.svg)](https://doi.org/10.5281/zenodo.15174551)
<!--[![CRAN status](https://www.r-pkg.org/badges/version/hexsession)](https://CRAN.R-project.org/package=hexsession)>
<!-- badges: end -->

`evalpheno` is a collection of evaluation functions and wrapper
functions that allow customized calibration phenology models coming from
the PhenoFlex modeling framework.

`evalpheno` aims to expand the calibration of the phenology model
PhenoFlex, part of the chillR package. By customizing evaluation
functions or wrapper functions input parameters can be fixed or
parameters can be replaxed by more narrowly defined intermediate
parameters. The evaluation functions make it easier to calibrate the
models with other global optimization algorithms. Also more structural
changes could be made to the model, like sharing chill and heat
accumulation submodel parameters across cultivars of the same species,
while still having cultivar-specific chill and heat requirements and
transition parameters. Or several phenological stages could be evaluated
in one model, instead of having seperate models for each stage.

## Installation

You can install the development version of evalpheno like so:

``` r
install.packages('devtools')
devtools::install_github('https://github.com/larscaspersen/eval_phenoflex')
```

## 
