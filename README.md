
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

## Population model

The population model is available directly from the installed package:

``` r
library(evalpheno)
weather <- data.frame(Temp = rep(15, 240), JDay = rep(1:10, each = 24))
par <- c(yc = 40, zc = 190, s1 = 0.5, Tu = 25,
         E0 = 3372.8, E1 = 9900.3, A0 = 6319.5, A1 = 5.939917e13,
         Tf = 4, Tc = 36, Tb = 4, slope = 1.6)
result <- phenoflex_population(weather, par, n = 100,
                              yc_sd = 2, zc_sd = 10, seed = 12345)
```

Use `jday_cut` to simulate forcing experiments and `population` to reuse an
explicit set of bud requirements across runs. The existing experimental call
`helper_run_pop_model(...)` is also available after `library(evalpheno)`.
`sourceCpp()` and sourcing files from `experimental/` are unnecessary for these
package functions. Normal sampling uses base R; skew-normal sampling requires `sn`.

Read `?phenoflex_population` for parameter order, hourly input requirements,
matrix dimensions, and the retained index conventions. `bloomindex = 0`
indicates no bloom; forcing output `exp` retains zero-based step indices and
the failure sentinel 9999. See `NEWS.md` for intentional behavior changes.

The `experimental/` directory remains a research workspace and is excluded
from package builds. Study-specific plots, data readers, and stage classifiers
are not yet part of the public package API.

### Parameter terminology

Use `phenoflex_parnames_characteristic` for theta_star, theta_c, tau and pie_c,
and `phenoflex_parnames_kinetic` for E0, E1, A0 and A1 (positions 5:8).
These are two parameterizations of the same chill model. Both vectors contain
all 12 names in model input order; conversions preserve the other eight values.

```r
kinetic <- characteristic_to_kinetic(characteristic)
characteristic <- kinetic_to_characteristic(kinetic)
```

The old/new name vectors and `convert_parameters()` /
`convert_parameters_old_to_new()` remain available as compatibility aliases.
The inverse conversion retains its existing numerical algorithm and can fail
for parameter sets outside its supported domain.

### Named sequential model interface

```r
model <- pheno_model(
  structure = "sequential",
  chill = chill_dynamic("characteristic"),
  heat = heat_gdh()
)
parameter_schema(model)
parameters <- default_parameters(model)
parameters["yc"] <- 45
prediction <- predict_phenology(model, weather = season, parameters = parameters)
```

`model` is an S3 list describing the algorithms. The named numeric parameter
vector is separate and may be reordered without changing the prediction.
`season` must contain complete consecutive hourly days with Temp, Year and JDay;
optional Hour must run from 0 to 23 each day. The result is a list with a one-based
`bloomindex` (zero means no bloom); `basic_output = FALSE` also returns `chill`
and `z`. Use `return_JDay(prediction$bloomindex, season$JDay, season$Year)`
for the fractional calendar result. Sequential, linear parallel, partial overlap
and PhenoFlex structures share Dynamic chill and scaled or unscaled GDH modules.

`chill_dynamic("kinetic")` selects E0/E1/A0/A1 instead. Kinetic defaults use the
original coefficients; characteristic defaults are a separate starting set and
are not their conversion. See `development/sequential_model_example.R` for a
runnable example with frozen station weather and equivalent parameter sets.

## Modular populations and forcing

```r
model <- population_pheno_model(
  structure = "phenoflex",
  chill = chill_dynamic("characteristic"),
  heat = heat_gdh("unscaled"),
  n = 100,
  sd = c(yc = 2, zc = 10),
  seed = 17
)
parameters <- default_parameters(model)
buds <- sample_population_parameters(model, parameters)
prediction <- predict_population_phenology(
  model, season, parameters,
  population = buds,
  cut_indices = c(240, 480),
  forcing_temperature = 23,
  max_hours_forcing = 1200
)
prediction$bloom_jday
prediction$forcing$hours_to_bloom
```

Use a season with at least 480 hourly rows for these example cuts. The population
specification holds its underlying single model in `model$model`. Named parameters
are population means; structure traits can vary between buds while chill and
potential heat are calculated once. For partial overlap, specify dispersion of
`b1` and `b2` explicitly instead of `zc`. Convert heat standard deviations as
well as means when changing GDH scaling. Zero dispersion reproduces a single bud.
Reuse an explicit buds-by-structure-parameters matrix for deterministic fitting
or correlated traits. Sampled parameters outside their model domains raise errors.

Cut indices identify one-based weather rows. Forcing retains accumulated field
heat and structure history while holding all chill pools at their cutting values.
Results are elapsed hours: zero means already bloomed at cutting and NA means
the requirement was not reached. Detailed output adds hours-by-buds heat matrices.
These units differ from the legacy `phenoflex_population()` forcing step indices;
that interface remains available with its original behavior.

