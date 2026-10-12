
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evalpheno <img src="fig/evalpheno.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15174551.svg)](https://doi.org/10.5281/zenodo.15174551)
<!-- badges: end -->

`evalpheno` lets you construct temperature-driven phenology models, set
their parameters, and simulate variation among individual buds. Models
combine chill and heat accumulation with sequential, parallel,
partial-overlap or PhenoFlex structures.

## Installation

Install the development version from GitHub:

``` r
install.packages("remotes")
remotes::install_github("larscaspersen/eval_phenoflex")
```

## Construct a model

Choose a model structure and the chill and heat submodels:

``` r
library(evalpheno)

model <- pheno_model(
  structure = "phenoflex",
  chill = chill_dynamic("characteristic"),
  heat = heat_gdh("unscaled")
)
```

Use `"sequential"`, `"parallel"` or `"partial_overlap"` for another
structure, or `chill_dynamic("kinetic")` for kinetic chill parameters.

## Parameterize and predict

Set parameters by name; omitted values use the selected model’s
defaults. `parameter_schema(model)` describes its parameters and
`model_parameters(model)` returns their current values.

``` r
model <- pheno_model(parameters = c(yc = 40, zc = 250))

# Prepare two seasons from the bundled daily weather data.
data("cka_weather", package = "evalpheno")
hourly <- chillR::stack_hourly_temps(
  cka_weather, latitude = attr(cka_weather, "latitude")
)$hourtemps
seasons <- chillR::genSeasonList(hourly, mrange = c(8, 6), years = 2008:2009)

predict_phenology(model, seasons)
#> [1] 121.0417 112.3750
```

Predictions are fractional days of the year; `NA` means the event was
not reached. Use `basic_output = FALSE` to inspect chill and heat
trajectories.

## Model a bud population

Share chill and potential heat accumulation across buds while allowing
their requirements to vary:

``` r
population <- population_pheno_model(
  structure = "phenoflex", n = 100,
  sd = c(yc = 2, zc = 10), seed = 17
)
means <- model_parameters(population)
means[c("yc", "zc")] <- c(40, 250)

buds <- predict_population_phenology(population, seasons[[1]], parameters = means)
quantile(buds$bloom_jday, c(0.1, 0.5, 0.9), na.rm = TRUE)
#>      10%      50%      90%
#> 119.1125 120.8750 123.0833
```

Population models also support cutting and constant-temperature forcing
experiments.

For full workflows, see [Getting
started](vignettes/getting-started.Rmd), [Model
calibration](vignettes/model-calibration_special-cases.Rmd), and [Bud
populations](vignettes/bud-population.Rmd).
