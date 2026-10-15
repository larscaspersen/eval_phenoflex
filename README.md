
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
  chill = chill_dynamic(),
  heat = heat_gdh("unscaled")
)
```

Use `"sequential"`, `"parallel"` or `"partial_overlap"` for another
structure, with `chill_dynamic()` for Dynamic Model chill.
Models use the original `E0`, `E1`, `A0` and `A1` coefficients by default.
Supplying `theta_star`, `theta_c`, `tau` or `pie_c` automatically selects
characteristic input; mixing it with kinetic input raises an error. Characteristic
values are converted at construction and retained in `model$chill$characteristic`;
`model$parameters` always contains kinetic values.
Alternative chill submodels are `chill_hours()`, `chill_utah()` and
`chill_positive_utah()`. `chill_none()` creates a heat-only model, bypassing
chill-dependent coupling. Chill-hours thresholds and Utah zone weights appear
in the parameter schema and can be calibrated. Positive Utah adds only positive
daily totals. Set chill requirements and search bounds in the selected model's
own units; they are not interchangeable with Dynamic Model portions.
Default `yc` is 1000 for chill hours, 1200 for Utah and positive Utah, and
40 for the Dynamic Model.

User-defined chill submodels can be supplied through `chill_custom()`:

```r
triangular_chill <- function(weather, parameters) {
  temp <- head(weather$Temp, -1L)
  increments <- pmax(0, 1 - abs(temp - parameters[["Topt"]]) / parameters[["width"]])
  matrix(c(0, cumsum(increments)), ncol = 1, dimnames = list(NULL, "y"))
}
custom <- chill_custom(
  triangular_chill, defaults = c(Topt = 5, width = 5),
  lower = c(Topt = 0, width = 1), upper = c(Topt = 10, width = 15),
  validate = function(p) p[["width"]] > 0, label = "Triangular chill"
)
custom_model <- pheno_model(chill = custom)
```

Construction checks the callback on a small hourly sample; prediction checks
every returned trajectory. `validate_chill_custom(custom, weather)` checks
supplied weather explicitly. See `help("chill_custom")` for the matrix contract,
optional sample weather and parameter-domain validation.

## Parameterize and predict

Set parameters by name; omitted values use the selected model’s
defaults. `parameter_schema(model)` describes its parameters and
`model_parameters(model)` returns their current values.

For characteristic calibration, request bounds with
`default_bounds(model, parameterization = "characteristic")`, or name selected
characteristic parameters in the bounds. The fitter infers the representation
from the bounds or explicit starting values. `model_parameters(model,
"characteristic")` supplies equivalent characteristic starting values. Fitted
models store kinetic coefficients; `fit$par` retains the optimizer's coordinates.

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

## Visualize the model

Plot the chill and heat temperature responses together with their transition:

``` r
model <- pheno_model("partial_overlap")
plot(model)
plot_temperature_response(model)
plot_model_structure(model)

# Overlay the season's effective chill/heat path and actual overlap cutoff.
plot_model_structure(model, weather = seasons[[1]])
plot_processes(model, seasons[[1]])
plot_processes(model, seasons)
plot_processes(model, seasons, ribbon = TRUE)
```

The partial-overlap plot marks minimum chilling `yc` and the initial heat
requirement `b1 + b2`, with the lower heat asymptote `b1` and overlap cutoff `b1 * ol`.
The region below minimum chilling is shaded: no amount of heat can trigger the event.
Dynamic Model chill axes stop at 120 portions by default; `chill_range` overrides this.
The latter is an accumulated-heat threshold; its chill coordinate depends on
the weather. Other structures show the chill-dependent fraction of effective
heat. Temperature responses use the native submodels and chosen heat scaling;
Dynamic Model chill is averaged over 1200 hours at constant temperature.

The component functions return editable ggplot objects, and `plot(model)`
returns a patchwork composition. Install `ggplot2` and `patchwork` to use them.
Observed-temperature histograms or gradients can be added with `hourtemps` and
`weather_freq_plot` in `plot_temperature_response()` or `plot_pheno_model()`.
In `plot(model, weather = seasons[[1]])`, the supplied season is used automatically.
Panel labels show the selected `chill_months` and `heat_months`, and the caption
identifies the weather source and relative frequencies.

`plot_processes()` follows chill accumulation, the transition, and effective heat
through time. With a season list, grey curves show individual seasons and black
curves show pointwise medians. A dot marks the median bloom date, with whiskers
spanning the earliest and latest successful predictions. Set `ribbon = TRUE`
to replace the individual curves with a shaded min-max range. Partial overlap shows
its changing heat requirement in the transition panel; heat-only models show
just heat accumulation.

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
started](vignettes/getting-started.Rmd),
[Model structures and submodels](vignettes/model-structures-and-submodels.Rmd), [Model
calibration](vignettes/model-calibration_special-cases.Rmd), and [Bud
populations](vignettes/bud-population.Rmd).
