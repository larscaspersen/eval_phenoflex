
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

model <- pheno_model()
```

By default, `pheno_model()` constructs a PhenoFlex model with the
Dynamic Model for chill accumulation and the Growing Degree Hour Model
for heat accumulation. You can select other model structures using
`"sequential"`, `"parallel"` or `"partial_overlap"`.

You can visualize the model structure:

``` r
plot(model)
```

![](README_files/figure-gfm/plot-structure-1.png)<!-- -->

Use `parameter_schema()` to view the model parameters and their
descriptions:

``` r
parameter_schema(model)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

name
</th>

<th style="text-align:left;">

component
</th>

<th style="text-align:right;">

default
</th>

<th style="text-align:left;">

description
</th>

<th style="text-align:left;">

unit
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

yc
</td>

<td style="text-align:left;">

structure
</td>

<td style="text-align:right;">

40
</td>

<td style="text-align:left;">

Chill requirement
</td>

<td style="text-align:left;">

Chill Portions
</td>

</tr>

<tr>

<td style="text-align:left;">

zc
</td>

<td style="text-align:left;">

structure
</td>

<td style="text-align:right;">

285.71
</td>

<td style="text-align:left;">

Heat requirement
</td>

<td style="text-align:left;">

h
</td>

</tr>

<tr>

<td style="text-align:left;">

s1
</td>

<td style="text-align:left;">

structure
</td>

<td style="text-align:right;">

0.5
</td>

<td style="text-align:left;">

Slope parameter for transition from chill accumulation to heat
accumulation
</td>

<td style="text-align:left;">

1 / (Chill Portions)
</td>

</tr>

<tr>

<td style="text-align:left;">

E0
</td>

<td style="text-align:left;">

chill
</td>

<td style="text-align:right;">

4153.5
</td>

<td style="text-align:left;">

time-independent activation energy of forming the precursor to the
dormancy-breaking factor
</td>

<td style="text-align:left;">

dimensionless
</td>

</tr>

<tr>

<td style="text-align:left;">

E1
</td>

<td style="text-align:left;">

chill
</td>

<td style="text-align:right;">

12888.8
</td>

<td style="text-align:left;">

time-independent activation energy of destroying the precursor to the
dormancy-breaking factor
</td>

<td style="text-align:left;">

dimensionless
</td>

</tr>

<tr>

<td style="text-align:left;">

A0
</td>

<td style="text-align:left;">

chill
</td>

<td style="text-align:right;">

139500
</td>

<td style="text-align:left;">

amplitude of the (hypothetical) process involved in forming the
precursor to the dormancy-breaking factor
</td>

<td style="text-align:left;">

dimensionless
</td>

</tr>

<tr>

<td style="text-align:left;">

A1
</td>

<td style="text-align:left;">

chill
</td>

<td style="text-align:right;">

2567000000000000000
</td>

<td style="text-align:left;">

amplitude of the (hypothetical) process involved in destroying the
precursor to the dormancy-breaking factor
</td>

<td style="text-align:left;">

dimensionless
</td>

</tr>

<tr>

<td style="text-align:left;">

Tf
</td>

<td style="text-align:left;">

chill
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:left;">

Transition temperature for irreversible chill conversion
</td>

<td style="text-align:left;">

degree C
</td>

</tr>

<tr>

<td style="text-align:left;">

slope
</td>

<td style="text-align:left;">

chill
</td>

<td style="text-align:right;">

1.6
</td>

<td style="text-align:left;">

Slope of temperature-dependent chill conversion
</td>

<td style="text-align:left;">

K^-1
</td>

</tr>

<tr>

<td style="text-align:left;">

Tb
</td>

<td style="text-align:left;">

heat
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:left;">

Base temperature for heat accumulation
</td>

<td style="text-align:left;">

degree C
</td>

</tr>

<tr>

<td style="text-align:left;">

Tu
</td>

<td style="text-align:left;">

heat
</td>

<td style="text-align:right;">

25
</td>

<td style="text-align:left;">

Optimal temperature for heat accumulation
</td>

<td style="text-align:left;">

degree C
</td>

</tr>

<tr>

<td style="text-align:left;">

Tc
</td>

<td style="text-align:left;">

heat
</td>

<td style="text-align:right;">

36
</td>

<td style="text-align:left;">

Upper critical temperature for heat accumulation
</td>

<td style="text-align:left;">

degree C
</td>

</tr>

</tbody>

</table>

Supply model parameters as a **named** vector when constructing the
model, or adjust them afterward:

``` r
model$parameters['yc'] <- 80

print(model)
#> Structure:      PhenoFlex
#> Chill submodel: Dynamic Model
#> Heat submodel:  GDH
#> Parameters:
#>                  yc                  zc                  s1                  E0
#>                  80              285.71                 0.5              4153.5
#>                  E1                  A0                  A1                  Tf
#>             12888.8              139500 2567000000000000000                   4
#>               slope                  Tb                  Tu                  Tc
#>                 1.6                   4                  25                  36
```

Alternative chill submodels are `chill_hours()`, `chill_utah()` and
`chill_positive_utah()`. `chill_none()` creates a heat-only model,
bypassing chill-dependent coupling. Chill-hours thresholds and Utah zone
weights appear in the parameter schema and can be calibrated.
User-defined chill submodels can be supplied through `chill_custom()`.

## Predict bloom dates

Prepare hourly temperature time series to predict bloom dates. The
[chillR
package](https://cran.r-project.org/web/packages/chillR/index.html)
(Luedeling et al., 2024) provides tools for acquiring and processing
temperature data:

``` r
data("KA_weather", package = "chillR")
hourly <- chillR::stack_hourly_temps(
  KA_weather, latitude = 50.4)$hourtemps
seasons <- chillR::genSeasonList(hourly, mrange = c(8, 6), years = 1999:2009)
names(seasons) <- as.character(1999:2009)
```

``` r
model <- pheno_model(parameters = c(yc = 40, zc = 250))

predict_phenology(model, seasons)
#>      1999      2000      2001      2002      2003      2004      2005      2006 
#> 111.16667 113.58333 112.91667  92.04167 103.33333 103.04167 102.87500 121.87500 
#>      2007      2008      2009 
#> 102.29167 117.95833 109.08333
```

## Visualize the processes

Plot chill and heat accumulation and the transition between them:

``` r
plot_processes(model, seasons, ribbon = TRUE)
```

![](README_files/figure-gfm/model-plots-1.png)<!-- -->

`plot_processes()` shows chill accumulation, the transition, and
effective heat over time. With a season list, grey curves show
individual seasons and black curves show pointwise medians. A dot marks
the median bloom date, with whiskers spanning the earliest and latest
successful predictions. Set `ribbon = TRUE` to replace the individual
curves with a shaded range from the minimum to the maximum.
Partial-overlap models show their changing heat requirement in the
transition panel; heat-only models show just heat accumulation.

## Calibrate model parameters

Use `fit_phenology()` to calibrate the model parameters. Only the
parameters specified in the lower and upper bounds are adjusted. Install
the [DEoptim
package](https://cran.r-project.org/web/packages/DEoptim/index.html)
Ardia et al. (2011) to use the default optimizer. Alternatively, install
the [GenSA
package](https://cran.r-project.org/web/packages/GenSA/index.html)
(Xiang et al., 2013) and supply `optimizer = "GenSA"` in the function
call.

``` r
require(DEoptim)
#> Loading required package: DEoptim
#> Loading required package: parallel
#> 
#> DEoptim package
#> Differential Evolution algorithm in R
#> Authors: D. Ardia, K. Mullen, B. Peterson and J. Ulrich
model <- pheno_model()
lower <- c(yc = 10)
upper <- c(yc = 80)

# Replace these synthetic predictions with observed bloom dates.
observed <- predict_phenology(model, seasons)

# Run the optimizer. Use more iterations for actual calibrations.
fit <- fit_phenology(model, evaluation = phenology_rss, observed = observed,
                     seasons = seasons, lower, upper, seed = 17, max_iterations = 2)
#> Iteration: 1 bestvalit: 2.612847 bestmemit:   38.552002
#> Iteration: 2 bestvalit: 2.612847 bestmemit:   38.552002

print(fit)
#> Calibrated phenology model (DEoptim)
#> Loss: 0 
#> Estimated parameters:
#> yc 
#> 40
```

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
#> 102.2417 111.0417 113.0458
```

Population models also support cutting and constant-temperature forcing
experiments.

Read the [online guide](https://larscaspersen.github.io/eval_phenoflex/)
for full workflows, with chapters on [Getting
started](https://larscaspersen.github.io/eval_phenoflex/getting-started.html),
[Model structures and
submodels](https://larscaspersen.github.io/eval_phenoflex/model-structures-and-submodels.html),
[Model
calibration](https://larscaspersen.github.io/eval_phenoflex/model-calibration_special-cases.html),
and [Bud
populations](https://larscaspersen.github.io/eval_phenoflex/bud-population.html).

## References

<div id="refs" class="references csl-bib-body hanging-indent"
data-entry-spacing="0">

<div id="ref-ardia_differential_2011" class="csl-entry">

Ardia, D., Boudt, K., Carl, P., Mullen, K.M., Peterson, B.G., 2011.
[Differential Evolution with DEoptim: An application to non-convex
portfolio optimization](https://doi.org/10.32614/RJ-2011-005). R Journal
3, 27–34.

</div>

<div id="ref-luedelingChillRStatisticalMethods2024" class="csl-entry">

Luedeling, E., Caspersen, L., Fernandez, E., 2024.
<span class="nocase">chillR</span>: Statistical methods for phenology
analysis in temperate fruit trees.
<https://CRAN.R-project.org/package=chillR>

</div>

<div id="ref-deoptim_package_2011" class="csl-entry">

Mullen, K., Ardia, D., Gil, D., Windover, D., Cline, J., 2011. [DEoptim:
An R package for global optimization by Differential
Evolution](https://doi.org/10.18637/jss.v040.i06). Journal of
Statistical Software 40, 1–26.

</div>

<div id="ref-xiang_generalized_2013" class="csl-entry">

Xiang, Y., Gubian, S., Suomela, B., Hoeng, J., 2013. Generalized
simulated annealing for efficient GlobalOptimization: The GenSA package
for r. The R Journal 5/1. <https://journal.r-project.org>

</div>

</div>
