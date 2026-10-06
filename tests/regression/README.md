# Weather regression baselines

These scripts record the existing implementations before modular refactoring.
The expected outputs have already been generated and saved as CSV fixtures
(include these files when committing the work). Ordinary regression runs execute
only the candidate implementation and compare it with the saved values. They
never execute a second copy of the old implementation or update the baseline.

## Run after a code change

From the package root in R, install the candidate into the local library first:

```r
dir.create(".local-library", showWarnings = FALSE)
install.packages(".", repos = NULL, type = "source", lib = ".local-library")
```

Then **restart R** so an already loaded namespace or DLL cannot mask the new build:

```r
.libPaths(c(normalizePath(".local-library"), .libPaths()))
source("tests/regression/run_regression.R")
```

On PowerShell, with Rscript on PATH, the equivalent comparison command is:

```powershell
$env:R_LIBS = Join-Path (Get-Location) '.local-library'
Rscript tests/regression/run_regression.R
```

The runner prints the package library it actually tests. It writes a comparison
to `development/regression-results/comparison.csv` and fails on differences.
Numerical tolerance is `1e-8 + 1e-10 * abs(expected)`. Missing values and output
shapes must match. The test is also included in the usual testthat/package checks.

## Saved results

Under `tests/testthat/fixtures/weather-regression/`:

- `baseline.csv`: all 3,858 expected values.
- `baseline_cka.csv`: Klein-Altendorf outputs.
- `baseline_quillota.csv`: Quillota outputs.
- `baseline_chillR_KA.csv`: outputs for two chillR example seasons.
- `inputs.rds`: frozen hourly seasons, complete parameter vectors and synthetic observations.
- `seasons.csv`, `parameters.csv`, `synthetic_observations.csv`: readable input summaries.
- `weather_sources.csv`, `baseline_sources.csv`, `session.txt`: source hashes and versions.

Rows identify the station, season, function, scenario, cultivar, stage, output
and vector element. Filter `output == "bloom_jday"` for predictions, `"F"` or
`"rss"` for objectives and `"g"` for constraint vectors. Daily chill/heat rows
contain end-of-day cumulative values in chronological order starting at the
season's start date. NA denotes no bloom, while evaluation functions may instead
return the existing 365-day penalty. Negative bloom days denote the preceding
calendar year under the current wrapper convention.

## Coverage and choices

| Requested scenario | Baseline cases |
| --- | --- |
| Successful bloom | `ordinary` for all six seasons |
| Crossing 31 December | Northern August-June seasons, `calendar_dec31` and `calendar_jan01` |
| Leap year | 2007/08 northern season contains 29 February; `calendar_feb29` tests actual hourly rows |
| No bloom | `no_bloom` and `no_bloom_penalty` |
| Two cultivars | `combined`, separate A/B predictions, summed objective and constraints |
| Multiple stages | Three heat thresholds, detailed trajectories and `three_stages` objective/constraints |
| Original/intermediate Dynamic Model parameters | `original_dynamic`, `ordinary` frozen kinetic equivalent, `characteristic` and `characteristic_evaluation` |
| Sequential, parallel, partial overlap | Each wrapper and its `eval_all_daoptim` objective |
| Fixed parameters | `fixed_parameters` and `fixed_sequential` |

The uploaded daily data cover 2007-2009. Klein-Altendorf uses August 2007-June
2008 and August 2008-June 2009. Quillota uses April-November 2008 and 2009,
reflecting the southern hemisphere. Two further seasons use chillR's KA_weather
for the northern periods. Its outputs currently match the uploaded CKA data;
these are separate provenance cases, not independent climates.

Hourly temperatures are generated once with chillR::stack_hourly_temps and the
supplied latitude. Longitude is recorded in the source coordinates but is not an
argument of that interpolator. The preparation script rejects missing daily
temperatures, date gaps, duplicate dates and Tmin > Tmax. Northern seasons use
genSeasonList; it requires a cross-year range, so southern seasons are selected
explicitly. Normal tests do not repeat weather interpolation.

Parameter sets and cultivar labels are test inputs, not calibrated cultivar
estimates. Evaluation observations are fixed synthetic dates, not measured bloom
data, and are independent of generated model predictions. Parameters intentionally
exercise different outputs; they need not produce agronomically realistic bloom
for every model. This suite does not yet exercise exchangeable Chill Hours or
Positive Utah submodels, optimizer searches, or a new sequential-stage fitting
API. Those need additional cases as the modular interfaces are introduced.

## Corrected three-stage problem

The evaluator previously indexed an hourly matrix with a fractional/possibly
negative Julian day instead of the hourly threshold-crossing row. This is now
fixed: the chill check uses the first hourly budburst index, and a missing
budburst yields the existing penalty safely.

Independent correctness tests in `tests/testthat/test-three-stage-index.R`
cover December/January, leap day, insufficient chill, missing budburst, and
agreement with separately evaluated stage dates for all six weather seasons.
Only five baseline values changed: three stage-objective values and the two
historical error indicators, now zero. All other values and constraints stayed
unchanged. `before-three-stage-fix/` preserves the old baseline, provenance and
`reviewed_changes.csv`. The historical case name `known_three_stage_index_error`
remains so the correction is directly comparable to the original snapshot.
## Sequential timing correction (28 September 2026)

Sequential heat now uses the chill state after the current interval's update,
matching the parallel kernel's timing. The interval reaching the chill requirement
contributes its full heat increment. An independent analytic test in
`tests/testthat/test-sequential-timing.R` covers threshold crossing, equality,
the final interval, insufficient chill and temperatures below the heat base.

Only two saved values changed: the `fixed_sequential` objectives for CKA and
chillR_KA, each from approximately 69401.25 to 69432.29. All other saved values remain
within the existing regression tolerance. The previous baselines and provenance
are preserved under `before-seq-timing/`, with the reviewed differences.

## Reviewed partial-overlap correction and modular implementation

The partial-overlap heat requirement now uses accumulated chill y, as corrected
by the user, instead of the labile x pool. The 12 resulting prediction/objective
changes across the six station seasons were reviewed during the modular migration.
They match the differences recorded before the modular change. The pre-refactor
and modular parallel/partial-overlap implementations produced identical x/y/z
trajectories and bloom indices in all 12 station/model comparisons.

`before-po-chill/` preserves the preceding baselines, provenance,
reviewed changes and the modular comparison. All other baseline values are unchanged.
Independent tests now check last-column selection, the use of frozen y in the
heat requirement and the intended freeze of residual x after the stopping point.

## Deliberately creating a new baseline

The additional-chill correction to partial overlap was checked against an
independent vector-based R calculation for all six station seasons. It changes
only the six partial-overlap predictions and their six evaluation scores. The
previous snapshots, reviewed differences and reference comparison are preserved
in `before-additional-chill/`. The first row reaching yc supplies the baseline
for Ca; hourly threshold overshoot is included in that baseline. Analytic tests
also cover zero additional chill, freezing, initial completion and invariance to
a shift of the total-chill origin. Heat-onset and overlap-stop timing were not
changed by this correction.

Only do this using a reviewed implementation that should define the new reference:

```powershell
# Initial generation (already completed here):
Rscript tests/regression/prepare_inputs.R
Rscript tests/regression/generate_baseline.R

# Deliberate replacement; these commands otherwise refuse to overwrite files:
Rscript tests/regression/prepare_inputs.R --overwrite
Rscript tests/regression/generate_baseline.R --overwrite
```

Keep the frozen weather unless inputs intentionally change. For a reviewed model
fix, only regenerate the output baseline, inspect the CSV differences, and retain
the old snapshot in version control. Never regenerate expected results merely
to make a failing regression test pass.

