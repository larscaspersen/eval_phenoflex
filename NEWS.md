# evalpheno 0.0.2.0

* Support parallel DEoptim fitting through native parallel controls and supplied
  clusters. Use the optimizer's returned optimum and evaluation count, preserving
  a better initial baseline. Load worker packages and stop fitter-owned clusters
  on both success and errors.

* Compact CV preparation to data, assignments and settings, deriving split
  indices from the inspected assignments. Compact fitted CV results to model,
  predictions and calibration, storing fold models once and retaining only
  independent validation weather. Add brief print/summary methods and optional
  keep_diagnostics and keep_training_predictions; scores and the optional refit
  now live under calibration.

* Accept single and CV fitter results directly in predict_phenology(). CV results
  contain a configured fold-model ensemble; an optional full-data refit remains
  separately selectable. Single fits group options under calibration_settings and
  run statistics under diagnostics. Store paired predicted/observed calibration
  dates for supplied RSS losses, with a prediction callback for custom losses.
  CV results retain held-out pairs and optionally training pairs.
  validate_phenology(cv) uses retained validation data without automatic holdout scoring.

* Add prepare_phenology_cv() with inspectable season-entry fold assignments,
  repeated CV and independent validation selected by original entry indices
  or a proportion of observed entries. Years are inspection metadata, allowing
  the same year at different locations to be withheld separately. Optional
  groups keeps related entries together within CV folds; stage seasons stay
  together automatically. fit_phenology_cv()
  consumes prepared splits unchanged, selects optimizer restarts by training loss,
  retains fold models and out-of-fold dates, and optionally refits all CV entries.
  Support single, combined cultivar and stage layouts, preserving NULL observations.
  Evaluate reserved entries separately with validate_phenology().

* Add pheno_ensemble(), predict_phenology_ensemble() and
  summarise_ensemble_predictions(). Aggregate independent fits with equal,
  explicit or inverse assessment-MSE weights, count/weight no-bloom votes,
  weighted descriptive spread and iterative weight caps. Keep cultivar/stage
  outputs separate. predict_phenology() accepts ensembles with stored parameters.

* Expose observed as an explicit named argument of fit_phenology(). Pass supplied
  observations unchanged to each evaluation, retaining modular RSS/custom losses
  and evaluators that capture observations when the argument is omitted.

* Accept named theta_star and theta_c values in 0--20 as Celsius using the
  native kernels' legacy +273 convention; retain Kelvin inputs. Normalize model
  construction, prediction overrides, calibration starts and bounds, including
  indexed collection parameters. Evaluators and fitted results use Kelvin.
  Add temperature_unit = "C" to default_bounds() for Celsius display.

* Add default_bounds() with model-specific suggested calibration ranges,
  reference conversion for scaled GDH and indexed bounds for cultivar/stage
  collections. fit_phenology() uses these when bounds are NULL; a single omitted
  bound fills only the explicitly selected parameters. Explicit partial bounds
  still fix all remaining values. Exclude numerical characteristic conversion
  failures with infinite loss while continuing to propagate other evaluator errors.

* Extend predict_phenology() to accept seasonlists for a single model and
  ordinary lists of independent models, including differing specifications.
  Model collections accept one common weather data frame, a common seasonlist,
  or nested seasonlists matched by model position. Preserve model and season
  labels and support detailed output and per-model parameter overrides.

* Allow named parameter subsets in pheno_model(). Omitted values use defaults
  for the selected model specifications; validate the complete resulting vector.

* Add pheno_model_list() to calibrate collections of ordinary simple models,
  stage_pheno_models() for cumulative ordered heat thresholds of one cultivar,
  and the explicitly supplied phenology_rss_stages() evaluator. NULL observation
  slots skip only their own stage/season pair. Only observed pairs are predicted.
  Add model_parameters() and as_pheno_models() to inspect stored calibration
  values and extract independent single models from collections, combined
  specifications and fitted results. The combined wrapper is now optional.

* Add combined_pheno_model() for cultivar-specific subsets of structure
  parameters and shared remaining parameters. Add cultivar_parameters(),
  predict_combined_phenology() and the explicitly supplied phenology_rss_combined()
  evaluator for nested cultivar seasonlists and observation vectors. Both
  optimizers support indexed bounds and preserve omitted parameters.

* Allow partial named bounds in fit_phenology(): parameters omitted from both
  bounds remain constant at their initial model values. Evaluators receive
  complete parameter vectors; lower and upper must name the same subset.

* Add phenology_rss() as a standalone example evaluation function: predict each
  season with candidate parameters and return RSS against observed dates. Pass
  the evaluator explicitly to fit_phenology().

* Add fit_phenology() for calibration with DEoptim or GenSA, a shared evaluation
  interface, named bounds (including fixed parameters), local seeds, iteration
  limits, a native GenSA time limit, and a fitted model specification. DEoptim
  time limits are unsupported; elapsed time is measured for reporting only.

* Add population_pheno_model(), sample_population_parameters() and
  predict_population_phenology() for modular bud populations with shared
  Dynamic chill/GDH inputs and bud-specific structure parameters. Support
  explicit populations, independent normal/skew-normal draws and local seeds.
  Constant-temperature forcing holds chill at cutting and retains field heat
  and structure history. Forcing reports elapsed hours (zero for already met,
  NA for not reached), separately from the legacy population kernel conventions.

* Add calculate_heat_gdh_unscaled() and apply_phenoflex_structure() for modular
  PhenoFlex coupling with interchangeable chill and heat inputs. The structure
  preserves PhenoFlex's start-of-interval chill timing and sigmoid limits.

* Partial overlap now returns continuing chill pools after overlap ends; only
  additional chill used in the heat requirement is held constant. Bloom
  predictions and zero-overlap behavior are preserved.

* Fix ol=0 in partial overlap: complete initial chilling,
  then accumulate heat toward b1+b2 with zero additional chill.

* Correct partial-overlap compensation to use additional chill after the first
  yc crossing, following Pope et al. (2014). The crossing row is the baseline,
  including hourly overshoot; total chill pools retain their existing values.
  b2 now represents the extra heat requirement at zero additional chill.

* Remove chill freezing from sequential, linear parallel and Landsberg structures.
  All pools continue updating; partial overlap holds only the coupling value.
  Landsberg's yc argument is retained for compatibility but no longer affects
  the calculation. This supersedes the earlier freezing descriptions below.

* Document linear parallel competence with Hänninen & Kramer (2007, Eq. B4b)
  and correct the meaning of kmin. Add apply_parallel_landsberg_structure()
  with exponential effectiveness from Landsberg (1974, Eq. 5), an explicit
  chill scale y0, and freezing of all chill pools at yc.

* Add apply_parallel_structure() and apply_partial_overlap_structure() using the
  shared precomputed chill and heat outputs. parallel_model() and po_model()
  now compose these modules. Both freeze all chill pools at their structure's
  stopping condition, with no subsequent conversion of residual x.

* Split the sequential model into calculate_chill_dynamic(), calculate_heat_gdh()
  and apply_sequential_structure(), implemented in separate C++ files. Chill
  outputs retain x/xs/y names with y last; the structure uses the last column.
  seq_model() now composes these modules and retains its calling/return interface.
  Once yc is reached, x and y are frozen, including any unconverted x remaining.

* Validate inputs directly in seq_model(), parallel_model() and po_model():
  matching temperature/time vectors with at least two finite values, consecutive
  one-hour intervals, usable temperatures and parameter domains. Invalid direct
  calls now raise R errors before vector indexing. Hourly spacing uses an absolute
  tolerance of 1e-8 hours.

* Preserve the first bloom index with stopatzc = FALSE in the sequential,
  parallel, partial-overlap and original population kernels. Both population
  kernels now compute complete trajectories when forcing cuts are supplied,
  including direct C++ interface calls with stopatzc = TRUE.

* Sequential heat accumulation now uses chill after the current interval's
  update, matching the parallel model's timing. The interval that reaches the
  chill requirement contributes its full heat increment. Sequential predictions
  and calibration objectives can therefore change.

* Add a sequential model specification, named parameter schema and validation, and predict_phenology() adapter for Dynamic Model + GDH.

* Fix three-stage evaluation: check chill at the hourly budburst index, with a safe penalty when budburst is absent. Preserve pre-fix regression results and add independent calendar and station correctness tests.

* Add characteristic_to_kinetic() and kinetic_to_characteristic().
* Add kinetic/characteristic parameter-name vectors with data() support.
* Preserve convert_parameters(), convert_parameters_old_to_new(), and old/new names as compatibility aliases.
* Move the existing LarsChill inverse conversion into evalpheno; its numerical method is unchanged.
# evalpheno 0.0.1.1

* `solve_nle()` now uses stable logarithmic residuals. Added
  `solve_nle_transformed()` with E0 = exp(z1), E1 = E0 + exp(z2).
  Extreme unrepresentable energy trials are rejected without an assertion
  aborting the line search. Invalid physical inputs to `solve_nle()` now error
  explicitly, including theta_c >= 297 K.
* All evalpheno calibration conversions use the transformed solver and stable
  A0/A1 reconstruction. A small solver step alone is no longer accepted as a
  root: both residuals must be within 1e-8 and physical coefficients usable.
  Consequently numerical trajectories and some failed-fit outcomes can change.
* Added `evalpheno::convert_parameters()` with the existing conversion failure
  policies. evalpheno no longer imports LarsChill. The LarsChill source and its
  inverse conversion have not been migrated or changed in this step.

* Added `phenoflex_population()` and `sample_bud_population()`, plus the
  experimental calling convention `helper_run_pop_model()`.
* The population C++ implementation is compiled during package installation.
  Zero-variance populations, explicit populations, reproducible sampling,
  fractional forcing temperatures, and input validation are supported.
* Seeded sampling restores the caller's RNG state. Skew-normal sampling now
  honours the supplied seed and does not reset it separately for each trait.
  Consequently, skew-normal samples differ from the experimental helper.
* Negative sampled requirements produce an error; they are not truncated or
  resampled. This is an explicit validation policy, not a new distribution model.
* Forcing runs compute complete heat trajectories even when `stop_at_zc = TRUE`,
  avoiding zero-filled trajectories at cuts after bloom. The legacy zero-based
  forcing step indices and failure sentinel 9999 remain unchanged.
* Shared calendar conversion fixes leap-year offsets and distinguishes the same
  day-of-year in different years. Existing scalar and detailed wrapper return
  formats are preserved.
* General optimizer adapters reject mismatched observations and non-scalar
  predictions rather than silently recycling values. Invalid constraints are
  checked before running the model. Existing strict constraint bounds and
  missing-prediction penalties are retained.
* Corrected the E0 upper Q10 check in `evaluation_function_gensa()` and exported
  the existing `eval_phenoflex_single()` function.
* Research scripts and manuscript artifacts are excluded from source packages.
  The temperature analysis moved from `R/` to `experimental/analysis/`;
  C++ prototypes moved to `experimental/src/`.

This is a software integration of the experimental population kernel, not an
independent scientific validation of its forcing assumptions or stage criteria.



