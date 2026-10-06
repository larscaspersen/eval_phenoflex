# evalpheno 0.0.2.0

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



