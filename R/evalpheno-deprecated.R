#' Deprecated evalpheno interfaces
#'
#' Legacy calibration and prediction interfaces are deprecated since version
#' 0.3.1. They remain exported and retain their existing behavior, but emit a
#' standard R deprecation warning. Numerical conversion residuals, low-level
#' model functions and specialized stage and dispersion evaluators remain usable.
#'
#' @section Calibration:
#' Use \code{\link{fit_phenology}} with an explicitly supplied evaluator:
#' \code{\link{phenology_rss}} for a single model or
#' \code{\link{phenology_rss_combined}} for a model collection. Construct models
#' with \code{\link{pheno_model}}, \code{\link{combined_pheno_model}} or
#' \code{\link{pheno_model_list}}. Use named parameters and omit fixed parameters
#' from both bound vectors, or give them equal bounds. GenSA remains supported
#' through \code{optimizer = "GenSA"}; DEoptim is the default.
#'
#' The new evaluator returns a scalar loss rather than MEIGO's
#' \code{list(F, g)}. The new fitters support DEoptim and GenSA, not MEIGO.
#' Legacy Q10 constraints and custom penalties must be supplied in a custom
#' evaluation function when required. In particular, the supplied RSS evaluators
#' return \code{Inf} for missing predictions rather than substituting a date
#' such as 365. Migration can therefore change the fitted optimum.
#'
#' @section Prediction:
#' Replace weather/parameter wrappers with \code{\link{predict_phenology}} and
#' a model specification containing named parameters. The sequential, parallel
#' and partial-overlap wrappers use \code{heat_gdh("scaled")}; the custom
#' PhenoFlex wrappers use \code{heat_gdh("unscaled")}. Match this setting when
#' reusing heat requirements. Default model parameters may also differ from old
#' scripts, so supply the original values explicitly.
#'
#' Detailed prediction uses \code{basic_output = FALSE}. Its return format
#' differs from \code{custom_PhenoFlex_GDHwrapper_v2()}'s
#' \code{list(JDay, chill_heat)}. Q10 checks enabled through a legacy wrapper's
#' \code{constraints} argument are not automatically applied by the new API.
#'
#' @section Compatibility aliases:
#' Replace \code{convert_parameters()} with
#' \code{\link{characteristic_to_kinetic}} and
#' \code{convert_parameters_old_to_new()} with
#' \code{\link{kinetic_to_characteristic}}. Their return values and failure
#' policies are unchanged. Replace \code{helper_run_pop_model()} with
#' \code{\link{phenoflex_population}}, using that function's argument order or
#' named arguments. The legacy population kernel remains available.
#'
#' @name evalpheno-deprecated
#' @keywords internal
NULL
