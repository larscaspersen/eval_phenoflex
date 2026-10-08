#' Suggested parameter search bounds for phenology calibration
#'
#' Starting search boxes for temperate fruit phenology, not universal biological
#' limits. Only parameters in the selected model schema are returned. Collections
#' repeat ranges for indexed parameters and retain one range for shared parameters.
#' Explicit subsets of bounds passed to [fit_phenology()] still keep omitted
#' parameters fixed. Stored starting values must lie inside the bounds being used;
#' ranges are never widened automatically to include a starting value.
#'
#' Unscaled heat requirements use zc = 150--400, b1 = 20--200 and b2 = 0--800.
#' For scaled GDH these ranges are multiplied by `heat_scale`, defaulting to
#' 25 - 4 = 21. This is a fixed reference conversion: if Tb or Tu is optimized,
#' the two rectangular search boxes are not exactly equivalent at every candidate.
#' Supply another reference difference, or edit the returned heat bounds.
#'
#' Other structure ranges are yc = 10--80, s1 = 0.1--1.5, kmin = 0--1,
#' b3 = 0--0.1 and ol = 0--1. b3 has inverse chill units and is not heat-scaled.
#' Heat temperature ranges are Tb = 0--10, Tu = 10--30, Tc = 20--40;
#' chill conversion ranges are Tf = 0--10 and slope = 0.1--5.
#'
#' Characteristic chill ranges are theta_star = 279--281 K, theta_c = 286--287 K,
#' tau = 16--48 hours and pie_c = 24--50 hours. The first three follow the
#' physiological ranges reported by Egea et al. (2021); pie_c extends their
#' 24--28 hour range to accommodate broader exploratory fitting.
#' Kinetic ranges are E0 = 2500--5500,
#' E1 = 8000--16000, A0 = 1e3--1e6 and A1 = 1e13--1e19. These exploratory
#' kinetic boxes cover both the package's defaults and the chillR vignette's
#' starting values; they are not the converted image of the characteristic box.
#' Amplitudes span orders of magnitude, so characteristic fitting is generally
#' easier to interpret. Bounds alone do not guarantee plausible chill responses.
#'
#' Model constraints still apply inside these boxes, including Tb < Tu < Tc
#' and ordered stage requirements. Characteristic-to-kinetic conversion may fail
#' for some candidates; fit_phenology() excludes these with an infinite loss.
#'
#' @param model A single, population, combined or model-list specification.
#'   An ordinary list is wrapped using [pheno_model_list()]'s default sharing.
#' @param parameter_names Optional character vector selecting parameters to fit.
#'   For collections, base names such as "zc" select all corresponding indexed
#'   parameters; individual indexed names such as "zc2" select only that member.
#'   NULL selects all parameters; character() selects none.
#' @param heat_scale Positive finite reference Tu - Tb used to convert heat
#'   requirement bounds to scaled GDH units. Ignored for unscaled GDH.
#' @param temperature_unit "K" (default) or "C" for returned theta_star and
#'   theta_c bounds, including indexed names. Celsius output subtracts the native
#'   kernels' legacy offset 273; all other temperatures remain Celsius.
#' @return A list with complete or selected named numeric vectors `lower` and
#'   `upper`, in model schema order. Modify these vectors for your data and stage.
#' @seealso [parameter_schema()], [fit_phenology()]
#' @references
#' \url{https://cran.r-project.org/web/packages/chillR/vignettes/PhenoFlex.html}
#'
#' Pope and DeJong (2017), Modeling spring phenology and chilling requirements
#' using the chill overlap framework. doi:10.17660/ActaHortic.2017.1160.26.
#' Egea et al. (2021), Reducing the uncertainty on chilling requirements for
#' endodormancy breaking of temperate fruits by data-based parameter estimation
#' of the dynamic model: a test case in apricot. doi:10.1093/treephys/tpaa054.
#' Published examples motivate the parameter meanings; these default ranges
#' are package suggestions rather than limits estimated by those studies.
#' @examples
#' model <- pheno_model("parallel")
#' bounds <- default_bounds(model, c("yc", "zc", "kmin"))
#' bounds$upper["zc"] <- 700
#' stages <- stage_pheno_models(pheno_model(), c(180, 250))
#' default_bounds(stages, "zc")
#' default_bounds(model, c("theta_star", "theta_c"), temperature_unit = "C")
#' @md
#' @export
default_bounds <- function(model, parameter_names = NULL, heat_scale = 21,
                            temperature_unit = c("K", "C")) {
  temperature_unit <- match.arg(temperature_unit)
  if (is.list(model) && !inherits(model, "pheno_model_list") && length(model) &&
      all(vapply(model, inherits, logical(1), what = "pheno_model")))
    model <- pheno_model_list(model)
  schema <- parameter_schema(model)
  base_names <- if ("base_name" %in% names(schema)) schema$base_name else schema$name
  base_model <- if (inherits(model, "pheno_model_list")) model[[1]] else
    if (inherits(model, "combined_pheno_model") || inherits(model, "population_pheno_model"))
      model$model else model
  if (!is.numeric(heat_scale) || length(heat_scale) != 1L ||
      !is.null(dim(heat_scale)) || !is.finite(heat_scale) || heat_scale <= 0)
    stop("heat_scale must be one positive finite reference Tu - Tb.", call. = FALSE)
  ranges <- list(
    lower = c(yc = 10, zc = 150, s1 = 0.1, kmin = 0, b1 = 20, b2 = 0,
              b3 = 0, ol = 0, theta_star = 279, theta_c = 286, tau = 16,
              pie_c = 24, E0 = 2500, E1 = 8000, A0 = 1e3, A1 = 1e13,
              Tf = 0, slope = 0.1, Tb = 0, Tu = 10, Tc = 20),
    upper = c(yc = 80, zc = 400, s1 = 1.5, kmin = 1, b1 = 200, b2 = 800,
              b3 = 0.1, ol = 1, theta_star = 281, theta_c = 287, tau = 48,
              pie_c = 50, E0 = 5500, E1 = 16000, A0 = 1e6, A1 = 1e19,
              Tf = 10, slope = 5, Tb = 10, Tu = 30, Tc = 40)
  )
  selected <- rep(TRUE, nrow(schema))
  if (!is.null(parameter_names)) {
    if (!is.character(parameter_names) || !is.null(dim(parameter_names)) ||
        anyNA(parameter_names) || anyDuplicated(parameter_names) ||
        !all(parameter_names %in% c(schema$name, base_names)))
      stop("parameter_names must contain distinct schema or base parameter names.", call. = FALSE)
    selected <- schema$name %in% parameter_names | base_names %in% parameter_names
  }
  if (!any(selected)) return(list(lower = numeric(), upper = numeric()))
  heat <- base_names %in% c("zc", "b1", "b2")
  lapply(ranges, function(values) {
    bounds <- stats::setNames(unname(values[base_names]), schema$name)
    if (identical(base_model$heat$scaling, "scaled")) bounds[heat] <- bounds[heat] * heat_scale
    if (temperature_unit == "C") {
      temperatures <- base_names %in% c("theta_star", "theta_c")
      bounds[temperatures] <- bounds[temperatures] - 273
    }
    if (any(!is.finite(bounds))) stop("Default bounds overflowed; use a smaller heat_scale.", call. = FALSE)
    bounds[selected]
  })
}
