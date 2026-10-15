#' Specify a phenology model
#'
#' A model specification is an S3 list storing algorithms and named parameter
#' values. It supports sequential, linear parallel,
#' partial-overlap and PhenoFlex coupling with the Dynamic Model, chill hours,
#' Utah or positive Utah chill and GDH heat. With chill_none(), only heat is
#' accumulated and the chosen coupling is bypassed.
#' Printing shows the structure, submodel names and parameters in a compact form.
#' @encoding UTF-8
#' @details
#' \subsection{Model structures}{
#' The structure controls how accumulated chill affects heat accumulation or the
#' heat requirement. The chill and heat submodels can be selected independently
#' of the structure, so a combination need not reproduce the original publication.
#' \describe{
#'   \item{\code{sequential}}{Heat accumulation starts once the chill requirement
#'   \code{yc} is met, following the sequential approach of Ashcroft et al. (1977).}
#'   \item{\code{parallel}}{Heat effectiveness increases linearly with accumulated
#'   chill, from \code{kmin} to full effectiveness at \code{yc}, following
#'   Hänninen and Kramer (2007), equation B4b.}
#'   \item{\code{partial_overlap}}{After \code{yc} is met, additional chill reduces
#'   the heat requirement until the overlap threshold is reached, following
#'   Pope et al. (2014).}
#'   \item{\code{phenoflex}}{Heat effectiveness increases sigmoidally with
#'   accumulated chill, following PhenoFlex (Luedeling et al., 2021).}
#' }
#' }
#' \subsection{Chill models}{
#' \describe{
#'   \item{\code{chill_dynamic()}}{The Dynamic Model represents chill accumulation
#'   as a two-step process and expresses chill in Chill Portions
#'   (Fishman et al., 1987a,b; Erez et al., 1990).}
#'   \item{\code{chill_hours()}}{The Chilling Hours Model counts hours within a
#'   temperature range, by default 0--7.2 degrees Celsius (Weinberger, 1950).}
#'   \item{\code{chill_utah()}}{The Utah Model assigns temperature-dependent
#'   weights, including negative contributions at high temperatures
#'   (Richardson et al., 1974).}
#'   \item{\code{chill_positive_utah()}}{The Positive Utah modification adds only
#'   positive daily chill totals (Linsley-Noakes et al., 1994).}
#'   \item{\code{chill_none()}}{A heat-only option that bypasses chill-dependent
#'   coupling.}
#' }
#' See \code{\link{chill_hours}} for implementation details of the alternative
#' chill models, and \code{\link{chill_custom}} for user-defined chill submodels.
#' }
#' \subsection{Dynamic Model parameters}{
#' The Dynamic Model has six parameters. Two, \code{Tf} and \code{slope}, control
#' the conversion of a reversible intermediate into irreversible Chill Portions.
#' The other four, \code{E0}, \code{E1}, \code{A0} and \code{A1}, govern the
#' formation and destruction of the intermediate. These are referred to here as
#' \emph{kinetic parameters} because they determine the temperature-dependent
#' rates of these processes.
#'
#' \strong{Why use characteristic parameters?}
#' Choosing values other than the defaults of Erez et al. (1990) can be difficult.
#' The kinetic parameters do not directly describe observable quantities such as
#' the optimal chilling temperature, and their values differ greatly in magnitude.
#' For example, the default value of \code{A1} is \code{2.567e18}.
#' This makes them difficult to interpret and can complicate calibration alongside
#' structural parameters such as chill and heat requirements.
#'
#' The model therefore also accepts four \emph{characteristic parameters}, which
#' describe responses to constant and alternating temperatures in the experiments
#' underlying the Dynamic Model (Fishman et al., 1987a,b):
#' \describe{
#'   \item{\code{theta_star}}{Optimal temperature for chill accumulation under
#'   constant conditions.}
#'   \item{\code{theta_c}}{Upper critical temperature, above which no chill
#'   accumulates under constant conditions.}
#'   \item{\code{tau}}{Time (hours) required to accumulate one Chill Portion at
#'   the optimal temperature.}
#'   \item{\code{pie_c}}{Critical period (hours) of the reference alternating
#'   temperature cycle at which chill accumulation ceases.}
#' }
#'
#' \strong{Supplying and storing parameters.}
#' Parameter names select the representation automatically. Supply either kinetic
#' or characteristic parameters; mixing the two sets raises an error. If neither
#' set is supplied, the original kinetic defaults are used. When characteristic
#' parameters are supplied, omitted values use the characteristic defaults.
#' These form a separate starting set, rather than an equivalent representation
#' of the original kinetic defaults.
#'
#' Characteristic inputs are converted to kinetic coefficients when the model is
#' constructed. The characteristic values are retained in
#' \code{model$chill$characteristic}, while \code{model$parameters} always stores
#' the kinetic coefficients.
#'
#' \strong{Temperature units.}
#' For \code{theta_star} and \code{theta_c}, values from 0 to 20 (inclusive) are
#' interpreted as degrees Celsius and converted to Kelvin using the native
#' implementation's legacy offset: \code{K = Celsius + 273}. Values outside this
#' range are interpreted as Kelvin. Both characteristic temperatures are stored
#' in Kelvin; all other temperature parameters use degrees Celsius.
#' }
#' \subsection{Choosing a chill model}{
#' Based on extensive model comparisons and a review of the relevant literature,
#' I do not recommend using the Chilling Hours or Utah Models, especially in warm
#' climates. The Dynamic Model (Chill Portions), although far from perfect, seems
#' more reliable.
#'
#' This recommendation is consistent with the comparisons by Luedeling and Brown
#' (2011) and Fernandez et al. (2020), which document strong differences among
#' chill models and limitations of Chilling Hours and Utah estimates in warm
#' climates. Model performance should still be evaluated for the species,
#' cultivar and climate of interest.
#' }
#' @references
#' Ashcroft, G. L., Richardson, E. A. and Seeley, S. D. (1977).
#' A statistical method of determining chill unit and growing degree hour
#' requirements for deciduous fruit trees. HortScience 12(4), 347--348.
#' \doi{10.21273/HORTSCI.12.4.347}.
#'
#' Hänninen, H. and Kramer, K. (2007). A framework for modelling the annual
#' cycle of trees in boreal and temperate regions. Silva Fennica 41(1), 167--205.
#' \doi{10.14214/sf.313}.
#'
#' Pope, K. S., Da Silva, D., Brown, P. H. and DeJong, T. M. (2014).
#' A biologically based approach to modeling spring phenology in temperate
#' deciduous trees. Agricultural and Forest Meteorology 198--199, 15--23.
#' \doi{10.1016/j.agrformet.2014.07.009}.
#'
#' Luedeling, E., Schiffers, K., Fohrmann, T. and Urbach, C. (2021).
#' PhenoFlex -- an integrated model to predict spring phenology in temperate
#' fruit trees. Agricultural and Forest Meteorology 307, 108491.
#' \doi{10.1016/j.agrformet.2021.108491}.
#'
#' Fishman, S., Erez, A. and Couvillon, G. A. (1987a). The temperature dependence
#' of dormancy breaking in plants: mathematical analysis of a two-step model
#' involving a cooperative transition. Journal of Theoretical Biology 124(4),
#' 473--483. \doi{10.1016/S0022-5193(87)80221-7}.
#' 
#' Fishman, S., Erez, A. and Couvillon, G. A. (1987b). The temperature dependence
#' of dormancy breaking in plants: computer simulation studied under controlled temperatures
#' Journal of Theoretical Biology 126(3),
#' 473--483. \doi{10.1016/S0022-5193(87)80237-0}.
#'
#' Erez, A., Fishman, S., Linsley-Noakes, G. C. and Allan, P. (1990).
#' The Dynamic Model for rest completion in peach buds. Acta Horticulturae 276,
#' 165--174. \doi{10.17660/ActaHortic.1990.276.18}.
#'
#' Weinberger, J. H. (1950). Chilling requirements of peach varieties.
#' Proceedings of the American Society for Horticultural Science 56, 122--128.
#'
#' Richardson, E. A., Seeley, S. D. and Walker, D. R. (1974). A model for
#' estimating the completion of rest for 'Redhaven' and 'Elberta' peach trees.
#' HortScience 9(4), 331--332. \doi{10.21273/HORTSCI.9.4.331}.
#'
#' Linsley-Noakes, G. C., Allan, P. and Matthee, G. (1994). Modification of rest
#' completion prediction models for improved accuracy in South African stone
#' fruit orchards. Journal of the Southern African Society for Horticultural
#' Sciences 4, 13--15.
#'
#' Luedeling, E. and Brown, P. H. (2011). A global analysis of the comparability
#' of winter chill models for fruit and nut trees. International Journal of
#' Biometeorology 55, 411--421. \doi{10.1007/s00484-010-0352-y}.
#'
#' Fernandez, E., Whitney, C. and Luedeling, E. (2020). The importance of chill
#' model selection -- a multi-site analysis. European Journal of Agronomy 119,
#' 126103. \doi{10.1016/j.eja.2020.126103}.
#' @param structure One of "sequential", "parallel", "partial_overlap" or "phenoflex".
#' @param chill Chill specification returned by chill_dynamic(), chill_hours(),
#'   chill_utah(), chill_positive_utah(), chill_none() or chill_custom().
#' @param heat Heat specification returned by heat_gdh().
#' @param parameters Named numeric vector specifying any subset of the model's
#' parameter names. Omitted parameters use default_parameters() for the selected
#' structure, chill and heat specifications. NULL uses all defaults. The complete
#' vector is validated after filling defaults.
#' See \emph{Dynamic Model parameters} in Details for the kinetic and
#' characteristic representations, their defaults, and temperature units.
#' @return pheno_model() returns a pheno_model list. Component constructors return
#' lists; validate_model_spec() invisibly returns TRUE or raises an error.
#' @examples
#' #construct a model (by default phenoflex)
#' model <- pheno_model()
#' 
#' #' #inspect model structure
#' print(model)
#' 
#' #visualize model structure and submodels
#' plot(model)
#' 
#' #inspect model parameters and their meaning
#' parameter_schema(model)
#' 
#' #' #adjusting model parameters in a named vector
#' model2 <- pheno_model(parameters = c('yc' = 70))
#' 
#' #choosing different model structure and chill submodel
#' parallel <- pheno_model(structure = "parallel", chill = chill_hours() )
#' 
#' #prediction
#' data("cka_weather", package = "evalpheno")
#' hourly <- chillR::stack_hourly_temps(
#' cka_weather, latitude = attr(cka_weather, "latitude")
#' )$hourtemps
#' years <- 2008:2009
#' seasons <- chillR::genSeasonList(hourly, mrange = c(8, 6), years = years)
#' names(seasons) <- as.character(years)
#' season <- seasons[[1]]
#' predict_phenology(model = model, weather = seasons)
#' 
#' @export
pheno_model <- function(structure = "phenoflex", chill = chill_dynamic(), heat = heat_gdh(),
                        parameters = NULL) {
  model <- list(structure = structure, chill = chill, heat = heat)
  validate_model_spec(model)
  class(model) <- "pheno_model"
  if (identical(chill$name, "dynamic"))
    model$chill$parameterization <- .dynamic_parameterization(names(parameters))
  defaults <- default_parameters(model)
  if (is.null(parameters)) {
    parameters <- defaults
  } else {
    if (!is.numeric(parameters) || !is.null(dim(parameters)) ||
        (length(parameters) > 0L && is.null(names(parameters))) ||
        anyNA(names(parameters)) || anyDuplicated(names(parameters)) ||
        !all(names(parameters) %in% names(defaults)))
      stop("parameters must be a named numeric vector using schema names, without duplicates.",
           call. = FALSE)
    if (length(parameters) < length(defaults)) {
      defaults[names(parameters)] <- parameters
      parameters <- defaults
    }
  }
  parameters <- .normalize_characteristic_temperatures(parameters)
  validate_parameters(model, parameters)
  model$parameters <- parameters
  .with_dynamic_parameterization(model, "kinetic")
}

#' @rdname pheno_model
#' @param x A model specification returned by pheno_model().
#' @param ... Additional arguments passed to print() for the parameter vector.
#' @export
print.pheno_model <- function(x, ...) {
  .print_pheno_components(x)
  cat("Parameters:\n")
  print(x$parameters, ...)
  invisible(x)
}

.print_pheno_components <- function(x) {
  structure_name <- if (identical(x$chill$name, "none")) "Heat only" else switch(x$structure,
    phenoflex = "PhenoFlex", sequential = "Sequential",
    parallel = "Parallel", partial_overlap = "Partial overlap")
  cat("Structure:      ", structure_name, "\n", sep = "")
  chill_name <- switch(x$chill$name, dynamic = "Dynamic Model", chill_hours = "Chill hours",
    utah = "Utah", positive_utah = "Positive Utah", none = "None", custom = x$chill$label)
  cat("Chill submodel: ", chill_name, "\n", sep = "")
  cat("Heat submodel:  GDH\n")
}

#' @rdname pheno_model
#' @param parameterization Legacy optional representation argument for
#' chill_dynamic(). Construction now infers the representation from parameter
#' names, so this argument does not select stored parameters. Select the fitting
#' representation in fit_phenology() instead.
#' @export
chill_dynamic <- function(parameterization = NULL) {
  if (!is.null(parameterization)) match.arg(parameterization, c("characteristic", "kinetic"))
  list(name = "dynamic", parameterization = "kinetic")
}

#' @rdname pheno_model
#' @param scaling "unscaled" returns Growing Degree Hour (GDH) units without scaling via (Tu-Tb).
#' "scaled" returns units multiplied by (Tu-Tb). The standard implementation of PhenoFlex
#' uses unscaled GDH, while most other studies use scaled GDH.
#' @export
heat_gdh <- function(scaling = c("unscaled", "scaled")) list(name = "gdh",
                            scaling = match.arg(scaling))

#' @rdname pheno_model
#' @param model Model specification.
#' @export
validate_model_spec <- function(model) {
  if (inherits(model, "stage_pheno_model"))
    return(validate_model_spec(.stage_model_list(model)))
  if (inherits(model, "pheno_model_list"))
    return(.validate_pheno_model_list(model))
  if (inherits(model, "combined_pheno_model"))
    return(.validate_combined_model_spec(model))
  if (inherits(model, "population_pheno_model"))
    return(validate_model_spec(model$model))
  if (!is.list(model) ||
      !(identical(sort(names(model)), sort(c("structure", "chill", "heat"))) ||
        identical(sort(names(model)), sort(c("structure", "chill", "heat", "parameters")))))
    stop("model must contain structure, chill and heat, with optional parameters.", call. = FALSE)
  if (!is.character(model$structure) ||
      length(model$structure) != 1L ||
      is.na(model$structure) ||
      !model$structure %in% c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    stop("Supported structures are sequential, parallel, partial_overlap and phenoflex.")
  }
  chill <- model$chill
  if (!is.list(chill) || !is.character(chill$name) || length(chill$name) != 1L ||
      is.na(chill$name) || !chill$name %in% c("dynamic", "chill_hours", "utah", "positive_utah", "none", "custom"))
    stop("Use chill_dynamic(), chill_hours(), chill_utah(), chill_positive_utah(), chill_none() or chill_custom().",
         call. = FALSE)
  if (chill$name == "dynamic") {
    if (anyDuplicated(names(chill)) || !all(c("name", "parameterization") %in% names(chill)) ||
        !all(names(chill) %in% c("name", "parameterization", "characteristic")) ||
        !(identical(chill$parameterization, "characteristic") ||
          identical(chill$parameterization, "kinetic")))
      stop("Use chill_dynamic() with characteristic or kinetic parameterization.", call. = FALSE)
    if (!is.null(chill$characteristic)) {
      p <- chill$characteristic
      if (!is.numeric(p) || !is.null(dim(p)) || length(p) != 4L ||
          !identical(names(p), .dynamic_characteristic_names) || any(!is.finite(p)))
        stop("chill$characteristic must contain the four finite characteristic parameters.", call. = FALSE)
      .validate_chill_params(unname(p))
    }
  } else if (chill$name == "custom") {
    .validate_custom_chill_spec(chill)
  } else if (!identical(names(chill), "name")) {
    stop("Alternative chill specifications must contain only name; set parameters in pheno_model().",
         call. = FALSE)
  }

  # Check the heat specification's structure
  if (!is.list(model$heat) ||
      !identical(sort(names(model$heat)), c("name", "scaling"))) {
    stop(
      "heat must contain exactly name and scaling; use heat_gdh().",
      call. = FALSE
    )
  }
  
  # Check the supported heat model
  if (!identical(model$heat$name, "gdh")) {
    stop("Only the GDH heat model is currently supported.", call. = FALSE)
  }
  
  # Check scaling
  if (!identical(model$heat$scaling, "scaled") &&
      !identical(model$heat$scaling, "unscaled")) {
    stop(
      "GDH scaling must be 'scaled' or 'unscaled'.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


#' Named parameters for a model specification
#'
#' Parameters use plain names; the schema records their component separately.
#' Defaults are illustrative starting
#' values, not calibrated estimates or optimization bounds. The kinetic defaults
#' use the original Dynamic Model coefficients; characteristic defaults are a
#' separate starting set, not an equivalent representation of those defaults.
#' Validation checks model domains, not experiment-specific calibration bounds.
#' @param model Model specification from pheno_model(), [population_pheno_model()]
#' or [combined_pheno_model()]. Combined specifications expand selected structure
#' parameters by cultivar while leaving other parameters shared.
#' A pheno_model_list() stores ordinary simple models and their parameter sharing.
#' A stage_pheno_models() specification stores shared parameters and stage heat
#' requirements once, expanding zc to indexed calibration names.
#' @param parameters Named numeric vector with exactly the names in the schema;
#' order is arbitrary. Matrices and arrays are not parameter vectors.
#' Complete Dynamic Model vectors may use either representation; names select
#' it automatically and mixing characteristic and kinetic names is an error.
#' Named theta_star and theta_c values in 0--20 are interpreted as Celsius
#' using the legacy +273 offset; Kelvin values remain supported. Validation
#' checks the converted values without modifying the supplied vector.
#' @return parameter_schema() returns a data frame with name, component, default,
#' description and unit columns. Combined models also include base_name and cultivar
#' (NA for shared parameters); model lists use member instead of cultivar.
#' default_parameters() returns a named numeric vector.
#' validate_parameters() invisibly returns TRUE or raises an error.
#' @details
#' Descriptions reflect the selected model structure. Chill requirements use
#' Chill Portions for the Dynamic Model, Chilling Hours for chill hours, and
#' Utah Chill Units for Utah and positive Utah. Heat requirements (zc, b1 and b2)
#' use degree C h for scaled GDH and h for unscaled GDH, which accumulates
#' temperature-weighted hours. The units of s1 and b3 are the inverse chill unit.
#' Dimensionless parameters are labelled explicitly.
#'
#' Characteristic temperatures theta_star and theta_c are reported in K, matching
#' their stored values; the constructor also accepts Celsius values as described
#' in pheno_model(). E0 and E1 are activation energies divided by the gas constant,
#' expressed in K. Other temperature thresholds use degrees Celsius.
#'
#' For user-defined chill submodels, yc uses custom chill units. Custom parameters
#' receive a generic description and NA units because their meaning and units
#' cannot be inferred from the callback.
#' @param parameterization Optional Dynamic Model representation for a calibration
#' schema or defaults. NULL uses the model's representation (kinetic for constructed
#' models); "characteristic" exposes characteristic calibration names.
#' @export
parameter_schema <- function(model, parameterization = NULL) {
  if (!is.null(parameterization)) model <- .with_dynamic_parameterization(model, parameterization, use_defaults = TRUE)
  if (inherits(model, "stage_pheno_model")) model <- .stage_model_list(model)
  if (inherits(model, "pheno_model_list"))
    return(.pheno_model_list_schema(model))
  if (inherits(model, "combined_pheno_model"))
    return(.combined_parameter_schema(model))
  if (inherits(model, "population_pheno_model"))
    return(parameter_schema(model$model))
  #check model
  validate_model_spec(model)
  

  # Starting requirements are illustrative in the chosen chill model's units.
  yc <- switch(model$chill$name, chill_hours = 1000, utah = 1200, positive_utah = 1200, 40)
  if(model$chill$name == 'dynamic'){
    
    chill <- if (model$chill$parameterization == "characteristic")
      c(.dynamic_characteristic_defaults, Tf = 4, slope = 1.6)
    else c(.dynamic_kinetic_defaults, Tf = 4, slope = 1.6)
  
  } else {
    chill <- switch(model$chill$name,
      chill_hours = c(CH_lower = 0, CH_upper = 7.2),
      utah = stats::setNames(c(0, 0.5, 1, 0.5, 0, -0.5, -1), paste0("utah_w", 1:7)),
      positive_utah = stats::setNames(c(0, 0.5, 1, 0.5, 0, -0.5, -1), paste0("pos_utah_w", 1:7)),
      none = c(heat_start_index = 1), custom = model$chill$defaults)
  }
  
  #check heat
  if(model$heat$name == 'gdh'){
    heat <- c(Tb = 4, Tu = 25, Tc = 36)
    heat_requirements <- c(zc = 6000, b1 = 1119, b2 = 8677)
    if(model$heat$scaling == 'unscaled'){
      heat_requirements <- heat_requirements /
        (heat[["Tu"]] - heat[["Tb"]])
    }
  } else {
    stop("Unknown heat submodel.")
  }
  
  # 1. Parameters belonging to the model structure
  structure_parameters <- if (model$chill$name == "none") c(zc = heat_requirements[["zc"]]) else switch(
    model$structure,
    sequential = c(yc = yc, zc = heat_requirements[["zc"]]),
    parallel   = c(yc = yc, zc = heat_requirements[["zc"]], kmin = 0.1),
    partial_overlap = c(yc = yc, b1 = heat_requirements[["b1"]], b2 = heat_requirements[["b2"]], b3 = 0.01119,
                        ol = 0.75),
    phenoflex = c(yc = yc, zc = heat_requirements[["zc"]], s1 = 0.5),
    stop("Unknown model structure.")
  )
  
  values <- c(structure_parameters, chill, heat)
  
  
  # Create one component label for every parameter.
  component <- c(
    rep("structure", length(structure_parameters)),
    rep("chill", length(chill)),
    rep("heat", length(heat))
  )
  
  metadata <- .parameter_metadata(model, names(values))
  
  data.frame(
    name = names(values),
    component = component,
    default = unname(values),
    description = unname(metadata$description),
    unit = unname(metadata$unit),
    stringsAsFactors = FALSE
  )
}

.parameter_metadata <- function(model, parameter_names) {
  chill_unit <- switch(model$chill$name,
    dynamic = "Chill Portions", chill_hours = "Chilling Hours",
    utah = "Utah Chill Units", positive_utah = "Utah Chill Units",
    custom = "custom chill units", none = NA_character_)
  heat_unit <- if (model$heat$scaling == "scaled") "degree C h" else "h"
  inverse_chill_unit <- paste0("1 / (", chill_unit, ")")

  description <- c(
    yc = switch(model$structure,
      phenoflex = "Chill requirement",
      parallel = "Chill requirement",
      "Chill requirement"),
    zc = "Heat requirement",
    s1 = "Slope parameter for transition from chill accumulation to heat accumulation",
    kmin = "Minimum relative heat effectiveness at zero chill",
    b1 = "Baseline heat requirement",
    b2 = "Additional heat requirement at the start of overlap",
    b3 = "Rate of heat-requirement decline with additional chill",
    ol = "Overlap multiplier for the baseline heat requirement",
    E0 = "time-independent activation energy of forming the precursor to the dormancy-breaking factor",
    E1 = "time-independent activation energy of destroying the precursor to the dormancy-breaking factor",
    A0 = "amplitude of the (hypothetical) process involved in forming the precursor to the dormancy-breaking factor",
    A1 = "amplitude of the (hypothetical) process involved in destroying the precursor to the dormancy-breaking factor",
    Tf = "Transition temperature for irreversible chill conversion",
    slope = "Slope of temperature-dependent chill conversion",
    theta_star = "Optimal temperature for chill accumulation",
    theta_c = "Upper critical temperature for chill accumulation",
    tau = "Time to accumulate one Chill Portion at the optimal temperature",
    pie_c = "Critical temperature-cycle period for chill negation",
    CH_lower = "Lower temperature threshold for Chilling Hours",
    CH_upper = "Upper temperature threshold for Chilling Hours",
    heat_start_index = "First hourly interval eligible for heat accumulation",
    Tb = "Base temperature for heat accumulation",
    Tu = "Optimal temperature for heat accumulation",
    Tc = "Upper critical temperature for heat accumulation")
  unit <- c(
    yc = chill_unit, zc = heat_unit, s1 = inverse_chill_unit,
    kmin = "dimensionless", b1 = heat_unit, b2 = heat_unit,
    b3 = inverse_chill_unit, ol = "dimensionless",
    E0 = "dimensionless", E1 = "dimensionless", A0 = "dimensionless", A1 = "dimensionless",
    Tf = "degree C", slope = "K^-1", theta_star = "K", theta_c = "K",
    tau = "h", pie_c = "h", CH_lower = "degree C", CH_upper = "degree C",
    heat_start_index = "row index", Tb = "degree C", Tu = "degree C", Tc = "degree C")

  if (model$chill$name %in% c("utah", "positive_utah")) {
    weights <- paste0(if (model$chill$name == "utah") "utah_w" else "pos_utah_w", 1:7)
    zones <- c("T <= 1.4", "1.4 < T <= 2.4", "2.4 < T <= 9.1",
      "9.1 < T <= 12.4", "12.4 < T <= 15.9", "15.9 < T <= 18", "18 < T <= 1000")
    description[weights] <- paste0("Chill weight for ", zones, " degrees C")
    unit[weights] <- "Utah Chill Units / h"
  } else if (model$chill$name == "custom") {
    custom_names <- names(model$chill$defaults)
    description[custom_names] <- paste("Custom chill parameter:", custom_names)
    unit[custom_names] <- NA_character_
  }

  list(description = description[parameter_names], unit = unit[parameter_names])
}

#' @rdname parameter_schema
#' @export
default_parameters <- function(model, parameterization = NULL) {
  schema <- parameter_schema(model, parameterization)
  stats::setNames(schema$default, schema$name)
}

#' @rdname parameter_schema
#' @export
validate_parameters <- function(model, parameters) {
  if (inherits(model, "stage_pheno_model")) model <- .stage_model_list(model)
  if (inherits(model, "pheno_model_list"))
    return(.validate_model_list_parameters(model, parameters))
  if (inherits(model, "combined_pheno_model"))
    return(.validate_combined_parameters(model, parameters))
  if (inherits(model, "population_pheno_model"))
    return(validate_parameters(model$model, parameters))
  if (identical(model$chill$name, "dynamic")) {
    representation <- .dynamic_parameterization(names(parameters), fallback = model$chill$parameterization)
    if (!identical(representation, model$chill$parameterization)) {
      # Validation only needs a schema; avoid converting the stored coefficients.
      model$chill$parameterization <- representation
    }
  }
  schema <- parameter_schema(model)
  if (!is.numeric(parameters) || !is.null(dim(parameters)) || is.null(names(parameters)) ||
      anyNA(names(parameters)) || anyDuplicated(names(parameters)) ||
      length(parameters) != nrow(schema) || !setequal(names(parameters), schema$name))
    stop("parameters must be a named numeric vector with exactly the schema names, without duplicates.", call. = FALSE)
  if (any(!is.finite(parameters))) stop("All parameter values must be finite.", call. = FALSE)
  p <- .normalize_characteristic_temperatures(parameters)
  if (model$chill$name == "none") {
    if (p[["zc"]] <= 0) stop("zc must be positive.", call. = FALSE)
    if (p[["heat_start_index"]] < 1)
      stop("heat_start_index must be at least one.", call. = FALSE)
  } else switch(
    model$structure,
    
    sequential = {
      if (p[["yc"]] <= 0 || p[["zc"]] <= 0)
        stop("yc and zc must be positive.", call. = FALSE)
    },
    
    parallel = {
      if (p[["yc"]] <= 0 || p[["zc"]] <= 0)
        stop("yc and zc must be positive.", call. = FALSE)
      
      if (p[["kmin"]] < 0 || p[["kmin"]] > 1)
        stop("kmin must be between 0 and 1.", call. = FALSE)
    },
    
    partial_overlap = {
      if (p[["yc"]] <= 0 || p[["b1"]] <= 0)
        stop("yc and b1 must be positive.", call. = FALSE)
      
      if (p[["b2"]] < 0 || p[["b3"]] < 0)
        stop("b2 and b3 must be non-negative.", call. = FALSE)
      
      if (p[["ol"]] < 0)
        stop("ol must be non-negative.", call. = FALSE)
      
      if (!is.finite(p[["b1"]] * p[["ol"]]))
        stop("b1 * ol must be finite.", call. = FALSE)
    },
    
    phenoflex = {
      if (p[["yc"]] <= 0 || p[["zc"]] <= 0 ||
          p[["s1"]] <= 0)
        stop("yc, zc and s1 must be positive.", call. = FALSE)
    },
    
    parallel_landsberg = {
      # Assumes your schema uses the proposed name chill_scale
      if (p[["chill_scale"]] <= 0 || p[["zc"]] <= 0)
        stop("chill_scale and zc must be positive.", call. = FALSE)
    },
    
    stop("Unknown model structure.", call. = FALSE)
  )
  
  if(model$chill$name == 'dynamic'){
    if (p[["slope"]] <= 0 || p[["Tf"]] <= -273)
      stop("slope must be positive and Tf must exceed -273 degrees C.", call. = FALSE)
    
    if (model$chill$parameterization == "characteristic") {
      .validate_chill_params(unname(p[c("theta_star", "theta_c", "tau", "pie_c")]))
    } else if (!.valid_chill_energies(unname(p[c("E0", "E1")])) ||
               p[["A0"]] <= 0 || p[["A1"]] <= 0) {
      stop("Kinetic parameters require 0 < E0 < E1, A0 > 0 and A1 > 0.", call. = FALSE)
    }
  } else if (model$chill$name == "chill_hours") {
    if (p[["CH_lower"]] <= -273 || p[["CH_lower"]] >= p[["CH_upper"]])
      stop("Chill-hours thresholds must satisfy -273 < CH_lower < CH_upper in Celsius.", call. = FALSE)
  } else if (model$chill$name == "custom") {
    .validate_custom_chill_parameters(model$chill, p[names(model$chill$defaults)])
  }
  
  if(model$heat$name == 'gdh'){
    if (!(p[["Tb"]] < p[["Tu"]] && p[["Tu"]] < p[["Tc"]]))
      stop("Heat parameters must satisfy Tb < Tu < Tc.", call. = FALSE)
    
    if (p[["Tb"]] <= -273)
      stop("Tb must exceed -273 degrees C.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Predict phenology for one or more models and seasons
#'
#' Calculates shared chill and potential heat and applies the selected C++
#' structure. Characteristic parameters use the shared conversion implementation.
#' Conversion failure raises an error; no bloom is represented by index zero.
#' Population specifications dispatch to [predict_population_phenology()].
#' Combined specifications use their cultivar-specific single models, with the
#' same weather formats as model lists.
#' Model lists accept one common seasonlist for all models, or nested seasonlists
#' matched by position to the models, and return results grouped by model.
#' Ordinary lists may contain models with different specifications and values;
#' prediction does not impose parameter sharing between them.
#' Ensemble specifications from [pheno_ensemble()] aggregate independent fits.
#' With basic_output = FALSE they return date/spread summaries and individual
#' member dates, rather than chill/heat trajectories.
#' @param model Model specification from pheno_model(), a non-empty ordinary list
#' of these models, a stage_pheno_models() multistage specification,
#' or a pheno_model_list() calibration collection.
#' Also accepts population_pheno_model() and combined_pheno_model() specifications.
#' A pheno_ensemble() uses its stored member parameters and aggregation policy;
#' explicit parameter overrides are unsupported for ensembles.
#' Fitted results from fit_phenology() and fit_phenology_cv() are accepted directly.
#' A CV result predicts with its stored fold-model ensemble, even when a separate
#' full-data refit is present in calibration$refit.
#' @param weather Data frame with Temp (degrees C), Year and JDay. Supply complete,
#' consecutive days, each with 24 hourly rows in chronological order. If Hour is
#' present it must be 0:23 for each day. Without Hour, within-day order is assumed.
#' A non-empty list of these data frames supplies multiple seasons. With a list
#' of models, a single data frame or flat seasonlist is evaluated by every model;
#' a nested list supplies one non-empty seasonlist per model, matched by position.
#' Combined models also accept these formats, with nested lists matched to cultivars.
#' @param parameters Named numeric vector; defaults to the model's stored parameters.
#' An explicit vector overrides stored values for this prediction only.
#' Complete Dynamic Model vectors may use either representation; input names
#' select it automatically and mixing the two sets is an error.
#' Named theta_star and theta_c values in 0--20 are interpreted as Celsius
#' using the same +273 convention as pheno_model(); other values are Kelvin.
#' For an ordinary model list, supply a list of complete parameter vectors in
#' model order, or one complete vector to use for every model. For calibration
#' collections, supply the complete indexed calibration vector.
#' @param stopatzc Stop at the first bloom if TRUE. FALSE completes trajectories
#' while retaining the first bloom index.
#' @param basic_output If TRUE, return a numeric fractional bloom Julian day
#' (NA for no bloom). FALSE returns the index and trajectories.
#' @return For a single model, a numeric fractional Julian day, or NA for no bloom.
#' With basic_output = FALSE, a list with bloomindex (one-based row; zero for no
#' bloom), chill and z. Conversion uses [return_JDay()]. Population specifications
#' return the result described in [predict_population_phenology()].
#' A single model with a seasonlist returns a numeric vector (or a list of
#' detailed results). A model list returns a list in model order, containing a
#' scalar per model for one data frame, or a prediction vector per model for
#' seasonlists. With basic_output = FALSE these contain detailed results instead.
#' Combined specifications use the same result shapes as model lists.
#' Model names label common-weather results; nested weather names label the outer
#' results when supplied. Season names are retained. Population seasonlists
#' return a list of population results.
#' @export
predict_phenology <- function(model, weather, 
                              parameters = model_parameters(model),
                              stopatzc = TRUE, 
                              basic_output = TRUE) {
  for (flag in list(stopatzc, basic_output)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("Output flags must be single non-missing logical values.", call. = FALSE)
  }
  if (inherits(model, "phenology_fit") || inherits(model, "phenology_cv")) {
    if (missing(parameters))
      return(predict_phenology(model$model, weather,
                              stopatzc = stopatzc, basic_output = basic_output))
    return(predict_phenology(model$model, weather, parameters,
                            stopatzc = stopatzc, basic_output = basic_output))
  }
  if (inherits(model, "pheno_ensemble")) {
    if (!missing(parameters))
      stop("Ensembles use their stored member parameters; parameter overrides are unsupported.", call. = FALSE)
    return(.predict_pheno_ensemble(model, weather, stopatzc, basic_output))
  }
  if (inherits(model, "pheno_model_list"))
    return(.predict_pheno_model_list(model, weather, parameters, stopatzc, basic_output))
  if (inherits(model, "stage_pheno_model"))
    return(.predict_pheno_model_list(model, weather, parameters, stopatzc, basic_output))
  if (inherits(model, "combined_pheno_model"))
    return(.predict_pheno_model_list(model, weather, parameters, stopatzc, basic_output))
  if (is.list(model) && !is.data.frame(model) &&
      !inherits(model, "pheno_model") && !inherits(model, "population_pheno_model")) {
    if (!length(model) || !all(vapply(model, inherits, logical(1), what = "pheno_model")))
      stop("model must be a specification or a non-empty list of pheno_model objects.",
           call. = FALSE)
    p <- if (missing(parameters)) lapply(model, model_parameters) else parameters
    if (is.numeric(p)) p <- rep(list(p), length(model))
    return(.predict_pheno_models(model, weather, p, stopatzc, basic_output))
  }
  if (is.list(weather) && !is.data.frame(weather)) {
    validate_parameters(model, parameters)
    return(.predict_pheno_seasons(model, weather, parameters, stopatzc, basic_output))
  }
  if (inherits(model, "population_pheno_model"))
    return(predict_population_phenology(model, weather, parameters,
      stopatzc = stopatzc, basic_output = basic_output))
  validate_parameters(model, parameters)
  .validate_sequential_weather(weather)
  inputs <- .calculate_pheno_inputs(model, weather, parameters)
  result <- .apply_pheno_structure(model, inputs$chill, inputs$heat, parameters,
                                  stopatzc, basic_output)
  if (basic_output)
    return(return_JDay(result$bloomindex, weather$JDay, weather$Year))
  result
}

.calculate_pheno_inputs <- function(model, weather, p) {
  p <- .normalize_characteristic_temperatures(p)
  
  if(model$chill$name == 'dynamic'){
    .dynamic_parameterization(names(p))
    p <- .convert_dynamic_parameters(p, "kinetic")
    
    chill_vec <- calculate_chill_dynamic(temp = weather$Temp, 
                                     times = seq_along(weather$Temp),
                                     E0 = p[["E0"]], E1 = p[["E1"]], 
                                     A0 = p[["A0"]], A1 = p[["A1"]],
                                     Tf = p[["Tf"]], slope = p[["slope"]])
  } else {
    times <- seq_along(weather$Temp)
    chill_vec <- switch(model$chill$name,
      chill_hours = calculate_chill_hours(weather$Temp, times,
        CH_lower = p[["CH_lower"]], CH_upper = p[["CH_upper"]]),
      utah = calculate_chill_utah(weather$Temp, times, weights = unname(p[paste0("utah_w", 1:7)])),
      positive_utah = calculate_chill_positive_utah(weather$Temp, times,
        weights = unname(p[paste0("pos_utah_w", 1:7)]),
        day_group = rep(seq_len(nrow(weather) / 24L), each = 24L)),
      none = calculate_chill_none(weather$Temp, times, heat_start_index = p[["heat_start_index"]]),
      custom = .calculate_custom_chill(model$chill, weather, p))
  }
  
  if(model$heat$name == 'gdh'){

    if(model$heat$scaling == "scaled"){
      heat_vec <- calculate_heat_gdh(temp = weather$Temp, 
                                     times = seq_along(weather$Temp), 
                                     Tb = p[["Tb"]],
                                     Tu = p[["Tu"]],
                                     Tc = p[["Tc"]])
    } else if(model$heat$scaling == "unscaled"){
      heat_vec <- calculate_heat_gdh_unscaled(temp = weather$Temp, 
                                              times = seq_along(weather$Temp), 
                                              Tb = p[["Tb"]],
                                              Tu = p[["Tu"]],
                                              Tc = p[["Tc"]])
    }
  }
  
  list(chill = chill_vec, heat = heat_vec)
}

.apply_pheno_structure <- function(model, chill_vec, heat_vec, p,
                                   stopatzc, basic_output) {
  if (model$chill$name == "none")
    return(apply_heat_only_structure(chill_vec, heat_vec, zc = p[["zc"]],
      stopatzc = stopatzc, basic_output = basic_output))
  #process the chill and heat with the model structure
  switch(
    model$structure,
    sequential = {
      apply_sequential_structure(chill = chill_vec, 
                                 heat = heat_vec, 
                                 yc = p[["yc"]], 
                                 zc = p[["zc"]], 
                                 stopatzc = stopatzc, 
                                 basic_output = basic_output)
    },
    parallel = {
      apply_parallel_structure(chill = chill_vec, 
                               heat = heat_vec, 
                               yc = p[["yc"]], 
                               zc = p[["zc"]], 
                               kmin = p[["kmin"]],
                               stopatzc = stopatzc, 
                               basic_output = basic_output)
    },
    partial_overlap = {
      apply_partial_overlap_structure(chill = chill_vec, 
                                      heat = heat_vec, 
                                      yc = p[["yc"]], 
                                      b1 = p[["b1"]], 
                                      b2 = p[["b2"]],
                                      b3 = p[["b3"]],
                                      ol = p[["ol"]],
                                      stopatzc = stopatzc, 
                                      basic_output = basic_output)
    },
    phenoflex = {
      apply_phenoflex_structure(chill = chill_vec, 
                                heat = heat_vec, 
                                yc = p[["yc"]], 
                                zc = p[["zc"]], 
                                s1 = p[["s1"]],
                                stopatzc = stopatzc, 
                                basic_output = basic_output)
    }
  )
}

.validate_sequential_weather <- function(weather) {
  required <- c("Temp", "Year", "JDay")
  if (!is.data.frame(weather) || !all(required %in% names(weather)) || anyDuplicated(names(weather)) ||
      nrow(weather) == 0L || nrow(weather) %% 24L != 0L ||
      !all(vapply(weather[required], function(x) is.numeric(x) && all(is.finite(x)), logical(1))))
    stop("weather must contain finite numeric Temp, Year and JDay for complete hourly days.", call. = FALSE)
  year <- weather$Year; day <- weather$JDay
  leap <- year %% 4 == 0 & (year %% 100 != 0 | year %% 400 == 0)
  if (any(year != floor(year) | year < 1 | year > 9999 | day != floor(day) | day < 1 | day > 365 + leap) ||
      any(weather$Temp <= -273))
    stop("weather contains invalid calendar values or temperatures <= -273 degrees C.", call. = FALSE)
  dates <- as.Date(sprintf("%04d-01-01", as.integer(year))) + day - 1
  starts <- seq.int(1L, nrow(weather), by = 24L)
  if (anyNA(dates) || !identical(dates, rep(dates[starts], each = 24)) ||
      any(diff(dates[starts]) != 1))
    stop("weather must have 24 chronological rows per day with no date gaps.", call. = FALSE)
  if ("Hour" %in% names(weather) &&
      (!is.numeric(weather$Hour) || anyNA(weather$Hour) ||
       !all(weather$Hour == rep(0:23, length(starts)))))
    stop("weather Hour must run from 0 to 23 each day.", call. = FALSE)
  invisible(TRUE)
}

