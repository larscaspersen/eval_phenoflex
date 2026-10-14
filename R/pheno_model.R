#' Specify a phenology model
#'
#' A model specification is an S3 list storing algorithms and named parameter
#' values. It supports sequential, linear parallel,
#' partial-overlap and PhenoFlex coupling with the Dynamic Model, chill hours,
#' Utah or positive Utah chill and GDH heat. With chill_none(), only heat is
#' accumulated and the chosen coupling is bypassed.
#' Printing shows the structure, submodel names and parameters in a compact form.
#' @param structure One of "sequential", "parallel", "partial_overlap" or "phenoflex".
#' @param chill Chill specification returned by chill_dynamic(), chill_hours(),
#'   chill_utah(), chill_positive_utah(), chill_none() or chill_custom().
#' @param heat Heat specification returned by heat_gdh().
#' @param parameters Named numeric vector specifying any subset of the model's
#' parameter names. Omitted parameters use default_parameters() for the selected
#' structure, chill and heat specifications. NULL uses all defaults. The complete
#' vector is validated after filling defaults.
#' Dynamic Model input names select the representation automatically: theta_star,
#' theta_c, tau or pie_c select characteristic input; E0, E1, A0 or A1 select
#' kinetic input. Mixing the two sets is an error. With neither set supplied,
#' the original kinetic coefficients are used. Characteristic inputs are converted
#' at construction and retained in chill$characteristic; parameters always stores
#' kinetic coefficients. Missing characteristic values use characteristic defaults.
#' For theta_star and theta_c, each value in 0--20 is interpreted as Celsius
#' and converted using the native kernels' legacy offset K = Celsius + 273.
#' Other values are interpreted as Kelvin. Stored values are always Kelvin;
#' all other temperature parameters remain Celsius.
#' @return pheno_model() returns a pheno_model list. Component constructors return
#' lists; validate_model_spec() invisibly returns TRUE or raises an error.
#' @examples
#' model <- pheno_model()
#' parameters <- default_parameters(model)
#' parameters["yc"] <- 45
#' model <- pheno_model(parameters = parameters)
#' model$parameters
#' parallel <- pheno_model("parallel", parameters = c(yc = 50, zc = 180))
#' model <- pheno_model(parameters = c(theta_star = 6, theta_c = 13.1))
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
#' @return parameter_schema() returns a data frame of names, components, defaults
#' and representations. Combined models also include base_name and cultivar
#' (NA for shared parameters); model lists use member instead of cultivar.
#' default_parameters() returns a named numeric vector.
#' validate_parameters() invisibly returns TRUE or raises an error.
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
    
    representation <- model$chill$parameterization
    chill <- if (representation == "characteristic")
      c(.dynamic_characteristic_defaults, Tf = 4, slope = 1.6)
    else c(.dynamic_kinetic_defaults, Tf = 4, slope = 1.6)
  
    # Tf and slope are native parameters in both representations.
    chill_representation <- rep("native", length(chill))
    chill_representation[1:4] <- representation
  } else {
    chill <- switch(model$chill$name,
      chill_hours = c(CH_lower = 0, CH_upper = 7.2),
      utah = stats::setNames(c(0, 0.5, 1, 0.5, 0, -0.5, -1), paste0("utah_w", 1:7)),
      positive_utah = stats::setNames(c(0, 0.5, 1, 0.5, 0, -0.5, -1), paste0("pos_utah_w", 1:7)),
      none = c(heat_start_index = 1), custom = model$chill$defaults)
    chill_representation <- rep("native", length(chill))
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
  
  representation <- c(
    rep("native", length(structure_parameters)),
    chill_representation,
    rep("native", length(heat))
  )
  
  
  data.frame(
    name = names(values),
    component = component,
    default = unname(values),
    representation = representation,
    stringsAsFactors = FALSE
  )
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

