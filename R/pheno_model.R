#' Specify a phenology model
#'
#' A model specification is an immutable-by-convention S3 list describing the
#' algorithms, not their parameter values. Supports sequential, linear parallel,
#' partial-overlap and PhenoFlex coupling with the Dynamic Model and GDH.
#' Other combinations are rejected as unsupported, not judged scientifically.
#' @param structure One of "sequential", "parallel", "partial_overlap" or "phenoflex".
#' @param chill Chill specification returned by chill_dynamic().
#' @param heat Heat specification returned by heat_gdh().
#' @return pheno_model() returns a pheno_model list. Component constructors return
#' lists; validate_model_spec() invisibly returns TRUE or raises an error.
#' @examples
#' model <- pheno_model()
#' parameters <- default_parameters(model)
#' parameters["yc"] <- 45
#' validate_parameters(model, parameters)
#' @export
pheno_model <- function(structure = "phenoflex", chill = chill_dynamic(), heat = heat_gdh()) {
  model <- list(structure = structure, chill = chill, heat = heat)
  validate_model_spec(model)
  class(model) <- "pheno_model"
  model
}

#' @rdname pheno_model
#' @param parameterization Dynamic Model representation: "characteristic" or "kinetic".
#' @export
chill_dynamic <- function(parameterization = c("characteristic", "kinetic")) {
  list(name = "dynamic", parameterization = match.arg(parameterization))
}

#' @rdname pheno_model
#' @param scaling "unscaled" returns Growing Degree Hour (GDH) units without scaling via (Tu-Tb).
#' "scaled" returns units multiplied by (Tu-Tb). Standard implementation of PhenoFlex
#' uses unscaled GDH, while most other studies use scaled GDH.
#' @export
heat_gdh <- function(scaling = c("unscaled", "scaled")) list(name = "gdh",
                            scaling = match.arg(scaling))

#' @rdname pheno_model
#' @param model Model specification.
#' @export
validate_model_spec <- function(model) {
  if (inherits(model, "population_pheno_model"))
    return(validate_model_spec(model$model))
  if (!is.list(model) || !identical(sort(names(model)), sort(c("structure", "chill", "heat"))))
    stop("model must contain exactly structure, chill and heat.", call. = FALSE)
  if (!is.character(model$structure) ||
      length(model$structure) != 1L ||
      is.na(model$structure) ||
      !model$structure %in% c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    stop("Supported structures are sequential, parallel, partial_overlap and phenoflex.")
  }
  if (!is.list(model$chill) ||
      !identical(sort(names(model$chill)), sort(c("name", "parameterization"))) ||
      !identical(model$chill$name, "dynamic") ||
      !(identical(model$chill$parameterization, "characteristic") ||
        identical(model$chill$parameterization, "kinetic")))
    stop("Use chill_dynamic() with characteristic or kinetic parameterization.", call. = FALSE)

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
#' @param model Model specification from pheno_model().
#' @param parameters Named numeric vector with exactly the names in the schema;
#' order is arbitrary. Matrices and arrays are not parameter vectors.
#' @return parameter_schema() returns a data frame of names, components, defaults
#' and representations. default_parameters() returns a named numeric vector.
#' validate_parameters() invisibly returns TRUE or raises an error.
#' @export
parameter_schema <- function(model) {
  if (inherits(model, "population_pheno_model"))
    return(parameter_schema(model$model))
  #check model
  validate_model_spec(model)
  

  #check the chill submodel
  if(model$chill$name == 'dynamic'){
    
    #typical chill requirement 
    yc = 40
    
    representation <- model$chill$parameterization
    chill <- if (representation == "characteristic")
      c(theta_star = 279, theta_c = 286.1, tau = 47.7, pie_c = 28, Tf = 4, slope = 1.6)
    else c(E0 = 4153.5, E1 = 12888.8, A0 = 139500, A1 = 2.567e18, Tf = 4, slope = 1.6)
  
    # Tf and slope are native parameters in both representations.
    chill_representation <- rep("native", length(chill))
    chill_representation[1:4] <- representation
  } else {
    stop("Unknown chill submodel.")
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
  structure_parameters <- switch(
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
default_parameters <- function(model) {
  schema <- parameter_schema(model)
  stats::setNames(schema$default, schema$name)
}

#' @rdname parameter_schema
#' @export
validate_parameters <- function(model, parameters) {
  if (inherits(model, "population_pheno_model"))
    return(validate_parameters(model$model, parameters))
  schema <- parameter_schema(model)
  if (!is.numeric(parameters) || !is.null(dim(parameters)) || is.null(names(parameters)) ||
      anyNA(names(parameters)) || anyDuplicated(names(parameters)) ||
      length(parameters) != nrow(schema) || !setequal(names(parameters), schema$name))
    stop("parameters must be a named numeric vector with exactly the schema names, without duplicates.", call. = FALSE)
  if (any(!is.finite(parameters))) stop("All parameter values must be finite.", call. = FALSE)
  p <- parameters
  switch(
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
  }
  
  if(model$heat$name == 'gdh'){
    if (!(p[["Tb"]] < p[["Tu"]] && p[["Tu"]] < p[["Tc"]]))
      stop("Heat parameters must satisfy Tb < Tu < Tc.", call. = FALSE)
    
    if (p[["Tb"]] <= -273)
      stop("Tb must exceed -273 degrees C.", call. = FALSE)
  }

  invisible(TRUE)
}

#' Predict one season with a model specification
#'
#' Calculates shared chill and potential heat and applies the selected C++
#' structure. Characteristic parameters use the shared conversion implementation.
#' Conversion failure raises an error; no bloom is represented by index zero.
#' Population specifications dispatch to [predict_population_phenology()].
#' @param model Model specification from pheno_model().
#' @param weather Data frame with Temp (degrees C), Year and JDay. Supply complete,
#' consecutive days, each with 24 hourly rows in chronological order. If Hour is
#' present it must be 0:23 for each day. Without Hour, within-day order is assumed.
#' @param parameters Named numeric vector; defaults to default_parameters(model).
#' @param stopatzc Stop at the first bloom if TRUE. FALSE completes trajectories
#' while retaining the first bloom index.
#' @param basic_output boolean. If 'TRUE', only the bloomindex is returned as 
#' a named element of the return list.
#' @return For a single model, a list with `bloomindex` (one-based weather row;
#' zero for no bloom). Detailed output adds `chill` and `z`. Use [return_JDay()]
#' to convert the index to a fractional Julian day. Population specifications
#' return the result described in [predict_population_phenology()].
#' @export
predict_phenology <- function(model, weather, 
                              parameters = default_parameters(model),
                              stopatzc = TRUE, 
                              basic_output = TRUE) {
  if (inherits(model, "population_pheno_model"))
    return(predict_population_phenology(model, weather, parameters,
      stopatzc = stopatzc, basic_output = basic_output))
  validate_parameters(model, parameters)
  .validate_sequential_weather(weather)
  inputs <- .calculate_pheno_inputs(model, weather, parameters)
  .apply_pheno_structure(model, inputs$chill, inputs$heat, parameters,
                         stopatzc, basic_output)
}

.calculate_pheno_inputs <- function(model, weather, p) {
  
  if(model$chill$name == 'dynamic'){
    if (model$chill$parameterization == "characteristic") {
      chill <- p[c("theta_star", "theta_c", "tau", "pie_c")]
      converted <- characteristic_to_kinetic(c(rep(0, 4), unname(chill), rep(0, 4)), failure_return = "NA")
      if (any(!is.finite(converted[5:8])))
        stop("Characteristic-to-kinetic conversion failed for these parameters.", call. = FALSE)
      p <- c(p, stats::setNames(converted[5:8], c("E0", "E1", "A0", "A1")))
    }
    
    chill_vec <- calculate_chill_dynamic(temp = weather$Temp, 
                                     times = seq_along(weather$Temp),
                                     E0 = p[["E0"]], E1 = p[["E1"]], 
                                     A0 = p[["A0"]], A1 = p[["A1"]],
                                     Tf = p[["Tf"]], slope = p[["slope"]])
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

