#' Specify a phenology model
#'
#' A model specification is an immutable-by-convention S3 list describing the
#' algorithms, not their parameter values. This first implementation supports
#' sequential coupling with the Dynamic Model and GDH on hourly temperatures.
#' Other combinations are rejected as unsupported, not judged scientifically.
#' @param structure Model structure; currently only "sequential".
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
pheno_model <- function(structure = "sequential", chill = chill_dynamic(), heat = heat_gdh()) {
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
#' @export
heat_gdh <- function() list(name = "gdh")

#' @rdname pheno_model
#' @param model Model specification.
#' @export
validate_model_spec <- function(model) {
  if (!is.list(model) || !identical(sort(names(model)), sort(c("structure", "chill", "heat"))))
    stop("model must contain exactly structure, chill and heat.", call. = FALSE)
  if (!is.character(model$structure) ||
      length(model$structure) != 1L ||
      is.na(model$structure) ||
      !model$structure %in% c("sequential", "parallel")) {
    stop("Supported structures are sequential and parallel.")
  }
  if (!is.list(model$chill) ||
      !identical(sort(names(model$chill)), sort(c("name", "parameterization"))) ||
      !identical(model$chill$name, "dynamic") ||
      !(identical(model$chill$parameterization, "characteristic") ||
        identical(model$chill$parameterization, "kinetic")))
    stop("Use chill_dynamic() with characteristic or kinetic parameterization.", call. = FALSE)
  if (!identical(model$heat, heat_gdh()))
    stop("Only heat_gdh() is currently supported.", call. = FALSE)
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
  #check model
  validate_model_spec(model)
  
  # 1. Parameters belonging to the model structure
  structure_parameters <- switch(
    model$structure,
    sequential = c(yc = 40, zc = 190),
    parallel   = c(yc = 40, zc = 190, kmin = 0.1),
    stop("Unknown model structure.")
  )
  
  #check the chill submodel
  if(model$chill$name == 'dynamic'){
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
  } else {
    stop("Unknown heat submodel.")
  }
  
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
  schema <- parameter_schema(model)
  if (!is.numeric(parameters) || !is.null(dim(parameters)) || is.null(names(parameters)) ||
      anyNA(names(parameters)) || anyDuplicated(names(parameters)) ||
      length(parameters) != nrow(schema) || !setequal(names(parameters), schema$name))
    stop("parameters must be a named numeric vector with exactly the schema names, without duplicates.", call. = FALSE)
  if (any(!is.finite(parameters))) stop("All parameter values must be finite.", call. = FALSE)
  p <- parameters
  if (p[["yc"]] <= 0 || p[["zc"]] <= 0)
    stop("yc and zc must be positive.", call. = FALSE)
  
  if (model$structure == "parallel") {
    if (p[["kmin"]] < 0 || p[["kmin"]] > 1) {
      stop("kmin must be between 0 and 1.")
    }
  }
  
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
  }

  invisible(TRUE)
}

#' Predict one season with a model specification
#'
#' Adapts named parameters to the existing sequential wrapper without changing
#' its calculations. Characteristic parameters are converted with the shared
#' conversion implementation. Conversion failure raises an error; failure to
#' reach bloom returns NA. No calibration penalty is applied here.
#' @param model Model specification from pheno_model().
#' @param weather Data frame with Temp (degrees C), Year and JDay. Supply complete,
#' consecutive days, each with 24 hourly rows in chronological order. If Hour is
#' present it must be 0:23 for each day. Without Hour, within-day order is assumed.
#' @param parameters Named numeric vector; defaults to default_parameters(model).
#' @return One numeric fractional Julian bloom day, or NA if bloom is not reached.
#' Dates in the preceding year retain the legacy negative-day convention.
#' @export
predict_phenology <- function(model, weather, parameters = default_parameters(model)) {
  validate_parameters(model, parameters)
  .validate_sequential_weather(weather)
  p <- parameters
  
  if(model$chill$name == 'dynamic'){
    if (model$chill$parameterization == "characteristic") {
      chill <- p[c("theta_star", "theta_c", "tau", "pie_c")]
      converted <- characteristic_to_kinetic(c(rep(0, 4), unname(chill), rep(0, 4)), failure_return = "NA")
      if (any(!is.finite(converted[5:8])))
        stop("Characteristic-to-kinetic conversion failed for these parameters.", call. = FALSE)
      p <- c(p, stats::setNames(converted[5:8], c("E0", "E1", "A0", "A1")))
    }
    
    chill_order <- c("E0", "E1", "A0", "A1", "Tf", "slope")
  }
  
  if(model$heat$name == 'gdh'){
    heat_order <- c( "Tb", "Tu", "Tc")
  }
  
  order_chill_heat <- c(chill_order, heat_order)
  
  switch(
    model$structure,
    sequential = {
      order <- c("yc", "zc", order_chill_heat)
      wrapper_seq_model(weather, unname(p[order]))
    },
    parallel = {
      order <- c("yc", "zc", "kmin", order_chill_heat)
      wrapper_parallel_model(weather, unname(p[order]))
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

