.dynamic_kinetic_names <- c("E0", "E1", "A0", "A1")
.dynamic_characteristic_names <- c("theta_star", "theta_c", "tau", "pie_c")
.dynamic_kinetic_defaults <- c(E0 = 4153.5, E1 = 12888.8, A0 = 139500, A1 = 2.567e18)
.dynamic_characteristic_defaults <- c(theta_star = 279, theta_c = 286.1, tau = 47.7, pie_c = 28)

# Inspect all supplied vectors together, including indexed collection names.
.dynamic_parameterization <- function(..., fallback = "kinetic") {
  keys <- unlist(list(...), use.names = FALSE)
  characteristic <- any(grepl("^(theta_star|theta_c|tau|pie_c)[0-9]*$", keys))
  kinetic <- any(grepl("^(E0|E1|A0|A1)[0-9]*$", keys))
  if (characteristic && kinetic)
    stop("Cannot mix characteristic and kinetic Dynamic Model parameters.", call. = FALSE)
  if (characteristic) "characteristic" else if (kinetic) "kinetic" else fallback
}

.convert_dynamic_parameters <- function(parameters, to) {
  from <- if (to == "kinetic") .dynamic_characteristic_names else .dynamic_kinetic_names
  target <- if (to == "kinetic") .dynamic_kinetic_names else .dynamic_characteristic_names
  if (!all(from %in% names(parameters))) return(parameters)
  input <- c(rep(0, 4), unname(parameters[from]), rep(0, 4))
  converted <- if (to == "kinetic") characteristic_to_kinetic(input, failure_return = "NA") else
    kinetic_to_characteristic(input)
  if (any(!is.finite(converted[5:8])))
    stop(structure(list(message = "Characteristic-to-kinetic conversion failed for these parameters.",
      call = NULL), class = c("phenology_conversion_error", "error", "condition")))
  positions <- match(from, names(parameters))
  parameters[positions] <- converted[5:8]
  names(parameters)[positions] <- target
  parameters
}

# Calibration views are temporary; public constructors and fitted models are kinetic.
.with_dynamic_parameterization <- function(model, parameterization, use_defaults = FALSE) {
  parameterization <- match.arg(parameterization, c("kinetic", "characteristic"))
  if (inherits(model, "phenology_fit")) model <- model$model
  if (is.list(model) && !inherits(model, "pheno_model_list") && length(model) &&
      all(vapply(model, inherits, logical(1), what = "pheno_model")))
    return(lapply(model, .with_dynamic_parameterization, parameterization = parameterization,
      use_defaults = use_defaults))
  if (inherits(model, "pheno_model_list")) {
    from <- if (parameterization == "kinetic") .dynamic_characteristic_names else .dynamic_kinetic_names
    to <- if (parameterization == "kinetic") .dynamic_kinetic_names else .dynamic_characteristic_names
    if (identical(model[[1]]$chill$name, "dynamic") &&
        !identical(model[[1]]$chill$parameterization, parameterization)) {
      shared <- attr(model, "shared")
      if (any(from %in% shared) && !all(from %in% shared))
        stop("Changing Dynamic Model parameterization requires sharing all four chill coefficients or none.", call. = FALSE)
      if (any(from %in% attr(model, "ordered")))
        stop("Ordered chill coefficients cannot be converted between parameterizations.", call. = FALSE)
      shared[match(from, shared, nomatch = 0L)] <- to[from %in% shared]
      attr(model, "shared") <- shared
    }
    for (i in seq_along(model)) model[[i]] <- .with_dynamic_parameterization(model[[i]], parameterization, use_defaults)
    return(model)
  }
  if (inherits(model, "population_pheno_model")) {
    model$model <- .with_dynamic_parameterization(model$model, parameterization, use_defaults)
    return(model)
  }
  if (inherits(model, "combined_pheno_model")) {
    old <- model$model
    model$model <- .with_dynamic_parameterization(old, parameterization, use_defaults)
    if (!identical(old$chill$parameterization, model$model$chill$parameterization) &&
        !is.null(model$parameters)) {
      from <- if (parameterization == "kinetic") .dynamic_characteristic_names else .dynamic_kinetic_names
      to <- if (parameterization == "kinetic") .dynamic_kinetic_names else .dynamic_characteristic_names
      if (use_defaults || identical(unname(model$parameters[from]), unname(old$parameters[from]))) {
        positions <- match(from, names(model$parameters))
        model$parameters[positions] <- model$model$parameters[to]
        names(model$parameters)[positions] <- to
      } else model$parameters <- .convert_dynamic_parameters(model$parameters, parameterization)
    }
    return(model)
  }
  if (!identical(model$chill$name, "dynamic")) return(model)
  if (!identical(model$chill$parameterization, parameterization)) {
    if (!is.null(model$parameters)) {
      p <- .normalize_characteristic_temperatures(model$parameters)
      if (use_defaults) {
        from <- if (parameterization == "kinetic") .dynamic_characteristic_names else .dynamic_kinetic_names
        values <- if (parameterization == "kinetic") .dynamic_kinetic_defaults else .dynamic_characteristic_defaults
        positions <- match(from, names(p))
        p[positions] <- values
        names(p)[positions] <- names(values)
        model$chill$characteristic <- NULL
      } else if (parameterization == "kinetic") {
        model$chill$characteristic <- p[.dynamic_characteristic_names]
        p <- .convert_dynamic_parameters(p, "kinetic")
      } else {
        cached <- model$chill$characteristic
        # Reuse exact supplied values only while they describe the stored kinetics.
        equivalent <- !is.null(cached) && isTRUE(all.equal(
          unname(.convert_dynamic_parameters(cached, "kinetic")),
          unname(p[.dynamic_kinetic_names]), tolerance = 1e-10))
        p <- if (equivalent) {
          positions <- match(.dynamic_kinetic_names, names(p))
          p[positions] <- cached
          names(p)[positions] <- .dynamic_characteristic_names
          p
        } else .convert_dynamic_parameters(p, "characteristic")
      }
      model$parameters <- p
    }
    model$chill$parameterization <- parameterization
  }
  model
}

.store_pheno_parameters <- function(model, parameters) {
  parameters <- .normalize_characteristic_temperatures(parameters)
  representation <- .dynamic_parameterization(names(parameters), fallback = model$chill$parameterization)
  if (identical(model$chill$name, "dynamic")) {
    model$chill$parameterization <- representation
    # Cached characteristic values must never describe an older fitted candidate.
    cached <- model$chill$characteristic
    if (representation == "characteristic" || is.null(cached) || !isTRUE(all.equal(
        unname(.convert_dynamic_parameters(cached, "kinetic")),
        unname(parameters[.dynamic_kinetic_names]), tolerance = 1e-10)))
      model$chill$characteristic <- NULL
  }
  model$parameters <- parameters
  .with_dynamic_parameterization(model, "kinetic")
}

.canonical_pheno_model <- function(model) {
  if (inherits(model, "stage_pheno_model") && identical(model$chill$name, "dynamic") &&
      identical(model$chill$parameterization, "kinetic"))
    return(.store_pheno_parameters(model, model$parameters))
  .with_dynamic_parameterization(model, "kinetic")
}
