#' Calibrate a list of ordinary phenology models
#'
#' This is a list of ordinary [pheno_model()] objects with attributes recording
#' which parameters are shared during calibration. Each `models[[i]]` remains
#' usable on its own. Unshared parameters receive indexed calibration names
#' (`yc1`, `yc2`, etc.); shared parameters retain their plain names.
#'
#' @param models Non-empty list of single models with identical structure, chill
#'   and heat specifications. Their stored parameter values may differ.
#' @param shared Names of shared parameters. NULL shares chill and heat submodel
#'   parameters and leaves structure parameters individual. Shared initial values
#'   must be equal across models. `character()` shares nothing.
#' @param ordered Names of unshared parameters whose values must be non-decreasing
#'   across list positions during fitting. Defaults to no order constraints.
#' @return A `pheno_model_list`, containing the original simple models.
#'   [model_parameters()] returns its complete calibration vector.
#'   [parameter_schema()] records base_name and member (NA for shared parameters).
#' @md
#' @export
pheno_model_list <- function(models, shared = NULL, ordered = character()) {
  if (!is.list(models) || is.data.frame(models) || !length(models) ||
      !all(vapply(models, inherits, logical(1), what = "pheno_model")))
    stop("models must be a non-empty list of single pheno_model objects.", call. = FALSE)
  models <- lapply(models, function(m) {
    m$parameters <- .normalize_characteristic_temperatures(m$parameters)
    m
  })
  schema <- parameter_schema(models[[1]])
  if (is.null(shared)) shared <- schema$name[schema$component != "structure"]
  class(models) <- c("pheno_model_list", "list")
  attr(models, "shared") <- shared
  attr(models, "ordered") <- ordered
  validate_model_spec(models)
  validate_parameters(models, model_parameters(models))
  models
}

.validate_pheno_model_list <- function(model) {
  if (!is.list(model) || is.data.frame(model) || !length(model) ||
      !all(vapply(model, inherits, logical(1), what = "pheno_model")))
    stop("Use a non-empty pheno_model_list of single models.", call. = FALSE)
  base <- model[[1]]
  for (m in model) {
    validate_parameters(m, m$parameters)
    chill <- m$chill; base_chill <- base$chill
    if (identical(chill$name, "dynamic")) {
      chill$characteristic <- NULL
      base_chill$characteristic <- NULL
    }
    if (!identical(m$structure, base$structure) || !identical(chill, base_chill) ||
        !identical(m$heat, base$heat))
      stop("All models must have identical structure, chill and heat specifications.", call. = FALSE)
  }
  parameter_names <- parameter_schema(base)$name
  shared <- attr(model, "shared")
  ordered <- attr(model, "ordered")
  for (setting in list(shared = shared, ordered = ordered)) {
    if (!is.character(setting) || anyNA(setting) || anyDuplicated(setting) ||
        !all(setting %in% parameter_names))
      stop("shared and ordered must contain distinct model parameter names.", call. = FALSE)
  }
  if (any(ordered %in% shared)) stop("Ordered parameters must be unshared.", call. = FALSE)
  for (name in shared) {
    values <- vapply(model, function(m) m$parameters[[name]], numeric(1))
    if (any(values != values[1]))
      stop("Shared parameter '", name, "' must have equal initial values across models.", call. = FALSE)
  }
  invisible(TRUE)
}

.pheno_model_list_schema <- function(model) {
  .validate_pheno_model_list(model)
  base <- parameter_schema(model[[1]])
  specific <- !base$name %in% attr(model, "shared")
  rows <- rep(seq_len(nrow(base)), ifelse(specific, length(model), 1L))
  schema <- base[rows, , drop = FALSE]
  schema$base_name <- schema$name
  schema$member <- unlist(lapply(seq_len(nrow(base)), function(i) {
    if (specific[i]) seq_along(model) else NA_integer_
  }), use.names = FALSE)
  indexed <- !is.na(schema$member)
  schema$name[indexed] <- paste0(schema$base_name[indexed], schema$member[indexed])
  rownames(schema) <- NULL
  schema
}

.model_list_member_parameters <- function(model, parameters, i) {
  base_names <- names(model[[i]]$parameters)
  keys <- ifelse(base_names %in% attr(model, "shared"), base_names, paste0(base_names, i))
  stats::setNames(unname(parameters[keys]), base_names)
}

.validate_model_list_parameters <- function(model, parameters) {
  representation <- .dynamic_parameterization(names(parameters),
    fallback = model[[1]]$chill$parameterization)
  if (!is.null(representation)) model <- .with_dynamic_parameterization(model, representation, use_defaults = TRUE)
  schema <- parameter_schema(model)
  if (!is.numeric(parameters) || !is.null(dim(parameters)) || is.null(names(parameters)) ||
      anyNA(names(parameters)) || anyDuplicated(names(parameters)) ||
      length(parameters) != nrow(schema) || !setequal(names(parameters), schema$name))
    stop("parameters must contain exactly the model-list schema names, without duplicates.", call. = FALSE)
  if (any(!is.finite(parameters))) stop("All parameter values must be finite.", call. = FALSE)
  parameters <- .normalize_characteristic_temperatures(parameters)
  for (i in seq_along(model)) {
    p <- .model_list_member_parameters(model, parameters, i)
    tryCatch(validate_parameters(model[[i]], p), error = function(e) {
      stop("Model ", i, ": ", conditionMessage(e), call. = FALSE)
    })
  }
  for (name in attr(model, "ordered")) {
    if (any(diff(parameters[paste0(name, seq_along(model))]) < 0))
      stop("Ordered parameter '", name, "' must be non-decreasing across models.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Get stored parameters for calibration
#'
#' @param model A single, population, combined, multistage or model-list specification, or a
#'   result from [fit_phenology()].
#' @return Complete named vector of stored kinetic parameters by default. For model
#'   lists this combines indexed member values and shared values.
#' @param parameterization Optional "characteristic" or "kinetic" representation
#'   for calibration starting values. NULL returns the stored kinetic values.
#' @md
#' @export
model_parameters <- function(model, parameterization = NULL) {
  if (inherits(model, "phenology_fit")) model <- model$model
  if (!is.null(parameterization)) model <- .with_dynamic_parameterization(model, parameterization)
  if (inherits(model, "stage_pheno_model")) model <- .stage_model_list(model)
  if (inherits(model, "population_pheno_model")) return(model$model$parameters)
  if (inherits(model, "pheno_model_list")) {
    schema <- parameter_schema(model)
    values <- vapply(seq_len(nrow(schema)), function(i) {
      member <- if (is.na(schema$member[i])) 1L else schema$member[i]
      model[[member]]$parameters[[schema$base_name[i]]]
    }, numeric(1))
    return(stats::setNames(values, schema$name))
  }
  if (inherits(model, "pheno_model") || inherits(model, "combined_pheno_model"))
    return(model$parameters)
  stop("Use a phenology model specification or fit result.", call. = FALSE)
}

#' Convert calibration specifications to ordinary single models
#'
#' Unpack a multistage specification into one ordinary model per stage, preserving
#' stage names and supplying each stage's heat requirement and all shared values.
#' Named combined specifications likewise retain their cultivar names.
#'
#' @param model A single model, [combined_pheno_model()], [stage_pheno_models()], [pheno_model_list()],
#'   or a fitted result from [fit_phenology()]. Population models are unsupported.
#' @param parameters Complete calibration vector, defaulting to stored/fitted values.
#'   Named theta_star and theta_c values in 0--20 are interpreted as Celsius
#'   using +273. Returned single models store kinetic coefficients, with supplied
#'   characteristic values retained in chill$characteristic in Kelvin.
#' @return A plain list of independent `pheno_model` objects, each storing a
#'   complete parameter vector including shared values. Calibration attributes
#'   are removed. Explicit parameters do not modify the supplied specification.
#' @md
#' @export
as_pheno_models <- function(model, parameters = model_parameters(model)) {
  force(parameters)
  if (inherits(model, "phenology_fit")) model <- model$model
  representation <- .dynamic_parameterization(names(parameters), fallback = NULL)
  if (!is.null(representation)) model <- .with_dynamic_parameterization(model, representation, use_defaults = TRUE)
  if (inherits(model, "stage_pheno_model")) model <- .stage_model_list(model)
  parameters <- .normalize_characteristic_temperatures(parameters)
  if (!inherits(model, "pheno_model") && !inherits(model, "combined_pheno_model") &&
      !inherits(model, "pheno_model_list"))
    stop("Use a single model, combined model or pheno_model_list.", call. = FALSE)
  validate_parameters(model, parameters)
  if (inherits(model, "pheno_model_list")) {
    result <- lapply(seq_along(model), function(i) {
      m <- model[[i]]
      .store_pheno_parameters(m, .model_list_member_parameters(model, parameters, i))
    })
    names(result) <- names(model)
    return(result)
  }
  if (inherits(model, "combined_pheno_model")) {
    base_names <- parameter_schema(model$model)$name
    result <- lapply(seq_len(model$n_cultivars), function(i) {
      m <- model$model
      .store_pheno_parameters(m, .combined_cultivar_parameters(model, parameters, i, base_names))
    })
    names(result) <- model$cultivar_names
    return(result)
  }
  list(.store_pheno_parameters(model, parameters))
}

#' Specify successive stages of one cultivar
#'
#' Store shared model settings and parameters once, with a named vector of
#' cumulative heat requirements in `stages` and `mode = "multistage"`.
#' All other parameters, including chill requirement, are shared.
#' Use [as_pheno_models()] to unpack the specification into ordinary single models.
#' Heat accumulation starts from the same season for every stage;
#' requirements are cumulative thresholds, not increments reset after each stage.
#' Requirements must remain non-decreasing during fitting.
#'
#' @param model Single [pheno_model()] containing a `zc` heat requirement
#'   (sequential, parallel or PhenoFlex). Its stored values initialize shared parameters.
#' @param heat_requirements Finite positive vector in stage order. Optional names
#'   label the stages and must be non-empty and distinct; otherwise labels are
#'   `stage1`, `stage2`, etc.
#' @return A `stage_pheno_model` S3 object containing `mode`, `structure`, `chill`,
#'   `heat`, shared `parameters` (excluding `zc`) and named `stages` heat requirements.
#'   Fitting and prediction accept this object directly. [model_parameters()]
#'   supplies calibration names `zc1`, `zc2`, etc., even for one stage.
#' @examples
#' models <- stage_pheno_models(pheno_model(), c(budbreak = 50, bloom = 100))
#' models
#' models$stages
#' individual_models <- as_pheno_models(models)
#' individual_models$bloom
#' model_parameters(models)
#' @md
#' @export
stage_pheno_models <- function(model = pheno_model(), heat_requirements) {
  if (!inherits(model, "pheno_model")) stop("Use a single pheno_model.", call. = FALSE)
  validate_parameters(model, model$parameters)
  if (!"zc" %in% names(model$parameters))
    stop("Stage heat requirements need a model with a zc parameter.", call. = FALSE)
  .validate_stage_heat_requirements(heat_requirements)
  if (is.null(names(heat_requirements)))
    names(heat_requirements) <- paste0("stage", seq_along(heat_requirements))
  model$parameters <- .normalize_characteristic_temperatures(model$parameters)
  result <- structure(list(mode = "multistage", structure = model$structure,
    chill = model$chill, heat = model$heat,
    parameters = model$parameters[setdiff(names(model$parameters), "zc")],
    stages = heat_requirements), class = "stage_pheno_model")
  validate_model_spec(result)
  result
}

.validate_stage_heat_requirements <- function(heat_requirements) {
  if (!is.numeric(heat_requirements) || !is.null(dim(heat_requirements)) ||
      !length(heat_requirements) || any(!is.finite(heat_requirements)) ||
      any(heat_requirements <= 0) || any(diff(heat_requirements) < 0))
    stop("heat_requirements must be positive, finite and non-decreasing in stage order.", call. = FALSE)
  labels <- names(heat_requirements)
  if (!is.null(labels) && (anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)))
    stop("Stage names must be non-empty and distinct.", call. = FALSE)
  invisible(TRUE)
}

# Adapt to the existing ordered-calibration machinery only when needed.
.stage_model_list <- function(model) {
  required <- c("mode", "structure", "chill", "heat", "parameters", "stages")
  if (!is.list(model) || !identical(sort(names(model)), sort(required)) ||
      !identical(model$mode, "multistage"))
    stop("Use a multistage specification from stage_pheno_models().", call. = FALSE)
  .validate_stage_heat_requirements(model$stages)
  if (is.null(names(model$stages)))
    stop("Multistage heat requirements must have stage names.", call. = FALSE)
  base <- structure(list(structure = model$structure, chill = model$chill,
    heat = model$heat), class = "pheno_model")
  parameter_names <- parameter_schema(base)$name
  shared <- setdiff(parameter_names, "zc")
  p <- model$parameters
  if (!"zc" %in% parameter_names)
    stop("Stage heat requirements need a model with a zc parameter.", call. = FALSE)
  if (!is.numeric(p) || !is.null(dim(p)) || is.null(names(p)) || anyNA(names(p)) ||
      anyDuplicated(names(p)) || length(p) != length(shared) || !setequal(names(p), shared))
    stop("Multistage parameters must contain exactly the shared schema names, excluding zc.", call. = FALSE)
  base$parameters <- c(p, zc = unname(model$stages[1]))
  base$parameters <- .normalize_characteristic_temperatures(base$parameters[parameter_names])
  models <- lapply(model$stages, function(zc) {
    m <- base
    m$parameters["zc"] <- zc
    m
  })
  # The calling schema, validation or prediction function validates this view.
  structure(models, class = c("pheno_model_list", "list"), shared = shared, ordered = "zc")
}

.stage_with_parameters <- function(model, parameters) {
  validate_parameters(model, parameters)
  parameters <- .normalize_characteristic_temperatures(parameters)
  model$parameters <- parameters[names(model$parameters)]
  model$stages[] <- parameters[paste0("zc", seq_along(model$stages))]
  model
}

#' @rdname stage_pheno_models
#' @param x A multistage model specification.
#' @param ... Additional arguments passed to print() for parameters and stages.
#' @export
print.stage_pheno_model <- function(x, ...) {
  cat("Mode:           multistage\n")
  print.pheno_model(x, ...)
  cat("Stages (cumulative heat requirements):\n")
  print(x$stages, ...)
  invisible(x)
}

.predict_pheno_model_list <- function(model, seasons, parameters, stopatzc, basic_output) {
  for (flag in list(stopatzc, basic_output)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("Output flags must be single non-missing logical values.", call. = FALSE)
  }
  simple <- as_pheno_models(model, parameters)
  .predict_pheno_models(simple, seasons, lapply(simple, model_parameters),
                        stopatzc, basic_output)
}
