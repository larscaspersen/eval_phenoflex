#' Specify and check a user-defined chill submodel
#'
#' Supply a deterministic R callback accepting `weather` and `parameters` by
#' name (or through `...`). It receives the hourly weather data frame and only
#' the named parameters declared in `defaults`. Return a finite numeric matrix
#' with one row per weather timestamp and at least one column. The final column
#' is the non-negative chill state consumed by the selected model structure;
#' preceding columns may contain signed diagnostics. Column names are optional,
#' but must be non-empty and distinct if supplied. Unnamed columns receive names.
#'
#' Row i describes the state at the start of interval i. The final weather row
#' is an endpoint, so there are n-1 temperature updates for n output rows.
#' An initial chill state is permitted; the state need not increase monotonically.
#'
#' Construction checks the metadata, validates default parameters and calls the
#' callback twice on a small varied hourly sample (or `check_weather`). Outputs
#' must be identical for the same inputs. [validate_chill_custom()] repeats this
#' check for supplied weather/parameters. Every prediction also checks output
#' shape, finiteness and the coupling state. These checks verify the interface,
#' not biological validity or behaviour on every possible weather series.
#'
#' Custom parameters enter [parameter_schema()] as chill parameters and are
#' shared by default in combined, stage and population models. Names must not
#' collide with structure/heat parameters, their indexed forms, or the package's
#' characteristic-temperature and heat-start names. Constructor bounds are
#' calibration search bounds, not parameter-domain constraints. Parameters
#' without supplied bounds stay fixed at their stored values in [default_bounds()].
#'
#' Use `validate` to reject parameter-domain violations. Invalid candidates are
#' excluded during fitting. Callback errors and malformed trajectories stop
#' fitting with an informative error instead of being silently treated as loss.
#' For parallel fitting, capture needed constants/helpers in the callback's
#' closure and use namespace-qualified calls to external packages.
#' @param fun Callback accepting weather and parameters.
#' @param defaults Finite numeric vector with distinct, non-empty parameter names.
#'   `numeric()` supports a callback with no adjustable chill parameters.
#' @param lower,upper Optional finite named bounds for the same subset of defaults.
#'   Supply both or neither. Defaults must lie within the supplied bounds.
#' @param validate Optional function accepting a parameter vector, returning TRUE
#'   or NULL on success, and FALSE or an error on invalid parameter values.
#' @param label Non-empty display label for compact model printing.
#' @param check_weather Optional hourly data frame for constructor checks. Supply
#'   this when the callback needs additional weather columns. It is not stored.
#' @param chill A specification returned by chill_custom().
#' @param weather Hourly weather data frame, or NULL for the built-in check sample.
#' @param parameters Complete vector of the custom chill parameters to check.
#' @return chill_custom() returns a custom chill specification.
#'   validate_chill_custom() invisibly returns TRUE or raises an error.
#' @examples
#' triangular <- function(weather, parameters) {
#'   temp <- head(weather$Temp, -1L)
#'   increments <- pmax(0, 1 - abs(temp - parameters[["Topt"]]) / parameters[["width"]])
#'   matrix(c(0, cumsum(increments)), ncol = 1, dimnames = list(NULL, "y"))
#' }
#' custom <- chill_custom(
#'   triangular, defaults = c(Topt = 5, width = 5),
#'   lower = c(Topt = 0, width = 1), upper = c(Topt = 10, width = 15),
#'   validate = function(p) p[["width"]] > 0, label = "Triangular chill"
#' )
#' model <- pheno_model(chill = custom)
#' model
#' parameter_schema(model)
#' @md
#' @export
chill_custom <- function(fun, defaults = numeric(), lower = NULL, upper = NULL,
                         validate = NULL, label = "Custom chill", check_weather = NULL) {
  chill <- list(name = "custom", fun = fun, defaults = defaults,
    lower = lower, upper = upper, validate = validate, label = label)
  validate_chill_custom(chill, weather = check_weather)
  chill
}

.custom_chill_numeric <- function(x, what, allowed = NULL) {
  if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x)) ||
      (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
        any(!nzchar(trimws(names(x)))) || anyDuplicated(names(x)))) ||
      (!is.null(allowed) && !all(names(x) %in% allowed)))
    stop(what, " must be a finite named numeric vector with distinct, non-empty parameter names.",
         call. = FALSE)
  invisible(TRUE)
}

.validate_custom_chill_spec <- function(chill) {
  fields <- c("name", "fun", "defaults", "lower", "upper", "validate", "label")
  if (!is.list(chill) || !identical(sort(names(chill)), sort(fields)) ||
      !identical(chill$name, "custom"))
    stop("Use chill_custom() to construct a custom chill specification.", call. = FALSE)
  if (!is.function(chill$fun)) stop("Custom chill fun must be a function.", call. = FALSE)
  arguments <- names(formals(chill$fun))
  if (!all(c("weather", "parameters") %in% arguments) && !"..." %in% arguments)
    stop("Custom chill fun must accept weather and parameters by name, or through ... .", call. = FALSE)
  if (!is.character(chill$label) || length(chill$label) != 1L || is.na(chill$label) ||
      !nzchar(trimws(chill$label))) stop("Custom chill label must be a non-empty string.", call. = FALSE)
  .custom_chill_numeric(chill$defaults, "Custom chill defaults")
  reserved <- "^(yc|zc|s1|kmin|b1|b2|b3|ol|chill_scale|Tb|Tu|Tc|theta_star|theta_c|heat_start_index)[0-9]*$"
  if (any(grepl(reserved, names(chill$defaults))))
    stop("Custom chill parameter names must not collide with reserved structure, heat or temperature names.",
         call. = FALSE)
  if (xor(is.null(chill$lower), is.null(chill$upper)))
    stop("Supply both lower and upper custom chill bounds, or neither.", call. = FALSE)
  if (!is.null(chill$lower)) {
    .custom_chill_numeric(chill$lower, "Custom chill lower bounds", names(chill$defaults))
    .custom_chill_numeric(chill$upper, "Custom chill upper bounds", names(chill$defaults))
    if (!setequal(names(chill$lower), names(chill$upper)))
      stop("Custom chill lower and upper must name the same subset of defaults.", call. = FALSE)
    keys <- names(chill$lower)
    if (!all(keys %in% names(chill$defaults)))
      stop("Custom chill bounds must use declared default names.", call. = FALSE)
    if (any(chill$lower > chill$upper[keys]))
      stop("Custom chill lower bounds must not exceed upper bounds.", call. = FALSE)
    if (any(chill$defaults[keys] < chill$lower | chill$defaults[keys] > chill$upper[keys]))
      stop("Custom chill defaults must lie within lower and upper bounds.", call. = FALSE)
  }
  if (!is.null(chill$validate) && !is.function(chill$validate))
    stop("Custom chill validate must be a function or NULL.", call. = FALSE)
  invisible(TRUE)
}

.validate_custom_chill_parameters <- function(chill, parameters) {
  .custom_chill_numeric(parameters, "Custom chill parameters", names(chill$defaults))
  if (length(parameters) != length(chill$defaults) ||
      !setequal(names(parameters), names(chill$defaults)))
    stop("Custom chill parameters must contain exactly the declared default names.", call. = FALSE)
  if (!is.null(chill$validate)) {
    valid <- tryCatch(chill$validate(parameters[names(chill$defaults)]), error = function(e)
      stop("Custom chill '", chill$label, "' parameter validation failed: ", conditionMessage(e), call. = FALSE))
    if (identical(valid, FALSE))
      stop("Custom chill '", chill$label, "' parameters were rejected by its validator.", call. = FALSE)
    if (!is.null(valid) && !isTRUE(valid))
      stop("Custom chill '", chill$label, "' parameter validator must return TRUE or NULL for valid parameters.",
           call. = FALSE)
  }
  invisible(TRUE)
}

.custom_chill_check_weather <- function() {
  data.frame(Temp = rep(c(-5, 0, 4, 5, 7.2, 8, 12, 20, 25, 30, 36, 40), 6),
    Year = 2008, JDay = rep(1:3, each = 24), Hour = rep(0:23, 3))
}

.calculate_custom_chill <- function(chill, weather, parameters) {
  result <- tryCatch(chill$fun(weather = weather, parameters = parameters[names(chill$defaults)]),
    error = function(e) stop("Custom chill '", chill$label, "' callback failed: ", conditionMessage(e), call. = FALSE))
  prefix <- paste0("Custom chill '", chill$label, "': ")
  if (!is.matrix(result) || !is.numeric(result) || nrow(result) != nrow(weather) || ncol(result) < 1L)
    stop(prefix, "return a numeric matrix with one row per weather timestamp and at least one column.", call. = FALSE)
  if (any(!is.finite(result))) stop(prefix, "all trajectory values must be finite.", call. = FALSE)
  if (any(result[, ncol(result)] < 0))
    stop(prefix, "the final coupling column must be non-negative; put signed diagnostics in earlier columns.", call. = FALSE)
  labels <- colnames(result)
  if (is.null(labels)) {
    labels <- paste0("diagnostic", seq_len(ncol(result)))
    labels[length(labels)] <- "y"
    colnames(result) <- labels
  } else if (anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)) {
    stop(prefix, "column names must be non-empty and distinct.", call. = FALSE)
  }
  storage.mode(result) <- "double"
  result
}

#' @rdname chill_custom
#' @export
validate_chill_custom <- function(chill, weather = NULL, parameters = chill$defaults) {
  .validate_custom_chill_spec(chill)
  .validate_custom_chill_parameters(chill, parameters)
  if (is.null(weather)) weather <- .custom_chill_check_weather()
  .validate_sequential_weather(weather)
  first <- .calculate_custom_chill(chill, weather, parameters)
  second <- .calculate_custom_chill(chill, weather, parameters)
  if (!identical(first, second))
    stop("Custom chill '", chill$label, "' callback must be deterministic for the same weather and parameters.",
         call. = FALSE)
  invisible(TRUE)
}
