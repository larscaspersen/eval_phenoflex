.is_pheno_seasonlist <- function(weather) {
  is.list(weather) && !is.data.frame(weather) && length(weather) > 0L &&
    all(vapply(weather, is.data.frame, logical(1)))
}

.predict_pheno_seasons <- function(model, seasons, parameters, stopatzc, basic_output) {
  if (!.is_pheno_seasonlist(seasons))
    stop("weather must be a non-empty seasonlist of hourly weather data frames.",
         call. = FALSE)
  predict <- function(weather) predict_phenology(model, weather, parameters,
    stopatzc = stopatzc, basic_output = basic_output)
  if (basic_output && inherits(model, "pheno_model"))
    vapply(seasons, predict, numeric(1))
  else lapply(seasons, predict)
}

.predict_pheno_models <- function(models, weather, parameters, stopatzc, basic_output) {
  if (!is.list(parameters) || is.data.frame(parameters) ||
      length(parameters) != length(models))
    stop("parameters must supply one complete named vector per model, or a common vector.",
         call. = FALSE)
  for (i in seq_along(models)) validate_parameters(models[[i]], parameters[[i]])

  single_season <- is.data.frame(weather)
  labels <- names(models)
  if (single_season || .is_pheno_seasonlist(weather)) {
    weather <- rep(list(weather), length(models))
  } else {
    if (!is.list(weather) || length(weather) != length(models) ||
        !all(vapply(weather, .is_pheno_seasonlist, logical(1))))
      stop("weather must be a data frame, a non-empty seasonlist, or one non-empty seasonlist per model.",
           call. = FALSE)
    if (!is.null(names(weather))) labels <- names(weather)
  }
  result <- lapply(seq_along(models), function(i) {
    if (single_season)
      predict_phenology(models[[i]], weather[[i]], parameters[[i]],
        stopatzc = stopatzc, basic_output = basic_output)
    else .predict_pheno_seasons(models[[i]], weather[[i]], parameters[[i]],
                                stopatzc, basic_output)
  })
  names(result) <- labels
  result
}
