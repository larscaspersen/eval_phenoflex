#' Evaluate phenology predictions using residual sum of squares
#'
#' Predict one phenology date per season with the candidate parameters and return
#' the sum of squared differences from the observed dates. This is a standalone
#' example evaluator; pass `evaluation = phenology_rss` to [fit_phenology()].
#' A missing or non-finite prediction
#' returns `Inf` so the optimizer excludes the candidate. Prediction errors
#' propagate.
#'
#' @param parameters Complete named numeric vector of candidate model parameters.
#' @param model A single [pheno_model()] specification. Population models require
#'   a custom evaluation function defining how bud predictions match observations.
#' @param seasons Non-empty list of hourly weather data frames accepted by
#'   [predict_phenology()], in the same order as `observed`.
#' @param observed Finite numeric vector with one observed phenology Julian day
#'   per season, or a list of finite numeric dates and NULL entries. A NULL
#'   skips that season. At least one observed date is required. Dates use the
#'   same units as [predict_phenology()].
#' @return One numeric RSS value to minimize, or `Inf` if a prediction is missing
#'   or non-finite. Observations are matched by position, not by name.
#' @examples
#' model <- pheno_model("sequential", chill_dynamic("kinetic"))
#' parameters <- model$parameters
#' parameters[c("yc", "zc")] <- c(0.1, 5)
#' weather <- data.frame(Temp = rep(8, 480), Year = 2008,
#'                       JDay = rep(50:69, each = 24))
#' observed <- predict_phenology(model, weather, parameters)
#' phenology_rss(parameters, model, seasons = list(weather), observed = observed)
#' @md
#' @export
phenology_rss <- function(parameters, model, seasons, observed) {
  if (!inherits(model, "pheno_model"))
    stop("phenology_rss requires a single pheno_model specification.", call. = FALSE)
  if (!is.list(seasons) || is.data.frame(seasons) || !length(seasons))
    stop("seasons must be a non-empty list of hourly weather data frames.", call. = FALSE)
  observed <- .phenology_observations(observed, length(seasons))
  evaluate <- !is.na(observed)
  if (!any(evaluate)) stop("At least one observed phenology date is required.", call. = FALSE)
  predicted <- vapply(seasons[evaluate], function(weather) {
    predict_phenology(model, weather, parameters)
  }, numeric(1))
  if (any(!is.finite(predicted))) return(Inf)
  sum((predicted - observed[evaluate])^2)
}

.phenology_observations <- function(observed, n, label = "observed") {
  invalid <- function() stop(label,
    " must contain one finite numeric phenology date per season, or a list with numeric dates and NULL for missing observations.",
    call. = FALSE)
  if (length(observed) != n || !is.null(dim(observed))) invalid()
  if (is.numeric(observed)) {
    if (any(!is.finite(observed))) invalid()
    return(observed)
  }
  if (!is.list(observed)) invalid()
  vapply(observed, function(x) {
    if (is.null(x)) return(NA_real_)
    if (!is.numeric(x) || length(x) != 1L || !is.null(dim(x)) || !is.finite(x)) invalid()
    as.numeric(x)
  }, numeric(1))
}
