#' Evaluate combined cultivar fitting using residual sum of squares
#'
#' Select each cultivar's structure parameters, predict its seasons with shared
#' submodel parameters and sum the squared residuals across all cultivars and
#' seasons. Cultivars may have different numbers of seasons. Each observation
#' contributes equally; cultivar totals are summed without averaging.
#' Pass this function explicitly as `evaluation = phenology_rss_combined` to
#' [fit_phenology()]. A missing or non-finite prediction returns `Inf`.
#'
#' @param parameters Complete expanded named parameter vector.
#' @param model A [combined_pheno_model()] specification or a [pheno_model_list()]
#'   of ordinary cultivar/stage models.
#' @param seasons Nested list with one seasonlist per cultivar.
#' @param observed List with one numeric observation vector per cultivar/model.
#'   Each vector must match its seasonlist's length and order. Cultivars are
#'   matched by outer position, not by names. Dates use [predict_phenology()] units.
#'   Inner lists may contain NULL to skip individual season observations. At least
#'   one observed date is required across all models. Only observed pairs are predicted.
#' @return One numeric RSS value, or `Inf` for a missing or non-finite prediction.
#' @md
#' @export
phenology_rss_combined <- function(parameters, model, seasons, observed) {
  if (!inherits(model, "combined_pheno_model") && !inherits(model, "pheno_model_list"))
    stop("Use a combined_pheno_model or pheno_model_list specification.", call. = FALSE)
  validate_model_spec(model)
  .validate_combined_seasons(model, seasons)
  n <- if (inherits(model, "pheno_model_list")) length(model) else model$n_cultivars
  if (!is.list(observed) || is.data.frame(observed) || length(observed) != n)
    stop("observed must contain one numeric observation vector per cultivar.", call. = FALSE)
  dates <- lapply(seq_len(n), function(i) {
    .phenology_observations(observed[[i]], length(seasons[[i]]), paste0("Cultivar ", i, ": observed"))
  })
  if (!any(vapply(dates, function(x) any(!is.na(x)), logical(1))))
    stop("At least one observed phenology date is required.", call. = FALSE)
  simple <- as_pheno_models(model, parameters)
  sum(vapply(seq_len(n), function(i) {
    evaluate <- !is.na(dates[[i]])
    if (!any(evaluate)) return(0)
    phenology_rss(simple[[i]]$parameters, simple[[i]], seasons[[i]][evaluate], dates[[i]][evaluate])
  }, numeric(1)))
}
