#' Evaluate successive phenological stages for one cultivar
#'
#' All stage models use the same seasonlist. Observations are ordered by stage,
#' then by season. A NULL observation skips only that stage/season pair; its
#' weather may still be evaluated for other stages. Missing predictions are
#' penalized only where observations exist. Pass this evaluator explicitly to
#' [fit_phenology()].
#'
#' @param parameters Complete model-list calibration vector.
#' @param model A model list returned by [stage_pheno_models()].
#' @param seasons Non-empty seasonlist of hourly weather data frames.
#' @param observed List with one entry per stage in heat-requirement order. Each
#'   entry is a numeric vector or a list as long as `seasons`. In an inner list,
#'   entries must be finite numeric dates or NULL. A fully unobserved stage can
#'   use `rep(list(NULL), length(seasons))`. At least one observation is required
#'   across the whole dataset. Names do not alter positional matching.
#' @return Total numeric RSS over observed stage/season pairs, or `Inf` for a
#'   non-finite prediction at an observed pair.
#' @md
#' @export
phenology_rss_stages <- function(parameters, model, seasons, observed) {
  if (!inherits(model, "pheno_model_list"))
    stop("Use stage_pheno_models() to create a list of simple stage models.", call. = FALSE)
  if (!is.list(seasons) || is.data.frame(seasons) || !length(seasons) ||
      !all(vapply(seasons, is.data.frame, logical(1))))
    stop("seasons must be a non-empty seasonlist of weather data frames.", call. = FALSE)
  nested_seasons <- rep(list(seasons), length(model))
  phenology_rss_combined(parameters, model, nested_seasons, observed)
}
