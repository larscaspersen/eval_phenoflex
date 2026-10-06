# Shared scalar prediction contract for optimizer adapters.
.phenology_rss <- function(par, modelfn, bloomJDays, SeasonList, na_penalty) {
  if (!is.function(modelfn) || !is.list(SeasonList) || !length(SeasonList)) {
    stop("Supply a model function and a non-empty SeasonList.", call. = FALSE)
  }
  if (!is.numeric(bloomJDays) || length(bloomJDays) != length(SeasonList) ||
      any(!is.finite(bloomJDays))) {
    stop("bloomJDays must contain one finite observation per season.", call. = FALSE)
  }
  if (!is.numeric(na_penalty) || length(na_penalty) != 1L || !is.finite(na_penalty)) {
    stop("na_penalty must be a finite scalar.", call. = FALSE)
  }
  predictions <- vapply(SeasonList, function(season) {
    prediction <- modelfn(season, par = par)
    if (!is.numeric(prediction) && !identical(prediction, NA)) {
      stop("The model must return one numeric prediction per season.", call. = FALSE)
    }
    if (length(prediction) != 1L) {
      stop("The model must return one numeric prediction per season.", call. = FALSE)
    }
    if (is.na(prediction)) na_penalty else as.numeric(prediction)
  }, numeric(1))
  sum((predictions - bloomJDays)^2)
}
