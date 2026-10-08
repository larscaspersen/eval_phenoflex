#' Construct an ensemble of fitted phenology models
#'
#' Separate independently calibrated fits from pheno_model_list(), which records
#' shared parameters within one calibration. Each ensemble entry can itself be a
#' combined/stage fit; predictions are aggregated separately for each output.
#'
#' @param models A phenology_cv result, or a non-empty list of phenology_fit or
#'   pheno_model/combined_pheno_model/pheno_model_list objects. Entries must have
#'   the same number and order of cultivar/stage outputs. Names identify ensemble
#'   members, not cultivars/stages.
#' @param weighting "equal" (default) or "inverse_mse". Inverse MSE requires a
#'   CV result and uses 1 / (assessment MSE + epsilon); failed members get zero
#'   weight. Fold assessment sets differ, so these scores also reflect their
#'   difficulty. Independent validation is needed to assess the weighted ensemble.
#' @param weights Optional explicit finite, non-negative weights, one per member.
#'   Named weights are matched by member name. Overrides weighting.
#' @param failure_threshold Failure fraction at or above which the ensemble date
#'   is NA. A number in (0, 1]; 0.5 means at least half the votes predict no bloom.
#' @param failure_vote "count" uses the unweighted fraction of all members;
#'   "weight" uses their base weight mass. This is not a calibrated probability.
#' @param max_weight Optional cap in (0, 1] on successful members' normalized
#'   weights. If too few positive-weight members survive to satisfy the cap,
#'   their successful weights are normalized without a cap and cap_relaxed is
#'   TRUE in the summary. Failure votes use base weights before redistribution.
#' @param epsilon Positive stabilizer in squared-day units for inverse MSE.
#' @return A pheno_ensemble containing models, normalized weights and aggregation
#'   settings. Use predict_phenology() or predict_phenology_ensemble() for new data.
#' @md
#' @export
pheno_ensemble <- function(models, weighting = c("equal", "inverse_mse"),
                            weights = NULL, failure_threshold = 0.5,
                            failure_vote = c("count", "weight"),
                            max_weight = NULL, epsilon = 1) {
  weighting <- match.arg(weighting)
  explicit_weights <- !is.null(weights)
  failure_vote <- match.arg(failure_vote)
  .ensemble_policy(failure_threshold, max_weight)
  .fit_positive_scalar(epsilon, "epsilon")
  metrics <- NULL
  if (inherits(models, "phenology_cv")) {
    metrics <- models$calibration$scores
    models <- models$model$models
  }
  if (!is.list(models) || is.data.frame(models) || !length(models))
    stop("models must be a non-empty list of model specifications or fitted results.", call. = FALSE)
  models <- lapply(models, function(m) {
    if (inherits(m, "phenology_fit")) m <- m$model
    if (!inherits(m, "pheno_model") && !inherits(m, "combined_pheno_model") &&
        !inherits(m, "pheno_model_list"))
      stop("Ensemble entries must be single, combined or model-list specifications.", call. = FALSE)
    validate_parameters(m, model_parameters(m))
    m
  })
  counts <- vapply(models, function(m) length(as_pheno_models(m)), integer(1))
  if (any(counts != counts[1]))
    stop("All ensemble members must have the same number of outputs.", call. = FALSE)
  output_names <- lapply(models, function(m) names(as_pheno_models(m)))
  if (!all(vapply(output_names, identical, logical(1), output_names[[1]])))
    stop("Ensemble cultivar/stage names must match in output order.", call. = FALSE)
  ids <- names(models)
  if (is.null(ids)) ids <- paste0("member", seq_along(models))
  .ensemble_ids(ids)
  names(models) <- ids
  if (is.null(weights)) {
    weights <- rep(1, length(models))
    if (weighting == "inverse_mse") {
      if (is.null(metrics)) stop("inverse_mse weighting requires a phenology_cv result.", call. = FALSE)
      mse <- metrics$mse[match(ids, metrics$split_id)]
      weights <- ifelse(is.finite(mse), 1 / (mse + epsilon), 0)
    }
  }
  weights <- .ensemble_weights(weights, ids)
  structure(list(models = models, weights = weights, weighting = if (explicit_weights) "explicit" else weighting,
                 failure_threshold = failure_threshold, failure_vote = failure_vote,
                 max_weight = max_weight, epsilon = epsilon), class = "pheno_ensemble")
}

#' Aggregate existing ensemble predictions with no-bloom handling
#'
#' Successful dates are averaged after renormalizing the surviving positive
#' weights. No successful weight always gives NA, even below the failure vote
#' threshold. ensemble_sd is sqrt(sum(w * (date - weighted_mean)^2)): descriptive
#' weighted spread, not a standard error or prediction interval. Zero spread
#' is returned for one surviving member. No-bloom is NA; infinite predictions
#' are rejected. The spread remains available even when the vote suppresses the
#' ensemble date, describing the successful members conditional on bloom.
#'
#' @param predicted Numeric matrix with cases in rows and ensemble members in
#'   columns; column names identify members. Alternatively, a long data frame
#'   with case_id, member_id and predicted, exactly one row per case/member pair.
#'   Missing pairs are errors; represent no bloom explicitly with NA.
#' @param weights,failure_threshold,failure_vote,max_weight As in pheno_ensemble().
#' @return A data frame in input case order with case_id, predicted, ensemble_sd,
#'   n_success (all finite predictions), n_fail, failure_fraction and cap_relaxed.
#' @md
#' @export
summarise_ensemble_predictions <- function(predicted, weights = NULL,
                                           failure_threshold = 0.5,
                                           failure_vote = c("count", "weight"),
                                           max_weight = NULL) {
  failure_vote <- match.arg(failure_vote)
  .ensemble_policy(failure_threshold, max_weight)
  if (is.data.frame(predicted)) {
    if (!all(c("case_id", "member_id", "predicted") %in% names(predicted)) || !nrow(predicted))
      stop("Long predictions require case_id, member_id and predicted columns.", call. = FALSE)
    if (anyNA(predicted$case_id) || anyNA(predicted$member_id))
      stop("Prediction IDs must be non-missing.", call. = FALSE)
    cases <- unique(predicted$case_id)
    ids <- as.character(unique(predicted$member_id))
    i <- match(predicted$case_id, cases)
    j <- match(as.character(predicted$member_id), ids)
    key <- (j - 1L) * length(cases) + i
    if (anyDuplicated(key) || length(key) != length(cases) * length(ids))
      stop("Supply exactly one prediction per case/member pair.", call. = FALSE)
    if (!is.numeric(predicted$predicted)) stop("Predicted dates must be numeric.", call. = FALSE)
    mat <- matrix(NA_real_, nrow = length(cases), ncol = length(ids))
    mat[cbind(i, j)] <- predicted$predicted
    colnames(mat) <- ids
  } else {
    mat <- predicted
    if (!is.matrix(mat) || !is.numeric(mat) || !nrow(mat) || !ncol(mat))
      stop("predicted must be a non-empty numeric matrix or long prediction data frame.", call. = FALSE)
    cases <- rownames(mat)
    if (is.null(cases)) cases <- seq_len(nrow(mat))
    ids <- colnames(mat)
    if (is.null(ids)) ids <- paste0("member", seq_len(ncol(mat)))
  }
  .ensemble_ids(ids)
  if (any(is.infinite(mat))) stop("Predicted dates must be finite or NA for no bloom.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, length(ids))
  base_weights <- .ensemble_weights(weights, ids)
  rows <- lapply(seq_len(nrow(mat)), function(i) {
    x <- mat[i, ]
    failed <- is.na(x)
    fraction <- if (failure_vote == "count") mean(failed) else sum(base_weights[failed])
    keep <- !failed & base_weights > 0
    date <- spread <- NA_real_
    relaxed <- FALSE
    if (any(keep)) {
      capped <- .ensemble_cap(base_weights[keep], max_weight)
      w <- capped$weights
      relaxed <- capped$relaxed
      date <- sum(x[keep] * w)
      spread <- sqrt(sum(w * (x[keep] - date)^2))
    }
    if (fraction >= failure_threshold) date <- NA_real_
    data.frame(predicted = date, ensemble_sd = spread, n_success = sum(!failed),
               n_fail = sum(failed), failure_fraction = fraction, cap_relaxed = relaxed)
  })
  cbind(data.frame(case_id = cases), do.call(rbind, rows))
}

.ensemble_ids <- function(ids) {
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Ensemble member names must be unique and non-empty.", call. = FALSE)
}

.ensemble_weights <- function(weights, ids) {
  if (!is.numeric(weights) || !is.null(dim(weights)) || length(weights) != length(ids) ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0))
    stop("weights must be finite and non-negative, with positive total weight and one value per member.", call. = FALSE)
  if (!is.null(names(weights))) {
    .ensemble_ids(names(weights))
    if (!setequal(names(weights), ids)) stop("Weight names must match ensemble member names.", call. = FALSE)
    weights <- weights[ids]
  }
  # Scaling first avoids overflow in the sum for large finite weights.
  weights <- weights / max(weights)
  stats::setNames(weights / sum(weights), ids)
}

.ensemble_policy <- function(threshold, cap) {
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.null(dim(x)) &&
    is.finite(x) && x > 0 && x <= 1
  if (!valid(threshold)) stop("failure_threshold must be in (0, 1].", call. = FALSE)
  if (!is.null(cap) && !valid(cap)) stop("max_weight must be NULL or in (0, 1].", call. = FALSE)
}

.ensemble_cap <- function(weights, cap) {
  weights <- weights / sum(weights)
  relaxed <- !is.null(cap) && length(weights) * cap < 1 - 1e-12
  if (is.null(cap) || relaxed) return(list(weights = weights, relaxed = relaxed))
  result <- numeric(length(weights))
  active <- seq_along(weights)
  mass <- 1
  while (length(active)) {
    proposal <- mass * weights[active] / sum(weights[active])
    exceed <- proposal > cap + 1e-12
    if (!any(exceed)) {
      result[active] <- proposal
      break
    }
    result[active[exceed]] <- cap
    mass <- mass - sum(exceed) * cap
    active <- active[!exceed]
  }
  list(weights = result, relaxed = FALSE)
}

#' Predict dates and spread for a phenology ensemble
#'
#' Each ensemble member predicts the same weather. Cultivars/stages remain
#' separate outputs; their dates are never averaged together.
#' @param model A pheno_ensemble from pheno_ensemble().
#' @param weather Weather accepted by predict_phenology(): one data frame, a
#'   seasonlist, or nested seasonlists for combined/stage outputs.
#' @param stopatzc Passed to each member's prediction function.
#' @param return_members If TRUE, return a list with summary and individual
#'   long-format predictions; otherwise return the summary data frame.
#' @return Summary columns from summarise_ensemble_predictions(), plus output
#'   (cultivar/stage label) and season. Individual predictions have case_id,
#'   member_id (ensemble fit ID) and predicted.
#'   predict_phenology(ensemble, weather) returns dates with the usual vector/list
#'   shape. With basic_output = FALSE it returns summary and individual predictions.
#' @md
#' @export
predict_phenology_ensemble <- function(model, weather, stopatzc = TRUE, return_members = FALSE) {
  .ensemble_prediction(model, weather, stopatzc, return_members)$result
}

.ensemble_prediction <- function(model, weather, stopatzc, return_members) {
  if (!inherits(model, "pheno_ensemble")) stop("Use a pheno_ensemble.", call. = FALSE)
  if (!is.logical(return_members) || length(return_members) != 1L || is.na(return_members))
    stop("return_members must be a single logical value.", call. = FALSE)
  predictions <- lapply(model$models, function(m) predict_phenology(m, weather, stopatzc = stopatzc))
  template <- predictions[[1]]
  shape <- function(p) {
    if (is.numeric(p)) return(list(lengths = length(p), labels = names(p), outputs = NULL))
    if (!is.list(p) || !all(vapply(p, is.numeric, logical(1))))
      stop("Ensemble members must predict numeric dates.", call. = FALSE)
    list(lengths = lengths(p), labels = lapply(p, names), outputs = names(p))
  }
  expected <- shape(template)
  if (!all(vapply(predictions, function(p) identical(shape(p), expected), logical(1))))
    stop("Ensemble members must predict the same cases in the same output order.", call. = FALSE)
  mat <- do.call(cbind, lapply(predictions, function(p) unname(unlist(p))))
  colnames(mat) <- names(model$models)
  summary <- summarise_ensemble_predictions(mat, model$weights, model$failure_threshold,
                                             model$failure_vote, model$max_weight)
  outputs <- if (is.numeric(template)) list(model = template) else template
  labels <- names(outputs)
  if (is.null(labels)) labels <- paste0("output", seq_along(outputs))
  summary$output <- rep(labels, lengths(outputs))
  summary$season <- unlist(lapply(outputs, function(p) {
    if (is.null(names(p))) as.character(seq_along(p)) else names(p)
  }), use.names = FALSE)
  individual <- data.frame(case_id = rep(summary$case_id, ncol(mat)),
                           member_id = rep(colnames(mat), each = nrow(mat)),
                           predicted = as.vector(mat))
  result <- if (return_members) list(summary = summary, individual_predictions = individual) else summary
  list(result = result, template = template, dates = summary$predicted)
}

.predict_pheno_ensemble <- function(model, weather, stopatzc, basic_output) {
  prediction <- .ensemble_prediction(model, weather, stopatzc, !basic_output)
  if (!basic_output) return(prediction$result)
  template <- prediction$template
  if (is.numeric(template)) return(stats::setNames(prediction$dates, names(template)))
  offset <- 0L
  lapply(template, function(p) {
    indices <- offset + seq_along(p)
    offset <<- offset + length(p)
    stats::setNames(prediction$dates[indices], names(p))
  })
}
