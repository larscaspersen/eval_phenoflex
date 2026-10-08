#' Fit models using prepared cross-validation splits
#'
#' Consumes splits made by prepare_phenology_cv(); never shuffles data or uses
#' independent validation entries. For each split, run the requested optimizer
#' restarts and retain the fit with the lowest training loss. Predict assessment
#' seasons using that fitted model and record held-out date errors.
#'
#' @param model A single pheno_model, combined_pheno_model, pheno_model_list,
#'   or plain list of simple models. Population models are unsupported because
#'   matching population predictions to observed dates needs a separate scorer.
#' @param data Prepared phenology_cv_data from prepare_phenology_cv().
#' @param evaluation Explicit training loss, accepting parameters, model,
#'   seasons, observed and any arguments in ... . Use phenology_rss(),
#'   phenology_rss_combined() or phenology_rss_stages() for the corresponding layout.
#' @param lower,upper,optimizer,control,max_iterations,max_seconds,parameters
#'   Passed to fit_phenology(). Parameters omitted from bounds remain fixed.
#' @param restarts Number of optimizer restarts per prepared split.
#' @param seed Non-negative integer or NULL; seeds individual fitting runs.
#'   Fold assignments are always taken unchanged from data.
#' @param refit If TRUE, additionally fit all CV entries, excluding independent
#'   validation. Stores its model separately in calibration$refit.
#' @param ... Additional fixed training-loss arguments. Season-dependent data
#'   must be supplied through seasons/observed in the preparation function;
#'   arguments in ... are not subset.
#' @param weighting,weights,failure_threshold,failure_vote,max_weight,epsilon
#'   Passed to pheno_ensemble() to configure the stored prediction ensemble.
#'   Equal weights are the default; independent validation is never used to set weights.
#' @param keep_diagnostics Store restart statistics, selected run settings and
#'   native optimizer results in calibration$diagnostics. Default FALSE.
#' @param keep_training_predictions Append training dates to predictions, marked
#'   with set = "training". Default FALSE retains assessment dates only.
#' @return A phenology_cv list with three elements, usable by predict_phenology():
#'   `model` is the fold-model ensemble; `predictions` pairs predicted and observed
#'   dates with case/split identifiers; `calibration` contains settings and scores,
#'   plus optional refit and diagnostics. Only reserved validation weather is
#'   retained, in calibration$validation, for validate_phenology(cv). Models are
#'   stored once in model$models. summary(cv) reports assessment scores and failures;
#'   printing previews six score rows. Independent validation remains unevaluated.
#'   MSE/RMSE are Inf if any prediction at
#'   an observed case is missing. Missing observations do not contribute to scores.
#'   These scores assess fold models, not a subsequently weighted ensemble.
#' @md
#' @export
fit_phenology_cv <- function(model, data, evaluation, lower = NULL, upper = NULL,
                             optimizer = c("DEoptim", "GenSA"), restarts = 1,
                             seed = 12345, control = list(), max_iterations = NULL,
                             max_seconds = NULL, parameters = model_parameters(model),
                             refit = FALSE, ...,
                             weighting = c("equal", "inverse_mse"), weights = NULL,
                             failure_threshold = 0.5, failure_vote = c("count", "weight"),
                             max_weight = NULL, epsilon = 1,
                             keep_diagnostics = FALSE, keep_training_predictions = FALSE) {
  .validate_cv_data(data)
  if (is.list(model) && !inherits(model, "pheno_model_list") && length(model) &&
      all(vapply(model, inherits, logical(1), what = "pheno_model")))
    model <- pheno_model_list(model)
  .cv_check_model(model, data)
  if (missing(evaluation) || !is.function(evaluation))
    stop("evaluation must be supplied as a function.", call. = FALSE)
  .fit_positive_scalar(restarts, "restarts", integer = TRUE)
  if (!is.null(seed)) .fit_positive_scalar(seed, "seed", integer = TRUE, zero = TRUE)
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit))
    stop("refit must be a single logical value.", call. = FALSE)
  for (option in c("keep_diagnostics", "keep_training_predictions")) {
    value <- get(option)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(option, " must be a single logical value.", call. = FALSE)
  }
  splits <- .cv_splits(data)
  optimizer <- match.arg(optimizer)
  weighting <- match.arg(weighting)
  failure_vote <- match.arg(failure_vote)
  .ensemble_policy(failure_threshold, max_weight)
  .fit_positive_scalar(epsilon, "epsilon")
  if (!is.null(weights)) .ensemble_weights(weights, names(splits))
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
      any(names(dots) %in% c("seasons", "observed", "prediction"))))
    stop("Additional evaluation arguments must be named; seasons/observed come from data and prediction is managed by CV.", call. = FALSE)
  n_runs <- (length(splits) + as.integer(refit)) * restarts
  seeds <- if (is.null(seed)) rep(list(NULL), n_runs) else
    as.list(.cv_with_seed(seed, function() sample.int(.Machine$integer.max, n_runs)))
  restart_metrics <- list()
  run_id <- 0L
  fit_subset <- function(indices, id) {
    subset <- .cv_subset(data, indices)
    if (data$settings$layout == "combined" && any(lengths(subset$seasons) == 0L))
      stop("Split ", id, " leaves a member without training seasons.", call. = FALSE)
    candidates <- lapply(seq_len(restarts), function(r) {
      run_id <<- run_id + 1L
      fit <- tryCatch(do.call(fit_phenology, c(list(model = model, evaluation = evaluation,
        lower = lower, upper = upper, optimizer = optimizer, seed = seeds[[run_id]],
        control = control, max_iterations = max_iterations, max_seconds = max_seconds,
        parameters = parameters, seasons = subset$seasons, observed = subset$observed,
        prediction = FALSE), dots)),
        error = function(e) stop("Split ", id, ", restart ", r, ": ", conditionMessage(e), call. = FALSE))
      if (keep_diagnostics) restart_metrics[[run_id]] <<- data.frame(split_id = id, restart = r,
        seed = if (is.null(seeds[[run_id]])) NA_integer_ else seeds[[run_id]],
        training_loss = fit$value, evaluations = fit$diagnostics$evaluations, elapsed = fit$diagnostics$elapsed)
      fit
    })
    best <- candidates[[which.min(vapply(candidates, function(f) f$value, numeric(1)))]]
    best
  }
  models <- predictions <- metrics <- vector("list", length(splits))
  selected <- list()
  names(models) <- names(splits)
  record <- function(fit, indices, split, set) {
    predicted <- .cv_predict_cases(fit$model, data, indices, skip_missing = set == "training")
    cbind(data.frame(split_id = split$id, repeat_id = split$repeat_id,
                     fold = split$fold, set = set), predicted)
  }
  save_diagnostics <- function(fit, id) {
    if (keep_diagnostics) selected[[id]] <<- list(settings = fit$calibration_settings,
                                                diagnostics = fit$diagnostics)
  }
  for (i in seq_along(splits)) {
    split <- splits[[i]]
    fit <- fit_subset(split$train, split$id)
    models[[i]] <- fit$model
    save_diagnostics(fit, split$id)
    predicted <- .cv_predict_cases(fit$model, data, split$assessment)
    predictions[[i]] <- cbind(data.frame(split_id = split$id, repeat_id = split$repeat_id,
                                       fold = split$fold, set = "assessment"), predicted)
    metrics[[i]] <- cbind(data.frame(split_id = split$id, repeat_id = split$repeat_id,
                                   fold = split$fold, training_loss = fit$value),
                          .cv_date_metrics(predicted$predicted, predicted$observed))
    if (keep_training_predictions)
      predictions[[i]] <- rbind(predictions[[i]], record(fit, split$train, split, "training"))
  }
  calibration <- list(settings = list(preparation = data$settings,
                    restarts = restarts, seed = seed, optimizer = optimizer,
                    lower = lower, upper = upper, control = control,
                    max_iterations = max_iterations, max_seconds = max_seconds,
                    parameters = parameters, evaluation = evaluation, refit = refit,
                    keep_diagnostics = keep_diagnostics,
                    keep_training_predictions = keep_training_predictions),
                    scores = do.call(rbind, metrics))
  reserved <- .cv_validation_indices(data)
  if (length(reserved)) calibration$validation <- .cv_case_records(data, reserved)
  if (refit) {
    indices <- setdiff(data$data$cases$case_id, reserved)
    full <- fit_subset(indices, "full")
    calibration$refit <- full$model
    calibration$refit_loss <- full$value
    save_diagnostics(full, "full")
    if (keep_training_predictions) predictions[[length(predictions) + 1L]] <-
      record(full, indices, list(id = "full", repeat_id = NA_integer_, fold = NA_integer_), "training")
  }
  if (keep_diagnostics) calibration$diagnostics <-
    list(restart_metrics = do.call(rbind, restart_metrics), selected = selected)
  result <- structure(list(model = list(models = models),
                           predictions = do.call(rbind, predictions), calibration = calibration),
                      class = "phenology_cv")
  result$model <- pheno_ensemble(result, weighting, weights, failure_threshold,
                                  failure_vote, max_weight, epsilon)
  rownames(result$predictions) <- NULL
  result
}

.cv_check_model <- function(model, data) {
  if (inherits(data, "phenology_cv_data")) data <- data$settings
  if (inherits(model, "phenology_cv")) model <- model$model
  if (inherits(model, "pheno_ensemble")) {
    for (m in model$models) .cv_check_model(m, data)
    return(invisible(TRUE))
  }
  if (inherits(model, "phenology_fit")) model <- model$model
  valid <- if (data$layout == "single") inherits(model, "pheno_model") else
    inherits(model, "pheno_model_list") ||
      (data$layout == "combined" && inherits(model, "combined_pheno_model"))
  if (!valid) stop("Model specification must match the prepared data layout.", call. = FALSE)
  validate_model_spec(model)
  n <- if (inherits(model, "pheno_model_list")) length(model) else
    if (inherits(model, "combined_pheno_model")) model$n_cultivars else 1L
  if (n != length(data$member_names)) stop("Model and data must have the same number of members.", call. = FALSE)
  if (inherits(model, "pheno_model_list") && data$explicit_member_names && !is.null(names(model)) &&
      !identical(names(model), data$member_names))
    stop("Named model members must match the prepared data's member order.", call. = FALSE)
  invisible(TRUE)
}

.cv_case_records <- function(data, indices) {
  cases <- .cv_cases(data)[indices, , drop = FALSE]
  first <- cases[!duplicated(cases$entry_id), , drop = FALSE]
  weather <- lapply(seq_len(nrow(first)), function(k) {
    row <- first[k, ]
    if (data$settings$layout == "combined") data$data$seasons[[row$member]][[row$season]] else
      data$data$seasons[[row$season]]
  })
  names(weather) <- as.character(first$entry_id)
  list(cases = cases, weather = weather, layout = data$settings$layout,
       member_names = data$settings$member_names,
       explicit_member_names = data$settings$explicit_member_names)
}

.cv_predict_cases <- function(model, data, indices, skip_missing = FALSE) {
  .cv_predict_records(model, .cv_case_records(data, indices), skip_missing)
}

.cv_predict_records <- function(model, data, skip_missing = FALSE) {
  if (inherits(model, "phenology_fit") || inherits(model, "phenology_cv")) model <- model$model
  cases <- data$cases
  if (inherits(model, "pheno_ensemble")) {
    individual <- lapply(model$models, function(m) .cv_predict_records(m, data, skip_missing)$predicted)
    mat <- do.call(cbind, individual)
    colnames(mat) <- names(model$models)
    cases$predicted <- summarise_ensemble_predictions(mat, model$weights,
      model$failure_threshold, model$failure_vote, model$max_weight)$predicted
    return(cases)
  }
  simple <- as_pheno_models(model)
  cases$predicted <- vapply(seq_len(nrow(cases)), function(k) {
    row <- cases[k, ]
    if (skip_missing && is.na(row$observed)) return(NA_real_)
    weather <- data$weather[[as.character(row$entry_id)]]
    predict_phenology(simple[[row$member]], weather)
  }, numeric(1))
  cases
}

.cv_date_metrics <- function(predicted, observed) {
  use <- !is.na(observed)
  n <- sum(use)
  failed <- sum(!is.finite(predicted[use]))
  mse <- if (!n) NA_real_ else if (failed) Inf else mean((predicted[use] - observed[use])^2)
  data.frame(n_observed = n, n_failed = failed, mse = mse, rmse = sqrt(mse))
}

.validate_cv_data <- function(data) {
  if (!inherits(data, "phenology_cv_data") || !is.data.frame(data$data$cases) ||
      !is.data.frame(data$assignments) || !is.list(data$settings))
    stop("data must be prepared by prepare_phenology_cv().", call. = FALSE)
  cases <- .cv_cases(data)
  n <- nrow(cases)
  assignments <- data$assignments
  invalid <- function() stop("Invalid CV split: assignments must cover each CV case once per repeat, with disjoint groups and validation entries excluded.", call. = FALSE)
  .fit_positive_scalar(data$settings$v, "v", integer = TRUE)
  .fit_positive_scalar(data$settings$repeats, "repeats", integer = TRUE)
  if (data$settings$v < 2 || !identical(cases$case_id, seq_len(n)) ||
      !all(c("case_id", "member", "season", "year", "set", "repeat_id", "fold") %in% names(assignments)) ||
      anyNA(assignments$case_id) || !all(assignments$case_id %in% cases$case_id) ||
      anyNA(assignments$set) || !all(assignments$set %in% c("cv", "validation"))) invalid()
  for (column in c("member", "season", "year"))
    if (!identical(assignments[[column]], cases[[column]][assignments$case_id])) invalid()
  validation <- .cv_validation_indices(data)
  cv_indices <- setdiff(seq_len(n), validation)
  reserved <- assignments[assignments$set == "validation", , drop = FALSE]
  rows <- assignments[assignments$set == "cv", , drop = FALSE]
  if (anyDuplicated(reserved$case_id) || any(!is.na(reserved$repeat_id)) ||
      any(!is.na(reserved$fold)) || anyNA(rows$repeat_id) || anyNA(rows$fold) ||
      !all(rows$repeat_id %in% seq_len(data$settings$repeats)) ||
      !all(rows$fold %in% seq_len(data$settings$v))) invalid()
  if (!is.character(cases$group) || length(cases$group) != n || anyNA(cases$group))
    stop("Invalid CV entry grouping.", call. = FALSE)
  if (data$settings$layout == "stages") {
    reserved_seasons <- cases$season[validation]
    if (any(cases$season[cv_indices] %in% reserved_seasons))
      stop("All stages of a validation season must be withheld together.", call. = FALSE)
  }
  for (r in seq_len(data$settings$repeats)) {
    ids <- rows$case_id[rows$repeat_id == r]
    if (anyDuplicated(ids) || !setequal(ids, cv_indices)) invalid()
  }
  for (s in .cv_splits(data)) {
    if (!any(!is.na(cases$observed[s$train])) || !any(!is.na(cases$observed[s$assessment])) ||
        any(cases$group[s$train] %in% cases$group[s$assessment]) ||
        (data$settings$layout == "stages" &&
         any(cases$entry_id[s$train] %in% cases$entry_id[s$assessment]))) invalid()
  }
  invisible(TRUE)
}

#' Evaluate a chosen model on reserved independent validation entries
#'
#' Call this after selecting the model, weights and failure policy using CV.
#' Repeatedly changing those choices using these results turns the reserved data
#' into model-selection data. fit_phenology_cv() never evaluates these entries.
#' @param model A single or CV fitted result, matching specification, or pheno_ensemble.
#' @param data Prepared data containing independent validation entries. NULL uses
#'   only the reserved validation records stored in a CV result's calibration.
#' @return A list with predictions (case records and predicted dates) and metrics
#'   (n_observed, n_failed, mse, rmse). Any failed prediction at an observed case
#'   makes MSE/RMSE infinite.
#' @md
#' @export
validate_phenology <- function(model, data = NULL) {
  if (is.null(data) && inherits(model, "phenology_cv")) {
    records <- model$calibration$validation
  } else {
    .validate_cv_data(data)
    records <- .cv_case_records(data, .cv_validation_indices(data))
  }
  if (is.null(records) || !nrow(records$cases))
    stop("No independent validation entries were reserved.", call. = FALSE)
  .cv_check_model(model, records)
  predicted <- .cv_predict_records(model, records)
  if (!any(!is.na(predicted$observed)))
    stop("Independent validation must contain at least one observed date.", call. = FALSE)
  list(predictions = predicted, metrics = .cv_date_metrics(predicted$predicted, predicted$observed))
}

#' @export
print.phenology_cv <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

#' @rdname fit_phenology_cv
#' @param object A fitted phenology_cv result.
#' @export
summary.phenology_cv <- function(object, ...) {
  assessment <- object$predictions[object$predictions$set == "assessment", , drop = FALSE]
  validation <- object$calibration$validation
  structure(list(n_models = length(object$model$models),
    n_validation = if (is.null(validation)) 0L else length(unique(validation$cases$entry_id)),
    restarts = object$calibration$settings$restarts,
    metrics = .cv_date_metrics(assessment$predicted, assessment$observed),
    scores = object$calibration$scores), class = "summary.phenology_cv")
}

#' @export
print.summary.phenology_cv <- function(x, ...) {
  cat("Fitted phenology CV:", x$n_models, "fold models;", x$restarts,
      "restart(s) per split;", x$n_validation, "independent validation entries\n")
  cat("Assessment predictions:", x$metrics$n_observed, "observed;",
      x$metrics$n_failed, "failed; pooled RMSE =", format(x$metrics$rmse), "days\n")
  print(utils::head(x$scores, 6L), row.names = FALSE, ...)
  if (nrow(x$scores) > 6L)
    cat("...", nrow(x$scores) - 6L, "more rows. Inspect $calibration$scores for all splits.\n")
  invisible(x)
}
