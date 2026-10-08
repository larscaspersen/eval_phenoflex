#' Calibrate a phenology model
#'
#' Minimize a user-supplied evaluation function using DEoptim or GenSA. The
#' evaluation function is called as `evaluation(parameters, model, ...)`, or
#' `evaluation(parameters, model, observed = observed, ...)` when observations
#' are supplied, with
#' a complete named parameter vector in schema order. It must return one numeric
#' loss to minimize; `Inf` excludes a candidate. Numerical characteristic-to-kinetic
#' conversion failures raised by predict_phenology() are excluded with `Inf`;
#' other errors in the evaluation function propagate. Candidates outside the
#' model's parameter domains are excluded before evaluation.
#'
#' @param model A [pheno_model()], [population_pheno_model()] or
#'   [combined_pheno_model()] specification, or a [pheno_model_list()]. A plain
#'   list of single models is wrapped with pheno_model_list()'s default sharing.
#' @param evaluation Function accepting parameters as its first argument and
#'   model as its second argument, followed by any arguments supplied in `...`.
#'   When observed is supplied, also accept a named `observed` argument or `...`.
#'   Required: explicitly pass the evaluation function appropriate to your
#'   calibration. [phenology_rss()] is an example for one observed date per season.
#'   Use [phenology_rss_combined()] for nested cultivar seasonlists and observations.
#'   Use [phenology_rss_stages()] for several stages on one common seasonlist.
#' @param observed Observed phenology dates, passed unchanged as the named
#'   `observed` argument to the supplied evaluation function on every evaluation.
#'   Supply a vector for single-model fitting or nested lists for cultivar/stage
#'   fitting, as required by the evaluator; NULL slots can represent missing
#'   observations for the example RSS evaluators. Matching predictions to dates
#'   and calculating the loss remain the evaluator's responsibility. This
#'   argument must be named. If omitted it is not passed, allowing evaluators
#'   that capture observations in a closure or use another interface.
#' @param prediction Optional function called as prediction(fitted_model),
#'   returning dates to store with the observations. NULL automatically predicts
#'   calibration dates for the three supplied RSS evaluators when seasons are
#'   supplied. Other evaluators store predicted = NULL unless a callback is given.
#'   FALSE disables date prediction. A callback can capture custom calibration
#'   data; observations and callback results are retained without reshaping.
#' @param lower,upper Finite named numeric vectors for the parameters to calibrate.
#'   NULL uses [default_bounds()]. If both are NULL, all parameters receive
#'   suggested bounds. If only one is NULL, defaults are used only for the names
#'   in the supplied vector. Review the suggested bounds for your species and stage.
#'   Named theta_star and theta_c bounds in 0--20 are interpreted as Celsius,
#'   including indexed collection names. Each endpoint is converted independently
#'   using +273; optimizers, evaluators and returned bounds use Kelvin.
#'   Both must name the same subset of model parameters; order is arbitrary.
#'   Omitted parameters remain fixed at their initial values in `parameters`
#'   (the model's stored values by default). Equal bounds also fix a parameter.
#'   Supply `numeric()` for both bounds to keep all parameters fixed.
#' @param optimizer Either `"DEoptim"` (default) or `"GenSA"`.
#' @param seed Optional non-negative integer seed. A supplied seed initializes
#'   both R's RNG and GenSA's internal RNG and preserves the caller's RNG state.
#'   It overrides `control$seed`. Time-limited runs need not be reproducible.
#' @param control Named list of optimizer-specific controls. See
#'   [DEoptim::DEoptim.control()] or [GenSA::GenSA()]. DEoptim supports
#'   `parallelType = "parallel"`, `"foreach"`, `"auto"` and `cluster`.
#'   For `"parallel"`, install parallelly to determine the available worker count;
#'   supply a PSOCK cluster to choose the count explicitly. The fitter loads
#'   evalpheno and `control$packages` on these workers using the caller's library
#'   paths. It stops only clusters it creates, including on errors. Export custom
#'   global dependencies with `control$parVar` or `parallel::clusterExport()`.
#'   For `"foreach"`, register and manage a backend yourself; its workers must
#'   be able to find evalpheno and any `foreachArgs$.packages`. Export custom
#'   dependencies with `foreachArgs$.export`. Deterministic DEoptim evaluators
#'   are reproducible with a supplied seed; stochastic evaluators can depend on
#'   the backend and worker count. GenSA defaults to `smooth = FALSE`.
#' @param max_iterations Optional positive integer iteration limit; overrides
#'   `control$itermax` for DEoptim or `control$maxit` for GenSA. Iterations are
#'   generations for DEoptim and annealing steps for GenSA, not evaluation counts.
#'   If NULL, the optimizer's control/default limit applies.
#' @param max_seconds Optional positive optimizer running-time limit in seconds,
#'   passed directly to GenSA as `control$max.time`. Overrides a value in `control`.
#'   DEoptim has no native time-limit control, so this argument is rejected for
#'   DEoptim; use `max_iterations` or its native stopping controls. The fitter
#'   measures elapsed time for reporting only. The baseline evaluation precedes
#'   optimization and is outside GenSA's time budget. An evaluation already
#'   running can finish after the deadline. The optimizer's stopping controls
#'   apply together.
#' @param parameters Initial complete named parameter vector, defaulting to the
#'   model's stored parameters. Also supplies the constant values for parameters
#'   omitted from the bounds. Bounded parameters must lie within the bounds. Used as GenSA's
#'   starting point and evaluated as a baseline for both optimizers. DEoptim
#'   initializes its population using its native controls.
#'   Initial theta_star and theta_c values also accept the automatic Celsius
#'   convention; evaluators and fitted parameters receive Kelvin values.
#' @md
#' @param ... Additional arguments passed to `evaluation`, such as `seasons`.
#'   Use a closure to adapt evaluation functions with other interfaces.
#' @return A `phenology_fit` list with `par` (best complete named parameters),
#'   `value` (loss), and `model` (specification storing fitted parameters).
#'   Pass the result directly to predict_phenology(). `predictions` contains
#'   `predicted` calibration dates and the original `observed` dates (NULL when
#'   unavailable). Automatic RSS predictions use NA for unobserved cases, keeping
#'   the original season/member positions without evaluating their weather.
#'   `calibration_settings` contains optimizer, complete lower/upper bounds,
#'   seed, effective control, requested limits, initial parameters and evaluation.
#'   `diagnostics` contains evaluations (including invalid candidates/baseline),
#'   elapsed (seconds including baseline), termination ("optimizer" or "fixed"),
#'   and optimizer_result (native result or NULL). An error is raised if no
#'   finite loss is found. Termination does not assert global optimality.
#'   For model lists, `model[[i]]` is an ordinary single model with its fitted
#'   parameters. [model_parameters()] returns the flat vector and [as_pheno_models()]
#'   materializes a plain list of single models from a fitted result.
#' @examples
#' model <- pheno_model()
#' lower <- c(yc = 10)
#' upper <- c(yc = 80)
#' loss <- function(parameters, model, target) (parameters[["yc"]] - target)^2
#' fit <- fit_phenology(model, loss, lower, upper, target = 45,
#'                      seed = 17, max_iterations = 10,
#'                      control = list(NP = 10, trace = FALSE))
#' fit$par
#' # Explicit observations with the example RSS evaluator:
#' model <- pheno_model("sequential", chill_dynamic("kinetic"),
#'                      parameters = c(yc = 0.1, zc = 5))
#' seasons <- list(data.frame(Temp = rep(8, 480), Year = 2008,
#'                            JDay = rep(50:69, each = 24)))
#' observed <- predict_phenology(model, seasons)
#' fit <- fit_phenology(model, evaluation = phenology_rss, observed = observed,
#'                      seasons = seasons, lower = numeric(), upper = numeric())
#' fit$value
#' predict_phenology(fit, seasons)
#' fit$predictions
#' @export
fit_phenology <- function(model, evaluation, lower = NULL, upper = NULL,
                          optimizer = c("DEoptim", "GenSA"), seed = 12345,
                          control = list(), max_iterations = NULL,
                          max_seconds = NULL,
                          parameters = model_parameters(model), ..., observed = NULL,
                          prediction = NULL) {
  if (is.list(model) && !inherits(model, "pheno_model_list") && length(model) &&
      all(vapply(model, inherits, logical(1), what = "pheno_model")))
    model <- pheno_model_list(model)
  if (!inherits(model, "pheno_model") && !inherits(model, "population_pheno_model") &&
      !inherits(model, "combined_pheno_model") && !inherits(model, "pheno_model_list"))
    stop("model must be a pheno_model, population_pheno_model, combined_pheno_model or pheno_model_list.",
         call. = FALSE)
  parameters <- .normalize_characteristic_temperatures(parameters)
  validate_parameters(model, parameters)
  if (missing(evaluation) || !is.function(evaluation))
    stop("evaluation must be supplied as a function.", call. = FALSE)
  has_observed <- !missing(observed)
  if (!is.null(prediction) && !identical(prediction, FALSE) && !is.function(prediction))
    stop("prediction must be NULL, FALSE or a function of the fitted model.", call. = FALSE)
  optimizer <- match.arg(optimizer)
  parameter_names <- parameter_schema(model)$name
  parameters <- parameters[parameter_names]
  if (is.null(lower) || is.null(upper)) {
    defaults <- default_bounds(model)
    if (is.null(lower) && is.null(upper)) {
      lower <- defaults$lower
      upper <- defaults$upper
    } else if (is.null(lower)) {
      upper <- .fit_bounds(upper, parameter_names, "upper")
      lower <- defaults$lower[names(upper)]
    } else {
      lower <- .fit_bounds(lower, parameter_names, "lower")
      upper <- defaults$upper[names(lower)]
    }
  }
  supplied_lower <- .fit_bounds(lower, parameter_names, "lower")
  supplied_upper <- .fit_bounds(upper, parameter_names, "upper")
  if (!setequal(names(supplied_lower), names(supplied_upper)))
    stop("lower and upper must name the same subset of model parameters.", call. = FALSE)
  lower <- upper <- parameters
  lower[names(supplied_lower)] <- supplied_lower
  upper[names(supplied_upper)] <- supplied_upper
  if (any(lower > upper)) stop("lower must not exceed upper.", call. = FALSE)
  if (any(parameters < lower | parameters > upper))
    stop("Initial parameters must lie within lower and upper.", call. = FALSE)
  .fit_control_list(control)
  if (!is.null(seed)) .fit_positive_scalar(seed, "seed", integer = TRUE, zero = TRUE)
  if (!is.null(max_iterations))
    .fit_positive_scalar(max_iterations, "max_iterations", integer = TRUE)
  if (!is.null(max_seconds)) .fit_positive_scalar(max_seconds, "max_seconds")

  if (optimizer == "DEoptim") {
    if (!is.null(max_seconds) || !is.null(control$max.time))
      stop("DEoptim has no native time limit; use max_iterations or its stopping controls.",
           call. = FALSE)
    if (!is.null(max_iterations)) control$itermax <- max_iterations
    if (!is.null(control$itermax))
      .fit_positive_scalar(control$itermax, "control$itermax", integer = TRUE)
  } else {
    if (is.null(control$smooth)) control$smooth <- FALSE
    if (!is.null(max_iterations)) control$maxit <- max_iterations
    if (!is.null(max_seconds)) control$max.time <- max_seconds
    if (!is.null(control$maxit))
      .fit_positive_scalar(control$maxit, "control$maxit", integer = TRUE)
    if (!is.null(control$max.time))
      .fit_positive_scalar(control$max.time, "control$max.time")
    if (!is.null(seed)) control$seed <- -as.integer(max(1, seed))
  }

  free <- lower < upper
  if (any(free) && !requireNamespace(optimizer, quietly = TRUE))
    stop("Install the '", optimizer, "' package to use this optimizer.", call. = FALSE)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  started <- proc.time()[["elapsed"]]
  elapsed <- function() proc.time()[["elapsed"]] - started
  best_par <- parameters
  best_value <- Inf
  evaluations <- 0L
  evaluate <- .fit_objective(parameters, free, lower, upper, model, evaluation,
                            if (has_observed) c(list(observed = observed), list(...)) else list(...))
  # GenSA uses local tracking; DEoptim's native result also works across workers.
  objective <- if (optimizer == "DEoptim") evaluate else function(x) {
    evaluations <<- evaluations + 1L
    value <- evaluate(x)
    if (value < best_value) {
      best_value <<- unname(value)
      best_par[free] <<- x
    }
    value
  }
  best_value <- evaluate(parameters[free])
  evaluations <- 1L
  termination <- if (any(free)) "optimizer" else "fixed"
  optimizer_result <- NULL
  if (any(free)) {
    if (optimizer == "DEoptim") {
      de_control <- do.call(DEoptim::DEoptim.control, control)
      if (!is.null(de_control$cluster) || de_control$parallelType == "parallel") {
        cl <- de_control$cluster
        if (is.null(cl)) {
          if (!requireNamespace("parallelly", quietly = TRUE))
            stop("Install 'parallelly' or supply control$cluster for parallel DEoptim.", call. = FALSE)
          cl <- parallel::makePSOCKcluster(parallelly::availableCores())
          on.exit(parallel::stopCluster(cl), add = TRUE)
        } else if (!inherits(cl, "cluster")) {
          stop("control$cluster must be a parallel cluster object.", call. = FALSE)
        }
        initialize_worker <- function(paths, packages) {
          .libPaths(paths)
          for (package in packages) library(package, character.only = TRUE)
          NULL
        }
        # Bootstrap library paths before a worker deserializes our namespace.
        environment(initialize_worker) <- baseenv()
        parallel::clusterCall(cl, initialize_worker, .libPaths(),
                              unique(c("evalpheno", de_control$packages)))
        if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
        de_control$cluster <- cl
        # A supplied cluster takes precedence. Keep lifecycle ownership here.
        de_control$parallelType <- "none"
      } else if (de_control$parallelType == "foreach") {
        de_control$foreachArgs$.packages <- unique(c("evalpheno", de_control$packages,
                                                    de_control$foreachArgs$.packages))
      }
      optimizer_result <- DEoptim::DEoptim(objective, lower[free], upper[free], control = de_control)
      evaluations <- 1L + optimizer_result$optim$nfeval
      if (optimizer_result$optim$bestval < best_value) {
        best_value <- unname(optimizer_result$optim$bestval)
        best_par[free] <- optimizer_result$optim$bestmem
      }
    } else {
      optimizer_result <- GenSA::GenSA(parameters[free], objective,
                                       lower[free], upper[free], control = control)
    }
  }
  if (!is.finite(best_value))
    stop("No candidate produced a finite evaluation loss.", call. = FALSE)
  fitted_model <- model
  if (inherits(model, "pheno_model_list")) {
    fitted_model <- as_pheno_models(model, best_par)
    attributes(fitted_model) <- attributes(model)
  } else if (inherits(model, "population_pheno_model")) fitted_model$model$parameters <- best_par
  else fitted_model$parameters <- best_par
  optimization_elapsed <- elapsed()
  dates <- if (is.function(prediction)) prediction(fitted_model) else
    if (is.null(prediction) && has_observed)
      .fit_calibration_dates(fitted_model, evaluation, list(...)$seasons, observed) else NULL
  structure(list(par = best_par, value = best_value, model = fitted_model,
                 predictions = list(predicted = dates, observed = if (has_observed) observed else NULL),
                 calibration_settings = list(optimizer = optimizer, lower = lower, upper = upper,
                   seed = seed, control = control, max_iterations = max_iterations,
                   max_seconds = max_seconds, parameters = parameters, evaluation = evaluation),
                 diagnostics = list(evaluations = evaluations, elapsed = optimization_elapsed,
                   termination = termination, optimizer_result = optimizer_result)), class = "phenology_fit")
}

# Keep the worker closure separate from the fitter's tracking and cluster state.
.fit_objective <- function(parameters, free, lower, upper, model, evaluation, arguments) {
  force(parameters); force(free); force(lower); force(upper)
  force(model); force(evaluation); force(arguments)
  function(x) {
    candidate <- parameters
    candidate[free] <- x
    if (any(!is.finite(candidate)) || any(candidate < lower | candidate > upper))
      return(Inf)
    valid <- tryCatch({validate_parameters(model, candidate); TRUE},
                      error = function(e) FALSE)
    if (!valid) return(Inf)
    value <- tryCatch(do.call(evaluation, c(list(candidate, model), arguments)),
                      phenology_conversion_error = function(e) Inf)
    if (!is.numeric(value) || length(value) != 1L || !is.null(dim(value)) ||
        is.na(value) || value == -Inf)
      stop("evaluation must return one numeric loss (finite or Inf).", call. = FALSE)
    unname(value)
  }
}

.fit_calibration_dates <- function(model, evaluation, seasons, observed) {
  if (is.null(seasons)) return(NULL)
  predict_member <- function(m, weather, dates) {
    normalized <- .phenology_observations(dates, length(weather))
    predicted <- rep(NA_real_, length(weather))
    names(predicted) <- names(weather)
    use <- which(!is.na(normalized))
    if (length(use)) predicted[use] <- predict_phenology(m, weather[use])
    predicted
  }
  if (identical(evaluation, phenology_rss))
    return(predict_member(model, seasons, observed))
  if (identical(evaluation, phenology_rss_combined) || identical(evaluation, phenology_rss_stages)) {
    simple <- as_pheno_models(model)
    weather <- if (identical(evaluation, phenology_rss_stages)) rep(list(seasons), length(simple)) else seasons
    result <- lapply(seq_along(simple), function(i) predict_member(simple[[i]], weather[[i]], observed[[i]]))
    names(result) <- if (!is.null(names(observed))) names(observed) else names(simple)
    return(result)
  }
  NULL
}

.fit_bounds <- function(x, parameter_names, name) {
  if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x)) ||
      (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
                    anyDuplicated(names(x)) || !all(names(x) %in% parameter_names))))
    stop(name, " must be a finite named numeric vector containing only model parameter names.",
         call. = FALSE)
  .normalize_characteristic_temperatures(x[parameter_names[parameter_names %in% names(x)]])
}

.fit_control_list <- function(control) {
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(names(control) == "") || anyDuplicated(names(control)))))
    stop("control must be a named list without duplicate names.", call. = FALSE)
}

.fit_positive_scalar <- function(x, name, integer = FALSE, zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.null(dim(x)) || !is.finite(x) ||
      x < 0 || (!zero && x == 0) ||
      (integer && (x != floor(x) || x > .Machine$integer.max)))
    stop(name, " must be a ", if (zero) "non-negative" else "positive",
         if (integer) " integer." else " finite number.", call. = FALSE)
}
