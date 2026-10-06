#' Specify a modular bud population model
#'
#' Wraps a [pheno_model()] specification with independent bud variation in named
#' structure parameters. Chill and potential heat submodel parameters are shared
#' by all buds. The parameter schema supplies population means; `sd` and `skew`
#' use the same parameter names. Unspecified standard deviations and shapes are
#' zero. Normal and skew-normal draws are validated without truncation or resampling.
#'
#' In partial overlap, variation is specified explicitly for `yc`, `b1`, `b2`,
#' `b3` or `ol`; there is no generic `zc` parameter. Standard deviations of heat
#' requirements must be converted along with their means when changing GDH scaling.
#' Supplied populations can represent correlated traits. Sampling is independent.
#' Cutting and constant-temperature forcing are supported by
#' [predict_population_phenology()]. The legacy [phenoflex_population()] retains
#' its own forcing index conventions.
#' @param structure,chill,heat Components passed to [pheno_model()].
#' @param n Positive integer number of buds.
#' @param sd Named numeric vector of non-negative structure-parameter standard deviations.
#' @param distribution Either `normal` or `normal_skewed`, used for varying traits.
#' @param skew Named numeric vector of finite skew-normal shape parameters.
#' @param seed Optional non-negative integer seed. Sampling with a seed preserves
#' the caller's random-number state.
#' @return A population_pheno_model specification containing `model`, `n`, `sd`,
#' `distribution`, `skew` and `seed`. [default_parameters()] and
#' [parameter_schema()] also accept this specification.
#' @examples
#' model <- population_pheno_model(
#'   structure = "phenoflex", n = 3, sd = c(yc = 1, zc = 2), seed = 17
#' )
#' parameters <- default_parameters(model)
#' buds <- sample_population_parameters(model, parameters)
#' weather <- data.frame(Temp = rep(8, 480), Year = 2008,
#'                       JDay = rep(1:20, each = 24))
#' result <- predict_population_phenology(model, weather, parameters, population = buds)
#' result$bloom_jday
#' @export
population_pheno_model <- function(structure = "phenoflex", chill = chill_dynamic(),
                                   heat = heat_gdh(), n = 100, sd = numeric(),
                                   distribution = c("normal", "normal_skewed"),
                                   skew = numeric(), seed = 12345) {
  model <- pheno_model(structure, chill, heat)
  .population_scalar(n, "n", positive = TRUE, integer = TRUE)
  if (!is.null(seed)) .population_scalar(seed, "seed", nonnegative = TRUE, integer = TRUE)
  distribution <- match.arg(distribution)
  schema <- parameter_schema(model)
  traits <- schema$name[schema$component == "structure"]
  sd <- .population_trait_settings(sd, traits, "sd", nonnegative = TRUE)
  skew <- .population_trait_settings(skew, traits, "skew")
  structure(list(model = model, n = as.integer(n), sd = sd,
                 distribution = distribution, skew = skew, seed = seed),
            class = "population_pheno_model")
}

.population_trait_settings <- function(x, traits, name, nonnegative = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x)) ||
      (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
                    anyDuplicated(names(x)) || !all(names(x) %in% traits))) ||
      (nonnegative && any(x < 0)))
    stop(name, " must be a finite named vector of structure parameters",
         if (nonnegative) " with non-negative values." else ".", call. = FALSE)
  result <- stats::setNames(rep(0, length(traits)), traits)
  result[names(x)] <- x
  result
}

#' Sample named structure parameters for a bud population
#'
#' Returns all structure parameters for every bud, including constant traits.
#' Explicit samples can be reused across objective evaluations to avoid changes
#' in random draws during fitting. With `seed = NULL`, sampling advances the RNG.
#' @param model Specification from [population_pheno_model()].
#' @param parameters Named population means and shared submodel parameters.
#' @return A numeric buds-by-structure-parameters matrix. Column names are the
#' structure parameter names; row order identifies buds.
#' @export
sample_population_parameters <- function(model, parameters = default_parameters(model)) {
  if (!inherits(model, "population_pheno_model"))
    stop("Use population_pheno_model().", call. = FALSE)
  validate_parameters(model, parameters)
  if (!is.null(model$seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(model$seed)
  }
  traits <- names(model$sd)
  population <- matrix(0, nrow = model$n, ncol = length(traits),
                       dimnames = list(NULL, traits))
  for (trait in traits) {
    population[, trait] <- .draw_population_requirement(model$n, parameters[[trait]],
      model$sd[[trait]], model$distribution, model$skew[[trait]])
  }
  .validate_modular_population(model, parameters, population)
  population
}

.validate_modular_population <- function(model, parameters, population) {
  traits <- names(model$sd)
  if (!is.matrix(population) || !is.numeric(population) ||
      nrow(population) != model$n || ncol(population) != length(traits) ||
      is.null(colnames(population)) || anyNA(colnames(population)) ||
      anyDuplicated(colnames(population)) || !setequal(colnames(population), traits) ||
      any(!is.finite(population)))
    stop("population must be a finite numeric matrix with n rows and exactly the structure parameter columns.",
         call. = FALSE)
  for (i in seq_len(model$n)) {
    bud <- parameters
    bud[colnames(population)] <- population[i, ]
    tryCatch(validate_parameters(model$model, bud), error = function(e)
      stop("Invalid population bud ", i, ": ", conditionMessage(e), call. = FALSE))
  }
  invisible(TRUE)
}

#' Predict phenology for a modular bud population
#'
#' Calculates shared chill and potential heat once, then applies the selected
#' C++ structure for each bud. [predict_phenology()] dispatches here when passed
#' a population specification. Each bud retains the first bloom index; zero
#' denotes no bloom and its calendar result is NA. All buds remain in the output.
#' Detailed output contains the full unrestricted shared chill trajectory and
#' hours-by-buds effective heat. Early stopping leaves zero-filled heat tails.
#' @param model Specification from [population_pheno_model()].
#' @param weather Complete hourly weather as required by [predict_phenology()].
#' @param parameters Named population means and shared submodel parameters.
#' @param population Optional numeric matrix from [sample_population_parameters()].
#' Supplying it bypasses sampling and permits correlated or externally sampled traits.
#' @param stopatzc Stop each bud at its first bloom if TRUE.
#' @param basic_output Omit chill and effective heat trajectories if TRUE.
#' @return A list with `bloomindex`, `bloom_jday` and `population`. Detailed output
#' additionally contains `chill` and `z`; the columns of `z` correspond to population rows.
#' With cuts, `forcing` contains `cut_indices`, cuts-by-buds matrices
#' `hours_to_bloom`, `z_at_cut` and `z_final`. Elapsed forcing hours start at zero
#' at cutting, unlike the legacy kernel's zero-based step indices: one means one
#' complete forcing interval, zero means already bloomed, and NA means not reached.
#' Detailed forcing output adds a list `z`, one (hours+1)-by-buds matrix per cut.
#' @param cut_indices Optional distinct one-based weather row indices at which
#' buds are cut. A cut uses the state stored at that row, before its next interval.
#' @param forcing_temperature Constant forcing temperature in degrees Celsius.
#' @param max_hours_forcing Positive integer number of forcing intervals.
#' During forcing all chill pools are held at their cutting values. Effective
#' heat and structure history continue from the field; partial overlap retains
#' its original chill-completion baseline and overlap state.
#' @export
predict_population_phenology <- function(model, weather,
                                         parameters = default_parameters(model),
                                         population = NULL, stopatzc = TRUE,
                                         basic_output = TRUE, cut_indices = NULL,
                                         forcing_temperature = 23,
                                         max_hours_forcing = 1200) {
  if (!inherits(model, "population_pheno_model"))
    stop("Use population_pheno_model().", call. = FALSE)
  for (flag in list(stopatzc = stopatzc, basic_output = basic_output))
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("Output flags must be single non-missing logical values.", call. = FALSE)
  validate_parameters(model, parameters)
  .validate_sequential_weather(weather)
  if (!is.null(cut_indices)) {
    if (!is.numeric(cut_indices) || !is.null(dim(cut_indices)) ||
        !length(cut_indices) || any(!is.finite(cut_indices)) ||
        any(cut_indices != floor(cut_indices)) || any(cut_indices < 1) ||
        any(cut_indices > nrow(weather)) || anyDuplicated(cut_indices))
      stop("cut_indices must be distinct valid one-based weather row indices.", call. = FALSE)
    .population_scalar(forcing_temperature, "forcing_temperature")
    if (forcing_temperature <= -273)
      stop("forcing_temperature must exceed -273 degrees C.", call. = FALSE)
    .population_scalar(max_hours_forcing, "max_hours_forcing", positive = TRUE, integer = TRUE)
  }
  if (is.null(population)) population <- sample_population_parameters(model, parameters)
  else .validate_modular_population(model, parameters, population)
  population <- population[, names(model$sd), drop = FALSE]
  inputs <- .calculate_pheno_inputs(model$model, weather, parameters)
  bloomindex <- numeric(model$n)
  z <- if (!basic_output) matrix(0, nrow(weather), model$n) else NULL
  for (i in seq_len(model$n)) {
    bud <- parameters
    bud[colnames(population)] <- population[i, ]
    result <- .apply_pheno_structure(model$model, inputs$chill, inputs$heat,
                                     bud, stopatzc, basic_output)
    bloomindex[i] <- result$bloomindex
    if (!basic_output) z[, i] <- result$z
  }
  out <- list(bloomindex = bloomindex,
    bloom_jday = vapply(bloomindex, return_JDay, numeric(1),
                       Jday_vec = weather$JDay, year_vec = weather$Year),
    population = population)
  if (!basic_output) {
    out$chill <- inputs$chill
    out$z <- z
  }
  if (!is.null(cut_indices))
    out$forcing <- .predict_population_forcing(model, parameters, population,
      inputs, cut_indices, forcing_temperature, max_hours_forcing, basic_output)
  out
}

.predict_population_forcing <- function(model, parameters, population, inputs,
                                        cuts, temperature, hours, basic_output) {
  # Replaying the field prefix preserves each structure's internal history.
  # Only potential heat is newly calculated; chill is held at the cutting state.
  heat_fn <- if (model$model$heat$scaling == "scaled")
    calculate_heat_gdh else calculate_heat_gdh_unscaled
  forcing_heat <- heat_fn(rep(temperature, hours + 1), seq_len(hours + 1),
    Tb = parameters[["Tb"]], Tu = parameters[["Tu"]], Tc = parameters[["Tc"]])
  elapsed <- z_cut <- z_final <- matrix(NA_real_, length(cuts), model$n)
  trajectories <- if (!basic_output) vector("list", length(cuts)) else NULL
  for (j in seq_along(cuts)) {
    cut <- cuts[j]
    chill <- rbind(inputs$chill[seq_len(cut), , drop = FALSE],
                   inputs$chill[rep(cut, hours), , drop = FALSE])
    heat <- c(inputs$heat[seq_len(cut - 1)], forcing_heat)
    if (!basic_output) trajectories[[j]] <- matrix(0, hours + 1, model$n)
    for (i in seq_len(model$n)) {
      bud <- parameters
      bud[colnames(population)] <- population[i, ]
      result <- .apply_pheno_structure(model$model, chill, heat, bud,
                                       stopatzc = FALSE, basic_output = FALSE)
      if (result$bloomindex > 0)
        elapsed[j, i] <- max(0, result$bloomindex - cut)
      z_cut[j, i] <- result$z[cut]
      z_final[j, i] <- tail(result$z, 1)
      if (!basic_output) trajectories[[j]][, i] <- result$z[cut:(cut + hours)]
    }
  }
  out <- list(cut_indices = cuts, hours_to_bloom = elapsed,
              z_at_cut = z_cut, z_final = z_final)
  if (!basic_output) out$z <- trajectories
  out
}
