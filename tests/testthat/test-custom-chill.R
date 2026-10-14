custom_chill_weather <- function(days = 4) {
  data.frame(Temp = rep(c(5, 8, 5, 20), length.out = days * 24), Year = 2008,
    JDay = rep(seq_len(days), each = 24), Hour = rep(0:23, days))
}

custom_window_chill <- function(lower = NULL, upper = NULL) {
  chill_custom(function(weather, parameters) {
    evalpheno::calculate_chill_hours(weather$Temp, seq_len(nrow(weather)),
      CH_lower = parameters[["tlo"]], CH_upper = parameters[["thi"]])
  }, defaults = c(tlo = 0, thi = 7.2), lower = lower, upper = upper,
    validate = function(p) p[["tlo"]] < p[["thi"]], label = "Custom window")
}

test_that("custom constructors and callbacks meet the hourly matrix contract", {
  custom <- custom_window_chill()
  weather <- custom_chill_weather()
  expect_true(validate_chill_custom(custom, weather))
  expect_true(validate_chill_custom(custom, weather, rev(custom$defaults)))
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    for (scaling in c("unscaled", "scaled")) {
      custom_model <- pheno_model(structure, custom, heat_gdh(scaling))
      builtin <- pheno_model(structure, chill_hours(), heat_gdh(scaling))
      p <- model_parameters(custom_model)
      p["yc"] <- 2
      p[intersect(c("zc", "b1", "b2"), names(p))] <- if (scaling == "scaled") 21 else 1
      builtin_p <- model_parameters(builtin)
      builtin_p[intersect(names(builtin_p), names(p))] <- p[intersect(names(builtin_p), names(p))]
      for (basic in c(TRUE, FALSE)) for (stop in c(TRUE, FALSE))
        expect_equal(predict_phenology(custom_model, weather, rev(p), basic_output = basic, stopatzc = stop),
          predict_phenology(builtin, weather, builtin_p, basic_output = basic, stopatzc = stop))
      expect_output(print(custom_model), "Chill submodel: Custom window")
    }
  }
  empty <- chill_custom(function(weather, parameters) {
    stopifnot(length(parameters) == 0)
    matrix(c(0L, rep(1L, nrow(weather) - 1)), ncol = 1)
  })
  model <- pheno_model(chill = empty)
  expect_false(any(parameter_schema(model)$component == "chill"))
  detail <- predict_phenology(model, weather, basic_output = FALSE, stopatzc = FALSE)
  expect_type(detail$chill, "double")
  expect_named(as.data.frame(detail$chill), "y")
})

test_that("custom metadata rejects malformed parameters, bounds and signatures", {
  callback <- function(weather, parameters) matrix(0, nrow(weather), 1)
  for (defaults in list(1, c(a = 1, a = 2), stats::setNames(1, ""), c(a = NA_real_),
                        matrix(1), stats::setNames(1, NA_character_)))
    expect_error(chill_custom(callback, defaults), "defaults")
  for (name in c("yc", "zc", "s1", "b1", "Tb", "theta_star", "theta_c2", "yc1", "heat_start_index"))
    expect_error(chill_custom(callback, stats::setNames(1, name)), "reserved")
  expect_error(chill_custom(1), "function")
  expect_error(chill_custom(function(x, y) matrix(0, 72, 1)), "weather and parameters")
  expect_error(chill_custom(function(weather, parameters, extra) extra), "callback failed")
  expect_error(chill_custom(callback, label = " "), "label")
  expect_error(chill_custom(callback, validate = 1), "validate")
  expect_error(chill_custom(callback, c(a = 1), lower = c(a = 0)), "both lower and upper")
  expect_error(chill_custom(callback, c(a = 1), lower = c(a = 0), upper = c(b = 2)), "bounds")
  expect_error(chill_custom(callback, c(a = 1), lower = c(a = 2), upper = c(a = 0)), "must not exceed")
  expect_error(chill_custom(callback, c(a = 1), lower = c(a = 2), upper = c(a = 3)), "lie within")
  expect_error(chill_custom(callback, lower = c(a = 0), upper = c(a = 1)), "declared default")
  spec <- chill_custom(callback, defaults = c(a = 1))
  spec$defaults <- c(yc = 1)
  expect_error(pheno_model(chill = spec), "reserved")
})

test_that("custom trajectory and determinism checks reject invalid callbacks", {
  callbacks <- list(
    function(weather, parameters) rep(0, nrow(weather)),
    function(weather, parameters) data.frame(y = rep(0, nrow(weather))),
    function(weather, parameters) matrix(0, nrow(weather) - 1, 1),
    function(weather, parameters) matrix(numeric(), nrow(weather), 0),
    function(weather, parameters) matrix(FALSE, nrow(weather), 1))
  for (callback in callbacks) expect_error(chill_custom(callback), "numeric matrix")
  for (value in c(NA_real_, Inf, NaN)) {
    callback <- local({ v <- value; function(weather, parameters) matrix(v, nrow(weather), 1) })
    expect_error(chill_custom(callback), "finite")
  }
  expect_error(chill_custom(function(weather, parameters) matrix(-1, nrow(weather), 1)), "non-negative")
  expect_error(chill_custom(function(weather, parameters) {
    matrix(0, nrow(weather), 2, dimnames = list(NULL, c("y", "y")))
  }), "column names")
  expect_error(chill_custom(function(weather, parameters) {
    matrix(stats::runif(nrow(weather)), ncol = 1)
  }), "deterministic")
  signed <- chill_custom(function(weather, parameters) {
    cbind(raw = -seq_len(nrow(weather)), y = rep(0, nrow(weather)))
  })
  expect_true(validate_chill_custom(signed))
})

test_that("extra weather requirements and runtime failures are checked explicitly", {
  callback <- function(weather, parameters) {
    if (!"Light" %in% names(weather)) stop("Light column is required")
    matrix(c(0, cumsum(head(weather$Light, -1))), ncol = 1)
  }
  expect_error(chill_custom(callback), "Light column")
  weather <- custom_chill_weather()
  weather$Light <- 1
  spec <- chill_custom(callback, check_weather = weather)
  expect_false("check_weather" %in% names(spec))
  model <- pheno_model(chill = spec, parameters = c(yc = 2, zc = 1))
  expect_true(is.finite(predict_phenology(model, weather)))
  weather$Light[4] <- NA
  expect_error(predict_phenology(model, weather), "finite")
  spec <- chill_custom(function(weather, parameters) matrix(0, 72, 1))
  model <- pheno_model(chill = spec)
  expect_error(predict_phenology(model, custom_chill_weather()), "one row per weather")
  expect_error(fit_phenology(model, phenology_rss, numeric(), numeric(),
    seasons = list(custom_chill_weather()), observed = 2), "one row per weather")
})

test_that("custom validators reject invalid domains before candidate evaluation", {
  callback <- function(weather, parameters) matrix(0, nrow(weather), 1)
  expect_error(chill_custom(callback, c(width = 0), validate = function(p) p[["width"]] > 0), "rejected")
  expect_error(chill_custom(callback, c(width = 0), validate = function(p) stop("width must be positive")),
               "width must be positive")
  expect_error(chill_custom(callback, c(width = 1), validate = function(p) NA), "TRUE or NULL")
  spec <- chill_custom(callback, c(width = 1), validate = function(p) stopifnot(p[["width"]] > 0))
  model <- pheno_model(chill = spec)
  p <- model_parameters(model)
  p["width"] <- -1
  expect_error(validate_parameters(model, p), "parameter validation failed")
  free <- names(p) == "width"
  lower <- upper <- model_parameters(model)
  lower["width"] <- -2; upper["width"] <- 2
  called <- FALSE
  objective <- .fit_objective(model_parameters(model), free, lower, upper, model,
    function(parameters, model) { called <<- TRUE; 0 }, list())
  expect_identical(objective(-1), Inf)
  expect_false(called)
  expect_identical(objective(1), 0)
  expect_true(called)
})

test_that("custom bounds expand correctly and unbounded parameters stay fixed", {
  spec <- custom_window_chill(lower = c(thi = 5), upper = c(thi = 15))
  model <- pheno_model(chill = spec, parameters = c(tlo = 2))
  for (m in list(model, stage_pheno_models(model, c(180, 250)), combined_pheno_model(model, 2))) {
    bounds <- default_bounds(m)
    expect_identical(bounds$lower[["tlo"]], 2)
    expect_identical(bounds$upper[["tlo"]], 2)
    expect_identical(bounds$lower[["thi"]], 5)
    expect_identical(bounds$upper[["thi"]], 15)
  }
  individual <- pheno_model_list(list(model, model), shared = character())
  bounds <- default_bounds(individual)
  expect_identical(unname(bounds$lower[c("tlo1", "tlo2")]), c(2, 2))
  expect_identical(unname(bounds$upper[c("thi1", "thi2")]), c(15, 15))
  p <- model_parameters(model)
  p["thi"] <- 20
  expect_true(validate_parameters(model, p)) # Search bounds are not a domain.
  unchanged <- model
  predict_phenology(model, custom_chill_weather(), p)
  expect_identical(model, unchanged)
})

test_that("custom callbacks integrate with fitting and collection/population workflows", {
  weather <- custom_chill_weather()
  other <- weather; other$Year <- 2009
  seasons <- list(first = weather, second = other)
  custom <- custom_window_chill()
  model <- pheno_model(chill = custom, parameters = c(yc = 2, zc = 1))
  observed <- predict_phenology(model, seasons)
  fit <- fit_phenology(model, phenology_rss, numeric(), numeric(), seasons = seasons, observed = observed)
  expect_equal(fit$value, 0)
  stages <- stage_pheno_models(model, c(early = 1, late = 2))
  expect_equal(predict_phenology(stages, seasons), predict_phenology(as_pheno_models(stages), seasons))
  combined <- combined_pheno_model(model, cultivar_names = c("Gala", "Fuji"))
  expect_equal(predict_phenology(combined, seasons), predict_phenology(as_pheno_models(combined), seasons))
  population <- population_pheno_model(chill = custom, n = 2)
  buds <- cbind(yc = c(2, 2), zc = c(1, 2), s1 = c(0.5, 0.5))
  result <- predict_population_phenology(population, weather, parameters = model_parameters(model),
    population = buds, cut_indices = c(1, 25), max_hours_forcing = 24)
  expect_equal(result$bloom_jday, unlist(predict_phenology(stages, weather), use.names = FALSE))
  expect_equal(dim(result$forcing$hours_to_bloom), c(2, 2))
  data <- prepare_phenology_cv(seasons, observed, v = 2)
  cv <- fit_phenology_cv(model, data, phenology_rss, numeric(), numeric())
  expect_equal(cv$calibration$scores$mse, c(0, 0))
  expect_equal(predict_phenology(pheno_ensemble(list(model, model)), seasons), observed)
})

test_that("custom chill parameters can be fitted with either optimizer", {
  weather <- custom_chill_weather()
  weather$Temp <- rep(c(8, 5, 20, 8), each = 24)
  model <- pheno_model(chill = custom_window_chill(), parameters = c(yc = 2, zc = 1))
  target <- model_parameters(model); target["thi"] <- 10
  seasons <- list(weather)
  observed <- predict_phenology(model, seasons, target)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 20, trace = FALSE) else list(max.call = 100)
    fit <- fit_phenology(model, phenology_rss, lower = c(tlo = -1, thi = 7.2), upper = c(tlo = 6, thi = 11),
      seasons = seasons, observed = observed, optimizer = optimizer, max_iterations = 5, seed = 17, control = control)
    expect_equal(fit$value, 0)
    expect_equal(predict_phenology(fit, seasons), observed)
  }
})

test_that("custom callbacks preserve captured constants when serialized", {
  callback <- function(weather, parameters) {
    increments <- rep(parameters[["rate"]] + offset, nrow(weather) - 1L)
    matrix(c(0, cumsum(increments)), ncol = 1)
  }
  captured <- new.env(parent = baseenv())
  captured$offset <- 0.2
  environment(callback) <- captured
  custom <- chill_custom(callback, defaults = c(rate = 1),
    validate = function(p) p[["rate"]] > 0)
  model <- pheno_model(chill = custom, parameters = c(yc = 2, zc = 1))
  weather <- custom_chill_weather()
  weather$Temp <- 25
  restored <- unserialize(serialize(model, NULL))
  expect_true(validate_chill_custom(restored$chill, weather))
  expect_equal(predict_phenology(restored, weather), predict_phenology(model, weather))
  expect_equal(get("offset", environment(restored$chill$fun)), 0.2)
})
