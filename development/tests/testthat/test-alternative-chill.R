alternative_chill_weather <- function() {
  data.frame(Temp = rep(c(5, 8, 20, 5), each = 24), Year = 2008,
    JDay = rep(1:4, each = 24), Hour = rep(0:23, 4))
}

test_that("chill requirements use submodel-specific defaults throughout model construction", {
  choices <- list(chill_dynamic(), chill_hours(), chill_utah(), chill_positive_utah())
  expected <- c(dynamic = 40, chill_hours = 1000, utah = 1200, positive_utah = 1200)
  for (chill in choices) {
    yc <- unname(expected[chill$name])
    for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
      model <- pheno_model(structure, chill = chill)
      expect_identical(model$parameters[["yc"]], yc)
      expect_identical(default_parameters(model)[["yc"]], yc)
      expect_identical(pheno_model(structure, chill = chill, parameters = c(yc = 500))$parameters[["yc"]], 500)
    }
    model <- pheno_model(chill = chill)
    population <- population_pheno_model(chill = chill, n = 2)
    expect_identical(model_parameters(population)[["yc"]], yc)
    combined <- combined_pheno_model(model, 2)
    expect_identical(unname(model_parameters(combined)[c("yc1", "yc2")]), rep(yc, 2))
    stages <- stage_pheno_models(model, c(180, 250))
    expect_identical(model_parameters(stages)[["yc"]], yc)
    expect_true(validate_parameters(model, model_parameters(model)))
  }
  expect_false("yc" %in% names(pheno_model(chill = chill_none())$parameters))
})

test_that("alternative chill requirement bounds include their default calibration starts", {
  skip_if_not_installed("GenSA")
  testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
    fn(par)
    list()
  }, .package = "GenSA")
  for (chill in list(chill_hours(), chill_utah(), chill_positive_utah())) {
    model <- pheno_model(chill = chill)
    bounds <- default_bounds(model, "yc")
    expected <- if (chill$name == "chill_hours") c(250, 2000) else c(300, 2400)
    expect_equal(c(unname(bounds$lower), unname(bounds$upper)), expected)
    fit <- fit_phenology(model, function(p, m) (p[["yc"]] - m$parameters[["yc"]])^2,
      optimizer = "GenSA", max_iterations = 1)
    expect_identical(fit$value, 0)
    expect_identical(fit$par[["yc"]], model$parameters[["yc"]])
    combined <- combined_pheno_model(model, 2)
    expect_equal(unname(default_bounds(combined, "yc")$lower), rep(expected[1], 2))
    expect_equal(unname(default_bounds(combined, "yc")$upper), rep(expected[2], 2))
  }
})

test_that("chill-hours bounds are inclusive and Celsius/Kelvin inputs agree", {
  temp <- c(-0.1, 0, 7.2, 7.21, 5)
  times <- seq_along(temp)
  expected <- matrix(c(0, 0, 1, 2, 2), ncol = 1, dimnames = list(NULL, "y"))
  expect_identical(calculate_chill_hours(temp, times), expected)
  expect_equal(calculate_chill_hours(temp + 273, times, deg_celsius = FALSE), expected)
  expect_identical(calculate_chill_hours(temp, times, CH_lower = -1, CH_upper = 8)[, "y"], 0:4 * 1.0)
  expect_error(calculate_chill_hours(temp, times, CH_lower = 8, CH_upper = 7), "less than")
  expect_error(calculate_chill_hours(temp, times, CH_lower = -273), "greater than")
})

test_that("Utah zone boundaries retain signed losses and recover before coupling", {
  temp <- c(1.4, 2.4, 9.1, 12.4, 15.9, 18, 19, 5)
  expected <- c(0, 0, 0.5, 1.5, 2, 2, 1.5, 0.5)
  chill <- calculate_chill_utah(temp, seq_along(temp))
  expect_identical(chill[, "y_raw"], expected)
  expect_identical(chill[, "y"], expected)
  expect_equal(calculate_chill_utah(temp + 273, seq_along(temp), deg_celsius = FALSE), chill)
  warm <- c(20, 20, 5, 5, 5, 5)
  chill <- calculate_chill_utah(warm, seq_along(warm))
  expect_identical(chill[, "y_raw"], c(0, -1, -2, -1, 0, 1))
  expect_identical(chill[, "y"], c(0, 0, 0, 0, 0, 1))
  expect_identical(calculate_chill_utah(warm, seq_along(warm), weights = rep(2, 7))[, "y_raw"],
                   c(0, 2, 4, 6, 8, 10))
})

test_that("positive Utah clips daily sums rather than individual hours", {
  temp <- c(rep(5, 12), rep(20, 12), rep(5, 24), 5)
  chill <- calculate_chill_positive_utah(temp, seq_along(temp))[, "y"]
  expect_identical(chill[1:48], rep(0, 48))
  expect_identical(chill[49], 24)
  expect_equal(calculate_chill_positive_utah(temp + 273, seq_along(temp), deg_celsius = FALSE)[, "y"], chill)
  custom <- calculate_chill_positive_utah(rep(5, 6), 1:6, day_group = c(1, 1, 2, 2, 3, 3))
  expect_identical(custom[, "y"], c(0, 0, 2, 2, 4, 5))
  expect_identical(unname(calculate_chill_positive_utah(rep(5, 48), 1:48)[48, "y"]), 47)
})

test_that("heat-only timing uses start-of-interval activation and retains first bloom", {
  temp <- rep(25, 8)
  heat <- calculate_heat_gdh_unscaled(temp, seq_along(temp))
  for (start in c(1, 3, 3.5, 8, 20)) {
    chill <- calculate_chill_none(temp, seq_along(temp), heat_start_index = start)
    result <- apply_heat_only_structure(chill, heat, zc = 2, stopatzc = FALSE, basic_output = FALSE)
    expected <- c(0, cumsum(as.numeric(1:7 >= start)))
    expect_identical(result$z, expected)
    crossing <- which(expected >= 2)
    expect_equal(result$bloomindex, if (length(crossing)) crossing[1] else 0)
    early <- apply_heat_only_structure(chill, heat, zc = 2)
    expect_identical(early$bloomindex, result$bloomindex)
  }
})

test_that("alternative chill kernels reject malformed input and overflow", {
  calculators <- list(calculate_chill_hours, calculate_chill_utah,
    calculate_chill_positive_utah, calculate_chill_none)
  for (fn in calculators) {
    expect_error(fn(c(5, NA, 5), 1:3), "finite")
    expect_error(fn(rep(5, 3), c(1, 2, 4)), "consecutive hourly")
    expect_error(fn(5, 1), "at least two")
    expect_error(fn(c(-273, 5), 1:2), "greater than")
  }
  for (fn in list(calculate_chill_utah, calculate_chill_positive_utah)) {
    expect_error(fn(rep(5, 3), 1:3, weights = 1:6), "seven")
    expect_error(fn(rep(5, 3), 1:3, weights = c(1:6, NA)), "finite")
    expect_error(fn(rep(5, 3), 1:3, weights = matrix(rep(1, 7))), "vector")
    expect_error(fn(rep(5, 3), 1:3, weights = rep(.Machine$double.xmax, 7)), "overflowed")
  }
  for (groups in list(c(1, 2), c(1, 2, 1), c(1, NA, 2), c(1, 1.5, 2), matrix(1:3)))
    expect_error(calculate_chill_positive_utah(rep(5, 3), 1:3, day_group = groups), "day_group")
  expect_error(calculate_chill_none(rep(5, 3), 1:3, heat_start_index = 0), "at least one")
  expect_error(apply_heat_only_structure(matrix(c(0, 0.5, 1), ncol = 1), c(1, 1)), "zero or one")
})

test_that("modular predictions match the supplied C++ reference engine", {
  fixture <- readRDS(test_path("fixtures", "chill-reference", "reference.rds"))
  for (id in names(fixture$records)) {
    record <- fixture$records[[id]]
    chill <- switch(record$chill, dynamic = chill_dynamic("kinetic"),
      chill_hours = chill_hours(), utah = chill_utah(), positive_utah = chill_positive_utah(), none = chill_none())
    model <- pheno_model(chill = chill, heat = heat_gdh(record$scaling), parameters = record$parameters)
    result <- predict_phenology(model, fixture$weather, stopatzc = record$stopatzc, basic_output = FALSE)
    raw <- if (record$chill == "utah") "y_raw" else "y"
    expect_equal(result$chill[, raw], record$y, tolerance = 1e-12, info = id)
    expect_equal(result$z, record$z, tolerance = 1e-12, info = id)
    expect_equal(result$bloomindex, record$first_bloomindex, info = id)
  }
})

test_that("schemas, bounds and predictions support every coupling and scaling", {
  weather <- alternative_chill_weather()
  for (chill in list(chill_hours(), chill_utah(), chill_positive_utah(), chill_none())) {
    for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
      predictions <- list()
      for (scaling in c("unscaled", "scaled")) {
        model <- pheno_model(structure, chill, heat_gdh(scaling))
        p <- model_parameters(model)
        if (chill$name == "none") {
          expect_named(p, c("zc", "heat_start_index", "Tb", "Tu", "Tc"))
          p[c("zc", "heat_start_index")] <- c(1, 25)
        } else {
          p["yc"] <- 2
          requirements <- intersect(c("zc", "b1", "b2"), names(p))
          p[requirements] <- 1
        }
        if (scaling == "scaled") p[intersect(c("zc", "b1", "b2"), names(p))] <-
          p[intersect(c("zc", "b1", "b2"), names(p))] * 21
        expect_true(validate_parameters(model, rev(p)))
        bounds <- default_bounds(model)
        expect_named(bounds$lower, parameter_schema(model)$name)
        expect_true(all(is.finite(unlist(bounds))))
        before <- model
        result <- predict_phenology(model, weather, p, basic_output = FALSE, stopatzc = FALSE)
        predictions[[scaling]] <- predict_phenology(model, weather, p)
        expect_equal(predictions[[scaling]], return_JDay(result$bloomindex, weather$JDay, weather$Year))
        expect_true(all(result$z >= 0))
        expect_identical(model, before)
      }
      expect_equal(predictions$unscaled, predictions$scaled)
    }
  }
  expect_error(pheno_model(chill = chill_hours(), parameters = c(CH_lower = 8, CH_upper = 7)), "CH_lower")
  expect_error(pheno_model(chill = chill_none(), parameters = c(yc = 1)), "schema names")
  expect_error(pheno_model(chill = chill_none(), parameters = c(heat_start_index = 0)), "at least one")
  expect_error(pheno_model(chill = list(name = "unknown")), "Use chill_dynamic")
  expect_error(pheno_model(chill = list(name = "utah", weights = rep(1, 7))), "only name")
  delayed <- pheno_model(chill = chill_none(), parameters = c(heat_start_index = 25))
  for (specification in list(delayed, stage_pheno_models(delayed, c(180, 250)),
                            combined_pheno_model(delayed, 2))) {
    expect_identical(default_bounds(specification)$lower[["heat_start_index"]], 25)
    expect_identical(default_bounds(specification)$upper[["heat_start_index"]], 25)
  }
})

test_that("new chill submodels work in fitting, stages, combined models and populations", {
  weather <- alternative_chill_weather()
  seasons <- list(a = weather, b = weather)
  for (chill in list(chill_hours(), chill_utah(), chill_positive_utah(), chill_none())) {
    parameters <- if (chill$name == "none") c(zc = 1) else c(yc = 2, zc = 1)
    model <- pheno_model(chill = chill, parameters = parameters)
    observations <- predict_phenology(model, seasons)
    fit <- fit_phenology(model, phenology_rss, numeric(), numeric(), seasons = seasons, observed = observations)
    expect_equal(fit$value, 0)
    stages <- stage_pheno_models(model, c(early = 1, late = 2))
    expect_identical(predict_phenology(stages, seasons), predict_phenology(as_pheno_models(stages), seasons))
    combined <- combined_pheno_model(model, cultivar_names = c("Gala", "Fuji"))
    expect_identical(predict_phenology(combined, seasons), predict_phenology(as_pheno_models(combined), seasons))
    population <- population_pheno_model(chill = chill, n = 2)
    buds <- if (chill$name == "none") matrix(c(1, 2), ncol = 1, dimnames = list(NULL, "zc")) else
      cbind(yc = c(2, 2), zc = c(1, 2), s1 = c(0.5, 0.5))
    out <- predict_population_phenology(population, weather, parameters = model_parameters(model),
      population = buds, basic_output = FALSE, cut_indices = c(1, 25), max_hours_forcing = 24)
    expect_equal(out$bloom_jday, unlist(predict_phenology(stages, weather), use.names = FALSE))
    expect_equal(dim(out$forcing$hours_to_bloom), c(2, 2))
    ensemble <- pheno_ensemble(list(model, model))
    expect_equal(predict_phenology(ensemble, seasons), observations)
    data <- prepare_phenology_cv(seasons, observations, v = 2)
    cv <- fit_phenology_cv(model, data, phenology_rss, numeric(), numeric())
    expect_equal(cv$calibration$scores$mse, c(0, 0))
  }
})

test_that("heat-only activation continues through forcing after an early cutting", {
  weather <- alternative_chill_weather()
  model <- population_pheno_model(chill = chill_none(), n = 2)
  p <- model_parameters(model)
  p[c("zc", "heat_start_index")] <- c(1, 25)
  buds <- matrix(c(1, 2), ncol = 1, dimnames = list(NULL, "zc"))
  out <- predict_population_phenology(model, weather, parameters = p, population = buds,
    cut_indices = c(1, 25), forcing_temperature = 25, max_hours_forcing = 48)
  expect_identical(out$forcing$hours_to_bloom[, 1], c(25, 1))
  expect_identical(out$forcing$hours_to_bloom[, 2], c(26, 2))
})

test_that("both optimizers can fit a chill-hours threshold through predictions", {
  weather <- alternative_chill_weather()
  weather$Temp <- rep(c(8, 5, 20, 8), each = 24)
  model <- pheno_model(chill = chill_hours(), parameters = c(yc = 2, zc = 1))
  p <- model_parameters(model)
  p["CH_upper"] <- 10
  seasons <- list(season = weather)
  observed <- predict_phenology(model, seasons, p)
  expect_true(is.finite(observed))
  expect_gt(phenology_rss(model_parameters(model), model, seasons, observed), 0)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 20, trace = FALSE) else list(max.call = 100)
    fit <- fit_phenology(model, phenology_rss,
      lower = c(CH_lower = -1, CH_upper = 7.2), upper = c(CH_lower = 6, CH_upper = 11),
      seasons = seasons, observed = observed, optimizer = optimizer, max_iterations = 5,
      seed = 17, control = control)
    expect_equal(fit$value, 0)
    expect_gte(fit$par[["CH_upper"]], 8)
    expect_equal(predict_phenology(fit, seasons), observed)
  }
})
