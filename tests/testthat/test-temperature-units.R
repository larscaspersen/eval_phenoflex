temperature_unit_weather <- function() {
  data.frame(Temp = rep(8, 480), Year = 2008, JDay = rep(50:69, each = 24))
}

test_that("model constructors normalize Celsius independently and preserve Kelvin", {
  reference <- pheno_model()
  integer_parameters <- pheno_model("sequential")$parameters
  integer_parameters[] <- round(integer_parameters)
  integer_parameters <- setNames(as.integer(integer_parameters), names(integer_parameters))
  expect_identical(pheno_model("sequential", parameters = integer_parameters)$parameters,
                   integer_parameters)
  for (values in list(c(theta_star = 6, theta_c = 13.1),
                      c(theta_star = 279, theta_c = 13.1),
                      c(theta_star = 6, theta_c = 286.1))) {
    model <- pheno_model(parameters = values)
    expect_identical(model, reference)
  }
  supplied <- rev(reference$parameters)
  supplied[c("theta_star", "theta_c")] <- c(6, 13.1)
  expect_identical(pheno_model(parameters = supplied)$parameters, rev(reference$parameters))
  expect_identical(supplied[c("theta_star", "theta_c")], c(theta_star = 6, theta_c = 13.1))
  expect_identical(pheno_model(parameters = c(theta_star = 0))$parameters[["theta_star"]], 273)
  expect_identical(pheno_model(parameters = c(theta_star = 19, theta_c = 20))$parameters[c("theta_star", "theta_c")],
                   c(theta_star = 292, theta_c = 293))
  expect_identical(pheno_model(parameters = c(theta_star = 21))$parameters[["theta_star"]], 21)
  expect_error(pheno_model(parameters = c(theta_star = 14, theta_c = 13)), "theta_star < theta_c")
  expect_error(pheno_model(parameters = c(theta_star = -1)), "theta_star")
  expect_error(pheno_model(parameters = c(theta_c = 300)), "theta_c")
  expect_error(pheno_model(parameters = c(theta_star = NA_real_)), "finite")
  # Other parameters use their existing units even when their values are small.
  model <- pheno_model(parameters = c(Tf = 6, Tb = 6, Tu = 20, Tc = 30, yc = 15))
  expect_identical(model$parameters[c("Tf", "Tb", "Tu", "Tc", "yc")],
                   c(Tf = 6, Tb = 6, Tu = 20, Tc = 30, yc = 15))
})

test_that("validation and predictions accept Celsius overrides without mutation", {
  weather <- temperature_unit_weather()
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    model <- pheno_model(structure)
    p <- model$parameters
    p[c("theta_star", "theta_c")] <- c(6, 13.1)
    before <- p
    expect_true(validate_parameters(model, p))
    for (basic in c(TRUE, FALSE)) {
      expect_identical(predict_phenology(model, weather, p, basic_output = basic),
                       predict_phenology(model, weather, basic_output = basic))
      expect_identical(predict_phenology(model, list(a = weather), p, basic_output = basic),
                       predict_phenology(model, list(a = weather), basic_output = basic))
    }
    expect_identical(as_pheno_models(model, p)[[1]]$parameters, model$parameters)
    expect_identical(p, before)
  }
  population <- population_pheno_model(n = 2)
  p <- model_parameters(population)
  p[c("theta_star", "theta_c")] <- c(6, 13.1)
  expect_identical(predict_phenology(population, weather, p), predict_phenology(population, weather))
})

test_that("collections normalize shared and indexed characteristic temperatures", {
  base <- pheno_model()
  combined <- combined_pheno_model(base, 2)
  p <- model_parameters(combined)
  p[c("theta_star", "theta_c")] <- c(6, 13.1)
  expect_identical(combined_pheno_model(base, 2, parameters = p), combined)
  expect_identical(cultivar_parameters(combined, 1, p), base$parameters)
  expect_identical(as_pheno_models(combined, p), as_pheno_models(combined))
  expect_identical(predict_phenology(combined, temperature_unit_weather(), p),
                   predict_phenology(combined, temperature_unit_weather()))
  shared <- setdiff(names(base$parameters), c("theta_star", "theta_c"))
  models <- pheno_model_list(list(a = base, b = base), shared = shared, ordered = "theta_star")
  p <- model_parameters(models)
  p[c("theta_star1", "theta_c1", "theta_star2", "theta_c2")] <- c(6, 286.1, 279, 13.1)
  expect_true(validate_parameters(models, p))
  expect_identical(as_pheno_models(models, p), as_pheno_models(models))
  input <- list(base, base)
  input[[2]]$parameters[c("theta_star", "theta_c")] <- c(6, 13.1)
  expect_identical(pheno_model_list(input), pheno_model_list(list(base, base)))
  expect_identical(input[[2]]$parameters[["theta_star"]], 6)
})

test_that("default bounds optionally display Celsius without changing other units", {
  for (model in list(pheno_model(), pheno_model(chill = chill_dynamic("kinetic")),
                     combined_pheno_model(pheno_model(), 2),
                     pheno_model_list(list(pheno_model(), pheno_model()), shared = character()))) {
    kelvin <- default_bounds(model)
    celsius <- default_bounds(model, temperature_unit = "C")
    expected <- lapply(kelvin, function(p) {
      temps <- grepl("^(theta_star|theta_c)[0-9]*$", names(p))
      p[temps] <- p[temps] - 273
      p
    })
    expect_identical(celsius, expected)
    expect_identical(default_bounds(model, temperature_unit = "K"), kelvin)
  }
  bounds <- default_bounds(pheno_model(), c("theta_star", "theta_c"), temperature_unit = "C")
  expect_identical(bounds, list(lower = c(theta_star = 6, theta_c = 13),
                                upper = c(theta_star = 8, theta_c = 14)))
  expect_error(default_bounds(pheno_model(), temperature_unit = "F"), "arg")
})

test_that("both optimizer interfaces receive Kelvin for Celsius bounds and starts", {
  model <- pheno_model()
  initial <- model$parameters
  initial[c("theta_star", "theta_c")] <- c(6, 13.1)
  bounds <- default_bounds(model, c("theta_star", "theta_c"), temperature_unit = "C")
  expected <- default_bounds(model, c("theta_star", "theta_c"))
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    received <- NULL
    testthat::local_mocked_bindings(DEoptim = function(fn, lower, upper, control) {
      received <<- list(lower = lower, upper = upper)
      candidate <- (lower + upper) / 2
      value <- fn(candidate)
      list(optim = list(bestmem = candidate, bestval = value, nfeval = 1L))
    }, .package = "DEoptim")
    testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
      received <<- list(lower = lower, upper = upper)
      fn((lower + upper) / 2)
      list()
    }, .package = "GenSA")
    fit <- fit_phenology(model, function(p, m) {
      expect_gt(p[["theta_star"]], 273)
      expect_gt(p[["theta_c"]], 273)
      (p[["theta_star"]] - 280)^2
    }, bounds$lower, bounds$upper, optimizer = optimizer, parameters = initial,
    max_iterations = 1)
    expect_identical(received, expected)
    expect_identical(fit$par[["theta_star"]], 280)
    expect_identical(fit$model$parameters, fit$par)
    expect_identical(fit$calibration_settings$lower[names(expected$lower)], expected$lower)
  }
  fit <- fit_phenology(model, function(p, m) sum(p), lower = c(theta_star = 6), upper = c(theta_star = 279))
  expect_identical(fit$diagnostics$termination, "fixed")
  expect_identical(fit$par, model$parameters)
})

test_that("partial default bounds and indexed fit bounds normalize before comparison", {
  skip_if_not_installed("GenSA")
  received <- NULL
  testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
    received <<- list(lower = lower, upper = upper)
    list(value = fn(par))
  }, .package = "GenSA")
  model <- pheno_model()
  loss <- function(p, m) sum(p)
  fit_phenology(model, loss, upper = c(theta_c = 14), optimizer = "GenSA")
  expect_identical(received, list(lower = c(theta_c = 286), upper = c(theta_c = 287)))
  fit_phenology(model, loss, lower = c(theta_c = 13), optimizer = "GenSA")
  expect_identical(received, list(lower = c(theta_c = 286), upper = c(theta_c = 287)))
  models <- pheno_model_list(list(model, model), shared = character())
  fit <- fit_phenology(models, loss, lower = c(theta_c1 = 13, theta_c2 = 286),
                       upper = c(theta_c1 = 287, theta_c2 = 14), optimizer = "GenSA")
  expect_identical(received, list(lower = c(theta_c1 = 286, theta_c2 = 286),
                                 upper = c(theta_c1 = 287, theta_c2 = 287)))
  expect_identical(fit$par, model_parameters(models))
  expect_error(fit_phenology(model, loss, lower = c(theta_c = 14), upper = c(theta_c = 13)),
                "lower must not exceed")
})
