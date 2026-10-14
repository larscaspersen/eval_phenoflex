test_that("construction infers Dynamic Model input and stores kinetic coefficients", {
  original <- c(E0 = 4153.5, E1 = 12888.8, A0 = 139500, A1 = 2.567e18)
  for (structure in c("phenoflex", "sequential", "parallel", "partial_overlap")) {
    model <- pheno_model(structure)
    expect_identical(model$parameters[names(original)], original)
    expect_identical(pheno_model(structure, chill_dynamic("characteristic")), model)
    expect_null(model$chill$characteristic)
    expect_identical(pheno_model(structure, parameters = c(E0 = 4000))$parameters[["E0"]], 4000)
    characteristic <- c(theta_star = 279, theta_c = 286.1, tau = 47.7, pie_c = 28)
    supplied <- c(theta_star = 6, theta_c = 13.1, yc = 45)
    converted <- pheno_model(structure, parameters = supplied)
    expect_identical(converted$chill$characteristic, characteristic)
    expected <- characteristic_to_kinetic(c(rep(0, 4), unname(characteristic), rep(0, 4)))[5:8]
    expect_equal(unname(converted$parameters[names(original)]), expected)
    expect_identical(converted$parameters[["yc"]], 45)
    expect_identical(model_parameters(converted, "characteristic")[names(characteristic)], characteristic)
    expect_true(validate_parameters(converted, converted$parameters))
    expect_false(any(names(characteristic) %in% names(converted$parameters)))
    expect_identical(supplied, c(theta_star = 6, theta_c = 13.1, yc = 45))
  }
  for (ch in c("theta_star", "theta_c", "tau", "pie_c"))
    for (ki in c("E0", "E1", "A0", "A1"))
      expect_error(pheno_model(parameters = setNames(c(1, 2), c(ch, ki))), "Cannot mix")
  expect_error(pheno_model(parameters = c(theta_star = 14, theta_c = 13)), "theta_star < theta_c")
  testthat::local_mocked_bindings(.fit_chill_parameters = function(...) list(termcd = 3L),
    .package = "evalpheno")
  expect_error(pheno_model(parameters = c(theta_star = 6)), class = "phenology_conversion_error")
})

test_that("characteristic calibration produces kinetic fitted models across layouts", {
  base <- pheno_model(parameters = c(theta_star = 6, theta_c = 13.1))
  models <- list(base, population_pheno_model(n = 2), combined_pheno_model(base, 2),
    stage_pheno_models(base, c(180, 250)), pheno_model_list(list(base, base)),
    pheno_model_list(list(base, base), shared = character()))
  weather <- data.frame(Temp = rep(8, 480), Year = 2008, JDay = rep(50:69, each = 24))
  for (model in models) {
    initial <- model_parameters(model, "characteristic")
    fit <- fit_phenology(model, function(p, m) {
      expect_true(any(grepl("^theta_star", names(p))))
      validate_parameters(m, p)
      sum(p)
    }, lower = numeric(), upper = numeric(), parameters = initial)
    expect_identical(fit$par, initial)
    expect_identical(fit$calibration_settings$parameterization, "characteristic")
    expect_true(any(grepl("^E0", names(model_parameters(fit)))))
    expect_false(any(grepl("^theta_star", names(model_parameters(fit)))))
    expect_true(validate_parameters(fit$model, model_parameters(fit)))
    if (!inherits(model, "population_pheno_model"))
      for (m in as_pheno_models(fit)) {
        expect_true(all(c("E0", "E1", "A0", "A1") %in% names(m$parameters)))
        expect_true(all(c("theta_star", "theta_c", "tau", "pie_c") %in% names(m$chill$characteristic)))
      }
    expect_equal(predict_phenology(fit, weather), predict_phenology(model, weather), tolerance = 1e-6)
  }
})

test_that("calibration bounds infer representation and mixed inputs fail early", {
  model <- pheno_model(parameters = c(theta_star = 6))
  fit <- fit_phenology(model, function(p, m) p[["theta_star"]],
    lower = c(theta_star = 6), upper = c(theta_star = 6))
  expect_identical(fit$par[["theta_star"]], 279)
  expect_identical(fit$model$chill$characteristic[["theta_star"]], 279)
  expect_identical(names(default_bounds(model)$lower), names(model$parameters))
  expect_identical(names(default_bounds(model, parameterization = "characteristic")$lower),
    names(model_parameters(model, "characteristic")))
  expect_identical(default_bounds(model, "theta_star", temperature_unit = "C")$lower,
    c(theta_star = 6))
  expect_error(fit_phenology(model, sum, lower = c(theta_star = 6), upper = c(E0 = 5000)), "Cannot mix")
  expect_error(fit_phenology(model, sum, lower = c(theta_star = 6), upper = c(theta_star = 8),
    parameters = model$parameters), "Cannot mix")
  expect_error(fit_phenology(model, sum, lower = c(E0 = 2500), upper = c(E0 = 5000),
    parameterization = "characteristic"), "conflicts")
  expect_error(default_bounds(model, c("theta_star", "E0")), "Cannot mix")
  expect_error(model_parameters(pheno_model_list(list(model, model), shared = "E0"), "characteristic"),
    "sharing all four")
})

test_that("different characteristic values remain valid collection members", {
  first <- pheno_model(parameters = c(theta_star = 6))
  second <- pheno_model(parameters = c(theta_star = 7))
  collection <- pheno_model_list(list(first, second), shared = character())
  expect_true(validate_model_spec(collection))
  p <- model_parameters(collection, "characteristic")
  expect_identical(p[c("theta_star1", "theta_star2")], c(theta_star1 = 279, theta_star2 = 280))
  expect_equal(as_pheno_models(collection, p)[[2]]$parameters, second$parameters)
  changed <- first
  changed$parameters["E0"] <- changed$parameters[["E0"]] + 1
  expect_false(identical(model_parameters(changed, "characteristic")["theta_star"],
    model_parameters(first, "characteristic")["theta_star"]))
})

test_that("changed characteristic candidates update kinetic coefficients and metadata", {
  skip_if_not_installed("GenSA")
  base <- pheno_model(parameters = c(theta_star = 6))
  models <- pheno_model_list(list(base, base), shared = character())
  testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
    fn(c(theta_star1 = 280, theta_star2 = 279.5))
    list()
  }, .package = "GenSA")
  fit <- fit_phenology(models, function(p, m) sum((p[c("theta_star1", "theta_star2")] - c(280, 279.5))^2),
    lower = c(theta_star1 = 6, theta_star2 = 6), upper = c(theta_star1 = 8, theta_star2 = 8),
    optimizer = "GenSA")
  expect_identical(fit$value, 0)
  for (i in 1:2) {
    expected <- pheno_model(parameters = c(theta_star = c(280, 279.5)[i]))
    expect_identical(fit$model[[i]]$parameters, expected$parameters)
    expect_identical(fit$model[[i]]$chill$characteristic, expected$chill$characteristic)
  }
  expect_false(any(grepl("theta", attr(fit$model, "shared"))))
})

test_that("characteristic CV stores kinetic fold models and refits", {
  x <- cv_example()
  x$model <- pheno_model("sequential", parameters = c(theta_star = 6, yc = 0.1, zc = 5))
  x$observed <- predict_phenology(x$model, x$seasons)
  data <- prepare_phenology_cv(x$seasons, x$observed, v = 2, seed = 17)
  cv <- fit_phenology_cv(x$model, data, phenology_rss, lower = numeric(), upper = numeric(),
    parameterization = "characteristic", refit = TRUE)
  expect_identical(cv$calibration$settings$parameterization, "characteristic")
  for (m in c(cv$model$models, list(cv$calibration$refit))) {
    expect_identical(m$parameters, x$model$parameters)
    expect_identical(m$chill$characteristic, x$model$chill$characteristic)
  }
  expect_equal(predict_phenology(cv, x$seasons), x$observed)
})

test_that("explicit characteristic inputs do not require inverse conversion", {
  model <- pheno_model()
  testthat::local_mocked_bindings(kinetic_to_characteristic = function(...) stop("must not invert"),
    .package = "evalpheno")
  initial <- default_parameters(model, "characteristic")
  expect_identical(parameter_schema(model, "characteristic")$name, names(initial))
  expect_identical(names(default_bounds(model, parameterization = "characteristic")$lower), names(initial))
  for (m in list(model, combined_pheno_model(model, 2), stage_pheno_models(model, c(180, 250)),
                pheno_model_list(list(model, model)))) {
    schema <- parameter_schema(m, "characteristic")
    p <- setNames(schema$default, schema$name)
    expect_true(validate_parameters(m, p))
    fit <- fit_phenology(m, function(p, m) sum(p), lower = numeric(), upper = numeric(), parameters = p)
    expect_true(validate_parameters(fit$model, model_parameters(fit)))
    expect_true(all(vapply(as_pheno_models(m, p), function(x) all(is.finite(x$parameters)), logical(1))))
  }
})

test_that("kinetic fitting preserves matching metadata and clears stale characteristic values", {
  base <- pheno_model(parameters = c(theta_star = 6))
  for (model in list(base, stage_pheno_models(base, c(180, 250)),
                     combined_pheno_model(base, 2), pheno_model_list(list(base, base)))) {
    initial <- model_parameters(model)
    same <- fit_phenology(model, function(p, m) sum(p), numeric(), numeric())
    for (m in as_pheno_models(same)) expect_identical(m$chill$characteristic, base$chill$characteristic)
    initial["E0"] <- initial[["E0"]] + 1
    changed <- fit_phenology(model, function(p, m) sum(p), numeric(), numeric(), parameters = initial)
    for (m in as_pheno_models(changed)) {
      expect_null(m$chill$characteristic)
      expect_identical(m$parameters[["E0"]], initial[["E0"]])
    }
  }
})
