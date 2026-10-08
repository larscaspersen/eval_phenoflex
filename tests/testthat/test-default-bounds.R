test_that("default bounds cover each schema and distinguish heat and chill choices", {
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    for (representation in c("characteristic", "kinetic")) {
      for (scaling in c("unscaled", "scaled")) {
        model <- pheno_model(structure, chill_dynamic(representation), heat_gdh(scaling))
        before <- model
        bounds <- default_bounds(model)
        expect_named(bounds, c("lower", "upper"))
        expect_identical(names(bounds$lower), parameter_schema(model)$name)
        expect_identical(names(bounds$upper), names(bounds$lower))
        expect_true(all(is.finite(bounds$lower) & is.finite(bounds$upper)))
        expect_true(all(bounds$lower < bounds$upper))
        expect_true(all(model$parameters >= bounds$lower & model$parameters <= bounds$upper))
        expect_identical(model, before)
      }
      unscaled <- default_bounds(pheno_model(structure, chill_dynamic(representation), heat_gdh("unscaled")))
      scaled <- default_bounds(pheno_model(structure, chill_dynamic(representation), heat_gdh("scaled")))
      requirements <- intersect(c("zc", "b1", "b2"), names(unscaled$lower))
      other <- setdiff(names(unscaled$lower), requirements)
      for (side in c("lower", "upper")) {
        expect_identical(scaled[[side]][requirements], 21 * unscaled[[side]][requirements])
        expect_identical(scaled[[side]][other], unscaled[[side]][other])
      }
    }
  }
  model <- pheno_model()
  bounds <- default_bounds(model)
  expect_identical(bounds$lower[c("yc", "zc", "s1", "Tb", "Tu", "Tc", "Tf", "slope")],
                   c(yc = 10, zc = 150, s1 = 0.1, Tb = 0, Tu = 10, Tc = 20, Tf = 0, slope = 0.1))
  expect_identical(bounds$upper[c("yc", "zc", "s1", "Tb", "Tu", "Tc", "Tf", "slope")],
                   c(yc = 80, zc = 400, s1 = 1.5, Tb = 10, Tu = 30, Tc = 40, Tf = 10, slope = 5))
  expect_identical(default_bounds(model, c("zc", "yc")),
                   lapply(bounds, function(p) p[c("yc", "zc")]))
  expect_identical(default_bounds(model, character()), list(lower = numeric(), upper = numeric()))
  expect_identical(default_bounds(pheno_model(heat = heat_gdh("scaled")), "zc", heat_scale = 25),
                   list(lower = c(zc = 3750), upper = c(zc = 10000)))
})

test_that("collection bounds repeat selected traits and retain shared names", {
  base <- pheno_model()
  combined <- combined_pheno_model(base, 3, cultivar_specific = c("yc", "zc"))
  bounds <- default_bounds(combined)
  expect_identical(names(bounds$lower), names(model_parameters(combined)))
  expect_identical(bounds$lower[1:7], c(yc1 = 10, yc2 = 10, yc3 = 10,
                                       zc1 = 150, zc2 = 150, zc3 = 150, s1 = 0.1))
  expect_identical(default_bounds(combined, "zc")$upper, c(zc1 = 400, zc2 = 400, zc3 = 400))
  expect_identical(default_bounds(combined, "zc2")$upper, c(zc2 = 400))
  models <- stage_pheno_models(base, c(180, 250))
  expect_identical(default_bounds(models, "zc")$lower, c(zc1 = 150, zc2 = 150))
  simple <- as_pheno_models(models)
  expect_identical(default_bounds(simple), default_bounds(pheno_model_list(simple)))
  expect_identical(default_bounds(population_pheno_model()), default_bounds(base))
})

test_that("bounds selections and reference scaling fail clearly when malformed", {
  model <- pheno_model()
  for (names in list(c("yc", "yc"), NA_character_, "unknown", 1, matrix("yc")))
    expect_error(default_bounds(model, names), "parameter_names")
  for (scale in list(0, -1, NA_real_, Inf, numeric(), c(20, 21), matrix(21)))
    expect_error(default_bounds(model, heat_scale = scale), "heat_scale")
  expect_error(default_bounds(pheno_model(heat = heat_gdh("scaled")), heat_scale = .Machine$double.xmax),
                "overflowed")
})

test_that("fitter defaults use the full schema or the explicitly supplied subset", {
  skip_if_not_installed("GenSA")
  received <- NULL
  testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
    received <<- list(lower = lower, upper = upper)
    list(par = par, value = fn(par))
  }, .package = "GenSA")
  model <- pheno_model()
  loss <- function(p, m) (p[["yc"]] - 45)^2
  fit <- fit_phenology(model, loss, optimizer = "GenSA", max_iterations = 1)
  expect_identical(received, default_bounds(model))
  expect_identical(fit$calibration_settings$lower, default_bounds(model)$lower)
  fit <- fit_phenology(model, loss, upper = c(yc = 70), optimizer = "GenSA")
  expect_identical(received, list(lower = c(yc = 10), upper = c(yc = 70)))
  fixed <- setdiff(names(model$parameters), "yc")
  expect_identical(fit$calibration_settings$lower[fixed], model$parameters[fixed])
  expect_identical(fit$calibration_settings$upper[fixed], model$parameters[fixed])
  fit_phenology(model, loss, lower = c(yc = 20), optimizer = "GenSA")
  expect_identical(received, list(lower = c(yc = 20), upper = c(yc = 80)))
  expect_identical(fit_phenology(model, loss, lower = numeric())$diagnostics$termination, "fixed")
  expect_identical(fit_phenology(model, loss, upper = numeric())$diagnostics$termination, "fixed")
  outside <- pheno_model(parameters = c(zc = 700))
  expect_error(fit_phenology(outside, loss), "Initial parameters")
  expect_identical(outside$parameters[["zc"]], 700)
  models <- list(model, model)
  fit <- fit_phenology(models, function(p, m) (p[["yc1"]] - 45)^2, optimizer = "GenSA")
  expect_identical(fit$calibration_settings$lower, default_bounds(models)$lower)
})

test_that("numerical conversion failures are excluded without hiding evaluation bugs", {
  skip_if_not_installed("GenSA")
  model <- pheno_model()
  bad <- model$parameters
  bad[c("theta_star", "theta_c", "tau", "pie_c")] <- c(279, 286.5, 16, 50)
  weather <- data.frame(Temp = rep(8, 480), Year = 2008, JDay = rep(50:69, each = 24))
  expect_true(validate_parameters(model, bad))
  expect_error(predict_phenology(model, weather, bad), class = "phenology_conversion_error")
  testthat::local_mocked_bindings(GenSA = function(par, fn, lower, upper, control) {
    expect_identical(fn(bad[names(par)]), Inf)
    list(par = par, value = fn(par))
  }, .package = "GenSA")
  fit <- fit_phenology(model, function(p, m) {
    predict_phenology(m, weather, p)
    (p[["yc"]] - 40)^2
  }, optimizer = "GenSA")
  expect_identical(fit$value, 0)
  expect_identical(fit$diagnostics$evaluations, 3L)
  expect_error(fit_phenology(model, function(p, m) stop("broken evaluator"), optimizer = "GenSA"),
                "broken evaluator")
})
