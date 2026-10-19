# Keep package checks small; the full regression suite is in development/tests/.
smoke_weather <- data.frame(
  Temp = rep(8, 480), Year = 2008,
  JDay = rep(50:69, each = 24), Hour = rep(0:23, 20)
)
smoke_model <- function() {
  pheno_model("sequential", chill_dynamic("kinetic"),
              parameters = c(yc = 0.1, zc = 5))
}

test_that("native chill, heat and coupling kernels retain known values", {
  A1 <- exp(2 / 281)
  chill <- calculate_chill_dynamic(rep(8, 5), 0:4,
    A0 = 10 * A1 * exp(-1 / 281), A1 = A1, E0 = 1, E1 = 2, Tf = 8)
  expect_equal(unname(chill[2, ]), c(5 * (1 - exp(-1)), 10, 5 * (1 - exp(-1))),
               tolerance = 1e-12)
  expect_equal(calculate_heat_gdh(c(4, 25, 36, 0), 0:3), c(0, 21, 0))
  out <- apply_sequential_structure(matrix(c(0, 0.4, 1, 1.2), ncol = 1),
    c(2, 3, 5), yc = 1, zc = 4, stopatzc = FALSE, basic_output = FALSE)
  expect_equal(out$z, c(0, 0, 3, 8))
  expect_equal(out$bloomindex, 4)
  expect_equal(return_JDay(24, rep(60, 24), rep(2024, 24)), 60.5)
})

test_that("Dynamic Model parameter conversion round trips", {
  characteristic <- c(40, 190, .5, 25, 279, 286.1, 47.7, 28, 4, 36, 4, 1.6)
  kinetic <- characteristic_to_kinetic(characteristic)
  expect_equal(kinetic_to_characteristic(kinetic), characteristic, tolerance = 1e-6)
  expect_identical(kinetic[-(5:8)], characteristic[-(5:8)])
})

test_that("all four structures construct and predict a short season", {
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    requirements <- if (structure == "partial_overlap") c(b1 = 5, b2 = 10) else c(zc = 5)
    model <- pheno_model(structure, chill_dynamic("kinetic"),
                          parameters = c(yc = 0.1, requirements))
    expect_true(validate_parameters(model, model_parameters(model)))
    expect_identical(names(model_parameters(model)), parameter_schema(model)$name)
    predicted <- predict_phenology(model, smoke_weather)
    expect_true(is.finite(predicted), info = structure)
    details <- predict_phenology(model, smoke_weather, basic_output = FALSE)
    expect_equal(predicted, return_JDay(details$bloomindex, smoke_weather$JDay,
                                       smoke_weather$Year))
  }
  model <- smoke_model()
  expect_identical(predict_phenology(model, list(a = smoke_weather)),
                   c(a = predict_phenology(model, smoke_weather)))
  expect_true(is.na(predict_phenology(model, smoke_weather,
    parameters = replace(model$parameters, "yc", 1e12))))
})

test_that("invalid model parameters and weather fail clearly", {
  model <- smoke_model()
  expect_error(pheno_model(parameters = c(yc = -1)), "positive")
  expect_error(validate_parameters(model, c(yc = 1)), "schema names")
  bad <- smoke_weather
  bad$Temp[1] <- NA_real_
  expect_error(predict_phenology(model, bad), "finite numeric")
  expect_error(predict_phenology(model, smoke_weather[-1, ]), "complete hourly")
})

test_that("RSS and fixed fits use the supplied observations", {
  model <- smoke_model()
  predicted <- predict_phenology(model, smoke_weather)
  seasons <- list(smoke_weather)
  expect_equal(phenology_rss(model$parameters, model, seasons, predicted + 2), 4)
  fit <- fit_phenology(model, phenology_rss, lower = numeric(), upper = numeric(),
                       seasons = seasons, observed = predicted)
  expect_s3_class(fit, "phenology_fit")
  expect_equal(fit$value, 0)
  expect_identical(fit$diagnostics$termination, "fixed")
  expect_identical(predict_phenology(fit, smoke_weather), predicted)
})

test_that("optimizer adapters run with a small evaluation budget", {
  model <- smoke_model()
  loss <- function(parameters, model) (parameters[["yc"]] - 0.15)^2
  for (optimizer in c("DEoptim", "GenSA")) {
    if (!requireNamespace(optimizer, quietly = TRUE)) next
    control <- if (optimizer == "DEoptim") list(NP = 10, trace = FALSE) else
      list(max.call = 20, simple.function = TRUE)
    fit <- fit_phenology(model, loss, lower = c(yc = 0.05), upper = c(yc = 0.2),
      optimizer = optimizer, seed = 17, max_iterations = 1, control = control)
    expect_equal(fit$value, loss(fit$par, model))
    expect_lte(fit$value, loss(model$parameters, model))
    expect_identical(fit$par[-1], model$parameters[-1])
    expect_identical(fit$model$parameters, fit$par)
  }
})

test_that("collections and constant populations retain single-model predictions", {
  model <- smoke_model()
  predicted <- predict_phenology(model, smoke_weather)
  combined <- combined_pheno_model(model, 2)
  expect_equal(unlist(predict_phenology(combined, smoke_weather), use.names = FALSE),
               rep(predicted, 2))
  stages <- stage_pheno_models(model, c(early = 5, late = 10))
  dates <- unlist(predict_phenology(stages, smoke_weather), use.names = FALSE)
  expect_true(all(is.finite(dates)))
  expect_equal(dates[1], predicted)
  expect_gte(dates[2], dates[1])
  population <- population_pheno_model("sequential", chill_dynamic("kinetic"), n = 2)
  expect_equal(predict_phenology(population, smoke_weather, parameters = model$parameters)$bloom_jday,
               rep(predicted, 2))
})

test_that("two-fold fixed cross-validation and ensembles predict", {
  model <- smoke_model()
  seasons <- lapply(2008:2011, function(year) {
    weather <- smoke_weather
    weather$Year <- year
    weather
  })
  observed <- predict_phenology(model, seasons)
  prepared <- prepare_phenology_cv(seasons, observed, v = 2, seed = 17)
  cv <- fit_phenology_cv(model, prepared, phenology_rss,
                         lower = numeric(), upper = numeric())
  expect_s3_class(cv, "phenology_cv")
  expect_equal(cv$calibration$scores$mse, c(0, 0))
  expect_equal(predict_phenology(pheno_ensemble(cv), smoke_weather), observed[1])
})
