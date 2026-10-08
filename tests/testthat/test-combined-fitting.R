combined_example <- function() {
  base <- pheno_model("sequential", chill_dynamic("kinetic"))
  p <- base$parameters
  p[c("yc", "zc")] <- c(0.1, 5)
  base <- pheno_model("sequential", chill_dynamic("kinetic"), parameters = p)
  model <- combined_pheno_model(base, n_cultivars = 3)
  weather <- data.frame(Temp = rep(8, 480), Year = 2008,
                        JDay = rep(50:69, each = 24), Hour = rep(0:23, 20))
  warmer <- weather
  warmer$Temp <- 10
  seasons <- list(a = list(cold = weather, warm = warmer),
                  b = list(cold = weather), c = list(warm = warmer))
  parameters <- model$parameters
  parameters[paste0("zc", 1:3)] <- c(5, 10, 15)
  list(base = base, model = model, seasons = seasons, parameters = parameters)
}

test_that("combined specifications expand structure traits and preserve stored initial values", {
  for (structure in c("sequential", "parallel", "partial_overlap", "phenoflex")) {
    base <- pheno_model(structure, chill_dynamic("kinetic"))
    base$parameters["yc"] <- 45
    model <- combined_pheno_model(base, 3)
    base_schema <- parameter_schema(base)
    traits <- base_schema$name[base_schema$component == "structure"]
    schema <- parameter_schema(model)
    expect_s3_class(model, "combined_pheno_model")
    expect_true(validate_model_spec(model))
    expect_true(validate_parameters(model, rev(model$parameters)))
    expect_length(model$parameters, length(base$parameters) + 2L * length(traits))
    expect_identical(unname(model$parameters[paste0("yc", 1:3)]), rep(45, 3))
    expect_identical(schema$name, names(model$parameters))
    expect_true(all(is.na(schema$cultivar[schema$component != "structure"])))
    for (i in 1:3) expect_identical(cultivar_parameters(model, i), base$parameters)
    expect_identical(unname(default_parameters(model)[paste0("yc", 1:3)]), rep(40, 3))
  }
  partial <- combined_pheno_model(pheno_model(), 3, cultivar_specific = c("yc", "zc"))
  expect_identical(names(partial$parameters)[1:7], c("yc1", "yc2", "yc3", "zc1", "zc2", "zc3", "s1"))
  shared <- combined_pheno_model(pheno_model(), 3, cultivar_specific = character())
  expect_identical(shared$parameters, shared$model$parameters)
  expect_true(all(is.na(parameter_schema(shared)$cultivar)))
  one <- combined_pheno_model(pheno_model(), 1, cultivar_specific = "yc")
  expect_true("yc1" %in% names(one$parameters))
})

test_that("combined predictions select each cultivar's traits and share submodel parameters", {
  x <- combined_example()
  x$parameters["Tb"] <- 2
  predicted <- predict_combined_phenology(x$model, x$seasons, rev(x$parameters))
  expect_named(predicted, names(x$seasons))
  expect_named(predicted[[1]], names(x$seasons[[1]]))
  expect_identical(lengths(predicted), c(a = 2L, b = 1L, c = 1L))
  for (i in 1:3) {
    p <- x$base$parameters
    p["zc"] <- c(5, 10, 15)[i]
    p["Tb"] <- 2
    expected <- vapply(x$seasons[[i]], function(weather) {
      predict_phenology(x$base, weather, p)
    }, numeric(1))
    expect_equal(predicted[[i]], expected)
  }
  expect_identical(predict_phenology(x$model, x$seasons, x$parameters), predicted)
  detail <- predict_combined_phenology(x$model, x$seasons, x$parameters, basic_output = FALSE)
  expect_named(detail[[1]][[1]], c("bloomindex", "chill", "z"))
  expect_identical(x$model$parameters[paste0("zc", 1:3)], c(zc1 = 5, zc2 = 5, zc3 = 5))
})

test_that("combined RSS sums all matched residuals and excludes missing bloom", {
  x <- combined_example()
  predicted <- predict_combined_phenology(x$model, x$seasons, x$parameters)
  expect_true(all(is.finite(unlist(predicted))))
  observed <- predicted
  observed[[1]] <- observed[[1]] + c(2, -3)
  observed[[2]] <- observed[[2]] + 4
  observed[[3]] <- observed[[3]] - 5
  expect_equal(phenology_rss_combined(x$parameters, x$model, x$seasons, predicted), 0)
  expect_equal(phenology_rss_combined(rev(x$parameters), x$model, x$seasons, observed), 54)
  separate <- vapply(1:3, function(i) {
    phenology_rss(cultivar_parameters(x$model, i, x$parameters), x$base,
                  x$seasons[[i]], observed[[i]])
  }, numeric(1))
  expect_equal(phenology_rss_combined(x$parameters, x$model, x$seasons, observed), sum(separate))
  x$parameters["zc2"] <- 1e12
  expect_identical(phenology_rss_combined(x$parameters, x$model, x$seasons, observed), Inf)
})

test_that("combined models validate selected traits, parameters and nested data", {
  x <- combined_example()
  for (n in list(0, 1.5, Inf, NA_real_, c(2, 3)))
    expect_error(combined_pheno_model(x$base, n), "positive integer")
  expect_error(combined_pheno_model(population_pheno_model(), 3), "single pheno_model")
  for (specific in list("unknown", "Tb", c("yc", "yc"), NA_character_))
    expect_error(combined_pheno_model(x$base, 3, specific), "structure parameter names")
  expect_error(validate_parameters(x$model, x$parameters[-1]), "expanded schema")
  bad <- x$parameters
  bad["yc2"] <- -1
  expect_error(validate_parameters(x$model, bad), "Cultivar 2.*positive")
  bad <- x$parameters
  bad["Tu"] <- bad["Tb"]
  expect_error(validate_parameters(x$model, bad), "Tb < Tu < Tc")
  custom <- combined_pheno_model(x$base, 3, parameters = rev(x$parameters))
  expect_identical(custom$parameters, x$parameters)
  for (i in list(0, 4, NA_real_, 1.5))
    expect_error(cultivar_parameters(x$model, i), "cultivar index")
  expect_error(predict_combined_phenology(x$model, x$seasons[-1]), "seasonlist")
  bad <- x$seasons
  bad[[2]] <- list()
  expect_error(predict_combined_phenology(x$model, bad), "seasonlist")
  observed <- predict_combined_phenology(x$model, x$seasons, x$parameters)
  expect_error(phenology_rss_combined(x$parameters, x$model, x$seasons, observed[-1]), "per cultivar")
  for (observation in list(c(55, 56), NA_real_, "55", matrix(55))) {
    bad <- observed
    bad[[2]] <- observation
    expect_error(phenology_rss_combined(x$parameters, x$model, x$seasons, bad), "Cultivar 2")
  }
})

test_that("both optimizers fit expanded vectors while preserving omitted shared values", {
  x <- combined_example()
  observed <- predict_combined_phenology(x$model, x$seasons, x$parameters)
  # Bounds built with R's repeated-vector naming syntax are accepted directly.
  lower <- c(zc = rep(1, 3))
  upper <- c(zc = rep(20, 3))
  baseline <- phenology_rss_combined(x$model$parameters, x$model, x$seasons, observed)
  for (optimizer in c("DEoptim", "GenSA")) {
    skip_if_not_installed(optimizer)
    control <- if (optimizer == "DEoptim") list(NP = 30, trace = FALSE) else
      list(max.call = 400)
    fit <- fit_phenology(x$model, evaluation = phenology_rss_combined,
                         lower = lower, upper = rev(upper),
                         seasons = x$seasons, observed = observed,
                         optimizer = optimizer, seed = 17, max_iterations = 20,
                         control = control)
    expect_lt(fit$value, baseline)
    expect_equal(fit$value, phenology_rss_combined(fit$par, x$model, x$seasons, observed))
    fixed <- setdiff(names(x$model$parameters), names(lower))
    expect_identical(fit$par[fixed], x$model$parameters[fixed])
    expect_identical(fit$model$parameters, fit$par)
    expect_equal(predict_combined_phenology(fit$model, x$seasons),
                  predict_combined_phenology(x$model, x$seasons, fit$par))
    expect_true(validate_parameters(fit$model, fit$par))
  }
})

test_that("expanded bounds can fix only part of a cultivar-specific parameter group", {
  x <- combined_example()
  observed <- predict_combined_phenology(x$model, x$seasons, x$model$parameters)
  fit <- fit_phenology(x$model, phenology_rss_combined, lower = c(zc2 = 5),
                       upper = c(zc2 = 5), seasons = x$seasons, observed = observed)
  expect_identical(fit$diagnostics$termination, "fixed")
  expect_equal(fit$value, 0)
  expect_identical(fit$model, x$model)
})
